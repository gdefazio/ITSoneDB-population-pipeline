#!/usr/bin/env python3
__author__ = 'Bruno Fosso'
__version__ = "1.0"

from librs.return_time import return_time
import argparse
import os
import shlex
import subprocess
import multiprocessing as mp
import sys
from random import choice
import argcomplete
from Bio import SeqIO
from gzip import open as gzopen


def input_options():
    parser = argparse.ArgumentParser(description="split fasta and tsv files", prefix_chars="-")
    parser.add_argument("-f", "--its1_fasta", type=str,
                        help="fasta file containing all ITS1 sequences",
                        action="store", required=True)
    parser.add_argument("-n", "--ncbi_taxonomy", type=str,
                        help="folder containing NBCI taxonomy dump files",
                        action="store", required=True)
    parser.add_argument("-o", "--output_folder", type=str, help="output folder",
                        action="store", required=False,
                        default="./etl")
    parser.add_argument("-p", "--processes", type=int, help="processes (min num 2)",
                        action="store", required=False,
                        default=2)
    argcomplete.autocomplete(parser)
    return parser.parse_args()


def generazione_dizionari_tassonomia(taxnomy_folder):
    if os.path.exists(os.path.join(taxnomy_folder, "nodes.dmp")):
        nodesfile = open(os.path.join(taxnomy_folder, "nodes.dmp"))
    else:
        sys.exit("No NODESFILE")

    if os.path.exists(
            os.path.join(taxnomy_folder, "names.dmp")):
        namesfile = open(
            os.path.join(taxnomy_folder, "names.dmp"))
    else:
        sys.exit("no namesfile")

    if os.path.exists(os.path.join(taxnomy_folder, "merged.dmp")):
        merged = open(os.path.join(taxnomy_folder, "merged.dmp"))
    else:
        sys.exit("no merged")

    if os.path.exists(os.path.join(taxnomy_folder, "delnodes.dmp")):
        deleted = open(os.path.join(taxnomy_folder, "delnodes.dmp"))
    else:
        sys.exit("no deleted")

    return_time("Generazione Dizionari")
    node2name = {}
    name2node = {}
    for linea in namesfile:
        linea = linea.strip()
        fields = list(map(str.strip, linea.split("|")))
        if fields[3] == "scientific name":
            nodeid, name = fields[0], fields[1]
            node2name[nodeid] = name
            name2node[name] = nodeid
    node2name["53267"] = "Trebouxia jamesii"

    node2parent = {}
    node2order = {}
    for linea in nodesfile:
        linea = linea.strip()
        fields = list(map(str.strip, linea.split("|")))
        nodeid, parentid, ordine = fields[0], fields[1], fields[2]
        node2parent[nodeid] = parentid
        node2order[nodeid] = ordine
    node2parent["53267"] = "13786"
    node2order["53267"] = "species"

    node2merged = {}
    for linea in merged:
        linea = linea.strip()
        fields = list(map(str.strip, linea.split("|")))
        nodo, new_id = fields[0], fields[1]
        node2merged[nodo] = new_id

    delnodes = []
    for linea in deleted:
        linea = linea.strip()
        fields = list(map(str.strip, linea.split("|")))
        delnodes.append(fields[0])
    return node2parent, node2order, node2name, node2merged, delnodes


def accession(identifier):
    parts = identifier.split(".")
    return parts[0]


def species2accession_lista(its1_seq, node2parent, node2merged, delnodes, node2order):
    """

    :param node2order: dictionary
    :param delnodes: list
    :param node2merged: dictionary
    :param node2parent: dictionary
    :param its1_seq: str
    :return species2acc: dictionary
    """
    acc2taxa = {}
    with gzopen(its1_seq, 'rt') as aa:
        for linea in aa:
            if linea.startswith(">"):
                s = list(map(str.strip, linea.split("|")))
                acc_num, node = s[0].lstrip(">"), s[1]
                # return_time(acc_num)
                acc2taxa[acc_num] = node
                # if node in node2parent.keys():
                #     acc2taxa[acc_num] = node
                # elif node in node2merged.keys():
                #     acc2taxa[acc_num] = node2merged[node]
                #     # node2order[acc] = "GB acc"
                # elif node in delnodes:
                #     print acc_num, node, "deleted"
                # else:
                #     print acc_num, node, "nothing"
    for acc_num, node in acc2taxa.items():
        if node2order[node] != "species":
            parent = node2parent[node]
            while node != parent:
                if node2order[parent] == "species":
                    acc2taxa[acc_num] = parent
                    node = parent
                else:
                    node = parent
                    parent = node2parent[node]
    species2acc = {}
    for acc_num in acc2taxa.keys():
        # print acc
        specie = acc2taxa[acc_num]
        species2acc.setdefault(specie, set())
        species2acc[specie].add(acc_num)
    return species2acc


# noinspection PyUnboundLocalVariable
def acc2seq(fasta_collection):
    """

    :rtype: dict
    """
    acctoseq = {}
    with gzopen(fasta_collection, 'rt') as aa:
        for linea in aa:
            if linea.startswith(">"):
                acc_num = list(map(str.strip, linea.split("|")))[0].lstrip(">")
                acctoseq.setdefault(acc_num, "")
                acctoseq[acc_num] += linea
            else:
                acctoseq[acc_num] += linea
    return acctoseq


def clustering_utility(fasta):
    cluster_file = os.path.join(species_collection_subfolder,
                                fasta.rstrip(".fa.gz") + ".usearch_clustering")
    cmd = shlex.split("vsearch --cluster_fast %s -id 0.97 -uc  %s" % (
        os.path.join(species_collection_subfolder, fasta), cluster_file))
    p = subprocess.Popen(cmd,
                         stdout=subprocess.DEVNULL,
                         stderr=subprocess.STDOUT)
    p.wait()
    with open(cluster_file, 'rt') as a:
        cluster_data = {}
        seq2len = {}
        for line in a:
            field = list(map(str.strip, line.split("\t")))
            if field[0] == "S":
                seq2len[field[-2]] = int(field[2])
                cluster_data.setdefault(field[-2], set())
                cluster_data[field[-2]].add(field[-2])
            elif field[0] == "H":
                seq2len[field[-2]] = int(field[2])
                cluster_data.setdefault(field[-1], set())
                cluster_data[field[-1]].add(field[-2])
    if len(cluster_data) == 0:
        return 'manual_control', fasta
    else:
        size = []
        for cluster in cluster_data.values():
            size.append(len(cluster))
        greatest = max(size)
        if size.count(greatest) == 1:
            len_max = []
            for cluster in cluster_data.keys():
                if len(cluster_data[cluster]) == greatest:
                    return "%s\t%s\n" % (fasta.split(".")[0], cluster)
        else:
            candidate = []
            for cluster in cluster_data.keys():
                if len(cluster_data[cluster]) == greatest:
                    candidate.append(cluster)
            rep = choice(candidate)

            return "%s\t%s\n" % (fasta.split(".")[0], rep)#.split("|")[0])


if __name__ == "__main__":
    param = input_options()
    its1_fa, output_folder, ncbi_folder, processes = param.its1_fasta, param.output_folder, param.ncbi_taxonomy, param.processes
    processes = processes//2
    output_folder = os.path.join(output_folder, 'Representative_computation')
    accession2sequence = acc2seq(its1_fa)
    parent_dict, order_dict, name_dict, merged_dict, deleted_nodes = generazione_dizionari_tassonomia(ncbi_folder)
    species2seq = species2accession_lista(its1_fa, parent_dict, merged_dict, deleted_nodes, order_dict)
    if os.path.exists(output_folder) is False:
        os.mkdir(output_folder)
    species_collection_subfolder = os.path.join(output_folder, "seq_to_species")
    if os.path.exists(species_collection_subfolder) is False:
        os.mkdir(species_collection_subfolder)

    scritte = 0
    for species, acc_list in species2seq.items():
        with gzopen(os.path.join(species_collection_subfolder, species + ".fa.gz"), "wt") as fasta:
            for acc in acc_list:

                fasta.write(accession2sequence[acc])
                scritte += 1
        if os.stat(os.path.join(species_collection_subfolder, "%s.fa.gz" % species))[6] == 0:
            os.remove(os.path.join(species_collection_subfolder, "%s.fa.gz" % species))
            return_time(os.path.join(species_collection_subfolder,"%s.fa.gz is empty" % species))

    return_time("%i sequences processed" % scritte)
    rep_file = open(os.path.join(output_folder, "rep_sequences.csv"), "wt")
    manual_control = set()
    to_cluster = list()
    for fasta in os.listdir(species_collection_subfolder):
        if fasta.endswith("fa.gz") and os.path.exists(os.path.join(species_collection_subfolder, fasta)):
            with gzopen(os.path.join(species_collection_subfolder, fasta), 'rt') as handle:
                parsed_fasta = list(SeqIO.parse(handle, "fasta"))
                if len(parsed_fasta) == 1:
                    # record = SeqIO.read(handle, "fasta")
                    if len(parsed_fasta[0].name) != 0:
                        rep_file.write("%s\t%s\n" % (fasta.split(".")[0], parsed_fasta[0].name))
                    else:
                        raise Exception("ParsingFastaError of %s" % fasta)
                else:
                    to_cluster.append(fasta)
        else:
            return_time("%s not in %s" % (fasta, species_collection_subfolder))
    # Parallelized clustering procedure added by gdefazio in 24/03/2021
    if len(to_cluster) > 0:
        return_time("Start clustering procedure for %s species" % len(to_cluster))
        if processes == 1:
            for el in to_cluster:
                clustering_utility(el)
        else:
            with mp.Pool(processes=processes) as p:
                results = p.map(
                                func=clustering_utility,
                                iterable=to_cluster,
                                chunksize=len(to_cluster)//processes
                                )
    else:
        return_time("There are not species sequences to cluster")
                # cluster_file = os.path.join(species_collection_subfolder,
                #                             fasta.rstrip(".fasta") + ".usearch_clustering")
                # cmd = shlex.split("vsearch --cluster_fast %s -id 0.97 -uc  %s" % (
                #     os.path.join(species_collection_subfolder, fasta), cluster_file))
                # p = subprocess.Popen(cmd)
                # p.wait()
                # with open(cluster_file, 'rt') as a:
                #     cluster_data = {}
                #     seq2len = {}
                #     for line in a:
                #         field = list(map(str.strip, line.split("\t")))
                #         if field[0] == "S":
                #             seq2len[field[-2]] = int(field[2])
                #             cluster_data.setdefault(field[-2], set())
                #             cluster_data[field[-2]].add(field[-2])
                #         elif field[0] == "H":
                #             seq2len[field[-2]] = int(field[2])
                #             cluster_data.setdefault(field[-1], set())
                #             cluster_data[field[-1]].add(field[-2])
                # if len(cluster_data) == 0:
                #     manual_control.add(fasta)
                # else:
                #     size = []
                #     for cluster in cluster_data.values():
                #         size.append(len(cluster))
                #     greatest = max(size)
                #     if size.count(greatest) == 1:
                #         len_max = []
                #         for cluster in cluster_data.keys():
                #             if len(cluster_data[cluster]) == greatest:
                #                 rep_file.write("%s\t%s\n" % (fasta.split(".")[0], cluster))
                #     else:
                #         candidate = []
                #         for cluster in cluster_data.keys():
                #             if len(cluster_data[cluster]) == greatest:
                #                 candidate.append(cluster)
                #         rep = choice(candidate)
                #         rep_file.write("%s\t%s\n" % (fasta.split(".")[0], rep.split("|")[0]))
    for el in results:
        if type(el) is str:
            rep_file.write(el)
        else:
            manual_control.add(el[1])
    rep_file.close()
    tmp = open(os.path.join(output_folder, "manual_controll_rep_species"), "wt")
    tmp.write("%s" % "\n".join(list(manual_control)))
    tmp.close()
    return_time('DONE')
