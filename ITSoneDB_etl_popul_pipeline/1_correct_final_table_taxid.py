#!/usr/bin/env python3
import os
import sys
import argparse
import argcomplete
from gzip import open as gzopen
from librs.return_time import return_time


def input_options():
    parser = argparse.ArgumentParser(
        description="This script corrects the final table in order to obtain the updated NCBI taxonomy ID",
        prefix_chars="-")
    parser.add_argument("-i", "--its_file", type=str,
                        help="tsv containing ITS1 location per accession",
                        action="store", required=True)
    parser.add_argument("-o", "--output_folder", type=str, help="output folder", action="store", required=False,
                        default="etl")
    parser.add_argument("-t", "--tax_folder", type=str, help="taxonomy_folder",
                        action="store", required=True)
    parser.add_argument("-m", "--mito_check", type=str,
                        help="log of ENA parsing containing also the info relative to organelle checking",
                        action="store", required=True)
    parser.add_argument("-p", "--previous_taxonomy", type=str,
                        help="path to the taxonomy data of the previous release",
                        action="store", required=True)
    argcomplete.autocomplete(parser)
    return parser.parse_args()


# noinspection DuplicatedCode
def generazione_dizionari_tassonomia(taxnomy_folder):
    if os.path.exists(os.path.join(taxnomy_folder, "nodes.dmp")):
        nodesfile = open(os.path.join(taxnomy_folder, "nodes.dmp"), 'rt')
    else:
        sys.exit("No NODESFILE")

    if os.path.exists(os.path.join(taxnomy_folder, "merged.dmp")):
        merged = open(os.path.join(taxnomy_folder, "merged.dmp"), 'rt')
    else:
        sys.exit("no merged")

    if os.path.exists(os.path.join(taxnomy_folder, "delnodes.dmp")):
        deleted = open(os.path.join(taxnomy_folder, "delnodes.dmp"), 'rt')
    else:
        sys.exit("no deleted")

    node2parent = {}
    node2order = {}
    for linea in nodesfile:
        linea = linea.strip()
        fields = list(map(str.strip, linea.split("|")))
        nodeid, parentid, ordine = fields[0], fields[1], fields[2]
        node2parent[nodeid] = parentid
        node2order[nodeid] = ordine

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
    if os.path.exists(os.path.join(taxnomy_folder, 'nucl_gb.accession2taxid.gz')):
        acc2taxid = {}
        with gzopen(os.path.join(taxnomy_folder, 'nucl_gb.accession2taxid.gz'), 'rt') as gz:
            gz.readline()
            for l in gz:
                s = l.split('\t')
                acc2taxid[s[0]] = s[2]
        return node2parent, node2merged, delnodes, acc2taxid
    else:
        return node2parent, node2merged, delnodes


def verifica_id(taxid, node2parent, merged, deleted):
    """
    This function verify the NCBI taxonomy ID associated to each sequences
    """
    if node2parent.get(taxid):
        pass
    elif merged.get(taxid):
        taxid = merged[taxid]
    elif taxid in deleted:
        taxid = None
    else:
        taxid = None
    path = [taxid]
    nodeid = taxid
    if taxid is None:
        taxid = "TO VERIFY"
        path = []
    else:
        parent = node2parent[taxid]
        while nodeid != parent:
            path.append(parent)
            nodeid = parent
            parent = node2parent[nodeid]
    return taxid, path


def mito_acc(log_data):
    """
    its1_ena.log
    """
    acc_list = []
    if os.path.exists(log_data):
        with open(log_data, 'rt') as l:
            for linea in l:
                if linea.find("organelle") >= 0:
                    s_i = list(map(str.strip, linea.split(" ")))
                    acc_list.append(s_i[-4])
    else:
        sys.exit("The %s doesn't exist" % log_data)
    return acc_list


def remap_deleted_nodes(taxid, node2parent, merged, deleted, ii):
    jj = 0
    if taxid in node2parent:
        while jj < ii:
            taxid = node2parent[taxid]
            jj += 1
    elif taxid in merged:
        taxid = merged[taxid]
        while jj < ii:
            taxid = node2parent[taxid]
            jj += 1
    else:
        taxid = "2759"
    return taxid


if __name__ == "__main__":
    return_time("Start taxonomy correction and mito/plast check")
    param = input_options()
    its1_tab, tax_folder, out, mito, old_tax = \
        param.its_file, param.tax_folder, param.output_folder, param.mito_check, param.previous_taxonomy
    correct = "%s_corrected.tsv" % os.path.split(its1_tab)[1].split(".")[0]
    parent_dict, merged_dict, deleted_list, acc2taxid = generazione_dizionari_tassonomia(tax_folder)
    old_parent_dict, old_merged_dict, old_deleted_list = generazione_dizionari_tassonomia(old_tax)
    organelle_acc = mito_acc(mito)
    # print organelle_acc
    with open(os.path.join(out, correct), "wt") as tmp, open(its1_tab, 'rt') as a:
        line = a.readline()
        tmp.write(line)
        for line in a:
            field = list(map(str.strip, line.split("\t")))
            # # print(field[0])
            node, tax_pat = verifica_id(field[2], parent_dict, merged_dict, deleted_list)
            if node == 'TO VERIFY':
                # verifica che il problema sia dovuto ad una
                # mancato aggiornamento del taxid da parte di ENA
                # rispetto a quello riportato da NCBI
                try:
                    node = acc2taxid[field[0]]
                    node, tax_pat = verifica_id(node, parent_dict, merged_dict, deleted_list)
                except KeyError:
                    pass
            i = 1
            # qui ho aggiunto una verifica che rende il codice un po' verboso ma ci permette di controllare la
            # cancellazionde dei nodi della tassonomia NCBI.
            # in pratica si usa la cartella della tassonomia relativa alla precedente release e su questa si
            # procede al check del nodo.
            # Se anche nella vecchia tassonomia c'e' nulla il nodo diventa None e la entry esclusa.
            while node == "TO VERIFY":
                return_time("to verify loop active for %s" % field[0])
                node = field[2]
                # print(node, i)
                node = remap_deleted_nodes(node, old_parent_dict, old_merged_dict, old_deleted_list, i)
                # print(node)
                node, tax_pat = verifica_id(node, parent_dict, merged_dict, deleted_list)
                # print(node)
                i += 1

            if field[0] not in organelle_acc and node is not None:
                field[2] = node
                tmp.write("%s\n" % "\t".join(field))
            else:
                if field[6] == "True":
                    # qui verifichiamo che abbiamo disponibile il macht per l'HMM
                    # nel qualcaso consideriamo un errore di annotazione
                    # l'organello
                    field[2] = node
                    field[-2] = field[
                                    -2] + ";" + "The ENA annotation if any has been discarded considering it is annotated as a non-nuclear sequence (organelle)"
                    tmp.write("%s\n" % "\t".join(field))
    return_time('DONE')
