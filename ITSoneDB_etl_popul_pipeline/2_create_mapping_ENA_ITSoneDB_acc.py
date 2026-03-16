#!/usr/bin/env python3
__author__ = 'Bruno Fosso'
__version__ = "1.0"
# creation date 28/02/2019
import os
import argparse
import argcomplete
from Bio import SeqIO
import sys
from gzip import open as gzopen
from librs.return_time import return_time

description_info = """This script maps the ENA accession to ITSoneDB accession.\n
The first version was developed during the update to the ENA release 138 and in order to accomplish the ELIXIR IP requirements.\n
It requires the file generate at the end of the ITSoneDB procedure for annotation extraction and HMM mapping.\n
This first file was created by Giuseppe Defazio and it is called "ITSoneDB_entries_generator_definitive.csv".\n
The generated mapping has the following structure:\n
\tENA_ACC: ENA entry accession number\n
\tENA_Version: ENA entry version. This version is inherited by the ITSoneDB entry\n
\tENA_Annotation_ITS1DB: ITSoneDB accession number\n
\tENA_start: ITS1 starting nucleotide on the sequence\n
\tENA_end: ITS1 ending nucleotide on the sequence\n
\tENA_Note: Notes associated to the sequence\n
\tHMM_Annotation_ITS1DB: ITSoneDB accession number\n 
\tHMM_start: ITS1 starting nucleotide on the sequence\n
\tHMM_end: ITS1 ending nucleotide on the sequence\n
\tHMM_Note: Notes associated to the sequence\n
If a previous version of the mapping file is available it can be used as input to correctly update the data\n
Prima di procedere a questa operazione va scaricata la tassonomia NCBI da :
ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
"""


def input_options():
    parser = argparse.ArgumentParser(description=description_info, prefix_chars="-")
    parser.add_argument("-i", "--its_file", type=str,
                        help="tsv containing ITS1 location per accession",
                        action="store", required=True)
    parser.add_argument("-o", "--output_folder", type=str, help="output folder", action="store", required=False,
                        default="etl")
    parser.add_argument("-r", "--ENA_release_version", type=str, help="ENA release version (e.g. -r 138)",
                        action="store", required=True)
    parser.add_argument("-f", "--fasta_index", type=str,
                        help="file containing the correspondence between ENA accession and fasta files",
                        action="store", required=True)
    parser.add_argument("-n", "--ncbi_taxonomy", type=str,
                        help="folder containing NBCI taxonomy dump files",
                        action="store", required=True)
    parser.add_argument("-m", "--previous_mapping_file", type=str, help="mapping file of the prevoius release",
                        action="store", required=False, default=None)
    argcomplete.autocomplete(parser)
    return parser.parse_args()


def starter_value(old_mapping):
    """
    Parse the mapping file of the previous release in order to annotate the data for the new one.
    Considering the release 138 the file is mapping_data_ENA_release_138.tsv in etl folder.
    :param old_mapping: str
    :return starter: int
    :return diz_info: dict
    """
    count_list = []
    diz_info = {}
    if old_mapping is not None:
        with open(old_mapping) as a:
            a.readline()
            for line in a:
                s = list(map(str.strip, line.split("\t")))
                if s[2] != "Null":
                    count_list.append(int(s[2].lstrip("ITS1DB")))
                if s[6] != "Null":
                    count_list.append(int(s[6].lstrip("ITS1DB")))
                diz_info.setdefault(s[0], {})
                diz_info[s[0]]["version"] = s[1]
                diz_info[s[0]]["ENA"] = s[2:6]
                diz_info[s[0]]["HMM"] = s[6:]
    starter = max(count_list) + 1
    return starter, diz_info


def acc2fasta_file(fa_index_file):
    """
    This function acquires the fa_index_file (tsv containing the association between accessions and fasta files) and
    converts it in a dictionary.
    :param fa_index_file: str
    :return acc2file: dict
    """
    acc2file = {}
    with gzopen(fa_index_file, 'rt') as a:
        a.readline()
        for line in a:
            s = list(map(str.strip, line.split("\t")))
            acc2file[s[0]] = s[1]
    return acc2file


def its1_seq(fa_dict, ena_acc, loc, note):
    """

    :param fa_dict: dict
    :param ena_acc: str
    :param loc: list
    :param note: str
    :return its1seq: str
    :return its1seq_flank: str
    """
    with gzopen(os.path.join(out.split('/etl')[0], fa_dict[ena_acc]), 'rt') as a:
        record = SeqIO.read(a, "fasta")
    loc = [int(i.replace("<", "").replace(">", "")) for i in loc]
    its1seq = record.seq[min(loc):max(loc)]
    st = min(loc) - 150
    en = max(loc) + 150
    if st < 0:
        st = 0
    if en > len(record):
        en = len(record)
    its1seq_flank = record.seq[st:en]
    if len({"COMPLEMENT", "RVS"}.intersection(set(note))) >= 1:
        its1seq.reverse_complement()
        its1seq_flank.reverse_complement()
    return str(its1seq), str(its1seq_flank)


def comparing_new_and_old_data(new_data_list, old_data_list):
    result = 0
    for n, o in zip(new_data_list, old_data_list):
        # print(n,o)
        if n != o:
            result += 1
    # print(result)
    return result


# noinspection PyPep8Naming
def create_accession(its1_tsv, old_mapping_file, fasta_diz):
    """
    This function creates the mapping file that will constitute the base for the database table generation.
    :param its1_tsv: str
    :param old_mapping_file: str
    :param fasta_diz: dict
    :return ena2ITSoneDB_accession: dict
    :return fasta_data: dict
    :return its1_db_data: dict
    """
    old_acc, new_acc = 0, 0
    ena2ITSoneDB_accession = []
    fasta_data = {}
    its1_db_data = {}
    starter, previous_info = starter_value(old_mapping_file)
    with open(its1_tsv, 'rt') as a:
        a.readline()
        for line in a:
            s = list(map(str.strip, line.split("\t")))
            acc_num, version, taxid, ena, hmm = s[0], s[1], s[2], s[3], s[6]
            lista_dat = [acc_num, version]
            hmm_note = []
            ena_note = []
            for note in list(map(str.strip, s[9].split("|"))):
                if note in ["MULTI", "COMPLEMENT", "JOIN"]:
                    ena_note.append(note)
                elif note in ["FWD", "RVS"]:
                    hmm_note.append(note)
                elif note == "The ENA annotation if any has been discarded considering it is annotated as a non-nuclear sequence (organelle)":
                    ena_note.append(note)
            if ena == "True":
                if acc_num in previous_info:
                    new_data = [acc_num, version, s[4], s[5]]
                    old_data = [acc_num, previous_info[acc_num]["version"], previous_info[acc_num]["ENA"][1],
                                previous_info[acc_num]["ENA"][2]]
                    if comparing_new_and_old_data(new_data, old_data) == 0:
                        its1_acc = previous_info[acc_num]["ENA"][0]
                        old_acc += 1
                    else:
                        its1_acc = "ITS1DB%s%s" % ("0" * (8 - len(str(starter))), str(starter))
                        new_acc += 1
                else:
                    its1_acc = "ITS1DB%s%s" % ("0" * (8 - len(str(starter))), str(starter))
                lista_dat.append(its1_acc)
                lista_dat.append(s[4])
                lista_dat.append(s[5])
                lista_dat.append(";".join(ena_note))
                starter += 1
                fasta_data[its1_acc] = its1_seq(fasta_diz, acc_num, [s[4], s[5]], ena_note)
                its1_db_data[its1_acc] = [taxid, "ENA"]
            else:
                lista_dat.extend(["Null", "Null", "Null", ";".join(ena_note)])
            if hmm == "True":
                if acc_num in previous_info:
                    new_data = [acc_num, version, s[7], s[8]]
                    old_data = [acc_num, previous_info[acc_num]["version"], previous_info[acc_num]["HMM"][1],
                                previous_info[acc_num]["HMM"][2]]
                    if comparing_new_and_old_data(new_data, old_data) == 0:
                        its1_acc = previous_info[acc_num]["HMM"][0]
                        old_acc += 1
                    else:
                        its1_acc = "ITS1DB%s%s" % ("0" * (8 - len(str(starter))), str(starter))
                        new_acc += 1
                else:
                    its1_acc = "ITS1DB%s%s" % ("0" * (8 - len(str(starter))), str(starter))
                    new_acc += 1
                lista_dat.append(its1_acc)
                lista_dat.append(s[7])
                lista_dat.append(s[8])
                lista_dat.append(";".join(hmm_note))
                starter += 1
                fasta_data[its1_acc] = its1_seq(fasta_diz, acc_num, [s[7], s[8]], hmm_note)
                its1_db_data[its1_acc] = [taxid, "HMM"]
            else:
                lista_dat.extend(["Null", "Null", "Null", "Null"])
            ena2ITSoneDB_accession.append(lista_dat)
    return_time("""--------------------
    Old ITSoneDB accession confirmed %s
    New ITSoneDB accession %s""" % (old_acc, new_acc))
    return ena2ITSoneDB_accession, fasta_data, its1_db_data


def generazione_dizionari_tassonomia(taxnomy_folder):
    # if os.path.exists(os.path.join(taxnomy_folder, "nodes.dmp")):
    #     nodesfile = open(os.path.join(taxnomy_folder, "nodes.dmp"))
    # else:
    #     sys.exit("No NODESFILE")

    if os.path.exists(
            os.path.join(taxnomy_folder, "names.dmp")):
        namesfile = open(
            os.path.join(taxnomy_folder, "names.dmp"))
    else:
        sys.exit("no namesfile")

    # if os.path.exists(os.path.join(taxnomy_folder, "merged.dmp")):
    #     merged = open(os.path.join(taxnomy_folder, "merged.dmp"))
    # else:
    #     sys.exit("no merged")
    #
    # if os.path.exists(os.path.join(taxnomy_folder, "delnodes.dmp")):
    #     deleted = open(os.path.join(taxnomy_folder, "delnodes.dmp"))
    # else:
    #     sys.exit("no deleted")

    # print("Generazione Dizionari")
    node2name = {}
    name2node = {}
    for line in namesfile:
        line = line.strip()
        fields = list(map(str.strip, line.split("|")))
        if fields[3] == "scientific name":
            nodeid, name = fields[0], fields[1]
            node2name[nodeid] = name
            name2node[name] = nodeid
    # node2name["53267"] = "Trebouxia jamesii"

    # node2parent = {}
    # node2order = {}
    # for linea in nodesfile:
    #     linea = linea.strip()
    #     fields = map(strip, linea.split("|"))
    #     nodeid, parentid, ordine = fields[0], fields[1], fields[2]
    #     node2parent[nodeid] = parentid
    #     node2order[nodeid] = ordine
    # node2parent["53267"] = "13786"
    # node2order["53267"] = "species"
    #
    # node2merged = {}
    # for linea in merged:
    #     linea = linea.strip()
    #     fields = map(strip, linea.split("|"))
    #     nodo, new_id = fields[0], fields[1]
    #     node2merged[nodo] = new_id
    #
    # delnodes = []
    # for linea in deleted:
    #     linea = linea.strip()
    #     fields = map(strip, linea.split("|"))
    #     delnodes.append(fields[0])
    return node2name


if __name__ == "__main__":
    param = input_options()
    its1, ena_version, old_map, out, fasta, taxonomy = param.its_file, param.ENA_release_version, param.previous_mapping_file, param.output_folder, param.fasta_index, param.ncbi_taxonomy
    data_dict, fasta_info, its12taxid = create_accession(its1, old_map, acc2fasta_file(fasta))
    name_dict = generazione_dizionari_tassonomia(taxonomy)
    tmp = open(os.path.join(out, "mapping_data_ENA_release_%s.tsv" % ena_version), "w")
    tmp.write(
        "ENA_ACC\tENA_Version\tENA_Annotation_ITS1DB\tENA_start\tENA_end\tENA_Note\tHMM_Annotation_ITS1DB\tHMM_start\tHMM_end\tHMM_Note\n")
    for lista in data_dict:
        tmp.write("%s\n" % "\t".join(lista))
    tmp.close()
    fasta_dir = os.path.join(out, "FASTA_STATICI")
    if os.path.exists(fasta_dir) is False:
        os.mkdir(fasta_dir)
    with gzopen(os.path.join(fasta_dir, "ITSoneDB_total_fasta_rel%s.fa.gz" % ena_version), "wt") as tmp, \
    gzopen(os.path.join(fasta_dir, "ITSoneDB_total_fasta_rel%s_with_flanking.fa.gz" % ena_version), "wt") as flank:
        for acc, seq in fasta_info.items():
            tmp.write(">%s|%s|%s|ITS1 located by %s annotation, %ibp \n%s\n" % (
                acc, its12taxid[acc][0], name_dict[its12taxid[acc][0]], its12taxid[acc][1], len(seq[0]), seq[0]))
            flank.write(">%s|%s|%s|ITS1 located by %s annotation, %ibp \n%s\n" % (
                acc, its12taxid[acc][0], name_dict[its12taxid[acc][0]], its12taxid[acc][1], len(seq[1]), seq[1]))
