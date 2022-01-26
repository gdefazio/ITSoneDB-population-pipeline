#!/usr/bin/env python3
# import os
import gzip
# import sys
# from string import strip
from librs.return_time import return_time
import argparse
import argcomplete


#lanciare in release_138

def split_options():
    parser = argparse.ArgumentParser(
        description="Merge the ITS1 table by HMM procedure with the ITS1 table by ENA procedure.",
        prefix_chars="-")
    parser.add_argument("-r", "--release", type=str,
                        help="the number of ena release",
                        action="store", required=True)
    parser.add_argument("-m", "--hmm", type=str,
                        help="The output file of HMM procedure",
                        action="store", required=True
                        )
    parser.add_argument("-e", "--ena", type=str,
                        help="The output file od ENA procedure",
                        action="store", required=True)
    parser.add_argument("-i", "--index", type=str, help="The path of index file",
                        action="store", required=True)
    parser.add_argument("-o", "--out_file", type=str, help="output file path name",
                         action="store", required=True)
    argcomplete.autocomplete(parser)
    return parser.parse_args()


# def inner_splitter(open_file):
#     "Inner data_parser function"
#     with gzip.open(open_file) as l:
#         acc2data = {}
#         data = l.readlines()
#     i = 0
#     while i < len(data):
#         if data[i].startswith("NEW ACC"):
#             acc = ""
#             acc = data[i+1].strip()
#             #print acc
#             acc2data.setdefault(acc, [])
#             while data[i].strip() != "--------------------------------------------------------":
#                 acc2data[acc].append(data[i])
#                 i += 1
#             else:
#                 i += 1
#     return acc2data


def index_reader(_index):
    """
    The Extraction of the version sequence, Taxon ID and .dat file name
    of ENA sequences by index file.
    """
    acc2taxon = {}
    with gzip.open(_index, "rt") as index:
        line = index.readline()
        while line:
            s = line.strip().split("\t")
            acc2taxon[s[0]] = (s[1], s[2], s[3], s[4])
            line = index.readline()
    return acc2taxon


def hmm_its_parser(_hmm):
    """
    Extraction of information about the start and the end of the ITS1 sequences
    inferred by the hmm alignments.
    """
    with open(_hmm, "r") as handle:
        readlines = handle.readlines()
    acc2data = {} #accession, version, tax_id, ITS1_from, ITS1_to, note
    for line in readlines[1:]:
        s = line.strip().split("\t")

        acc = s[0].split('.')[0]
        vers = s[1]
        taxid = s[2]
        # print(acc)
        # print(acc[1].split('|'))
        # vers, taxid = s[0].split('.')[1].split('|')
        acc2data.setdefault(acc, {})
        if "+" in s[-2]:
            its1_from = int(s[8]) + 1
            its1_to = int(s[12]) - 1
            strand = "FWD"
        else:
            its1_from = int(s[8]) - 1
            its1_to = int(s[12]) + 1
            strand = "RVS"
        acc2data[acc] = [vers, taxid, its1_from, its1_to, strand]
    return acc2data


# def hmm_acc_data_dct(_hmm, _index):
#     "To associate the Taxon ID to the ITS1 sequences inferred by the HMM alignments."
#     hmm_acc2data = {}
#     with open(_hmm, "r") as handle:
#         lines = handle.readlines()
#     div2acc2data = hmm_its_parser(lines)
#     acc2taxon = index_reader(_index)
#     #folder = "tsv_138"
#     for divk in div2acc2data.keys(): #os.listdir(folder):
#         #file_name = "rel_%s_r138_output.tsv.gz" % divk
#         #acc2taxon = data_parser(inner_splitter(os.path.join(folder, file_name)))
#         for acc in div2acc2data[divk]:
#                         #accession, version, tax_id, ITS1_from, ITS1_to, note
#             hmm_acc2data[acc] = [acc, acc2taxon[acc][0], acc2taxon[acc][1], div2acc2data[divk][acc][0],
#                                  div2acc2data[divk][acc][1], div2acc2data[divk][acc][2]]
#             # if len(hmm_acc2data[acc]) != 5:
#             #     print acc, len(hmm_acc2data[acc]), hmm_acc2data[acc]
#     return hmm_acc2data, acc2taxon


def ena_acc_data_dct(_ena):
    """Extraction of the information about the ITS1 sequences inferred by the ENA
    annotations"""
    ena_acc2data = {}
    with open(_ena, "r") as handle:
        lines = handle.readlines()
    for line in lines[1:]:
        s = line.strip().split("\t")
        ena_acc2data.setdefault(s[0], [])
        ena_acc2data[s[0]].append(s)
        # if len(ena_acc2data[s[0]]) != 4:
        #     print(ena_acc2data[s[0]])
    return ena_acc2data


if __name__ == '__main__':
    opts = split_options()
    release, index, hmm_table, \
    ena_table, _output = opts.release, opts.index, opts.hmm, opts.ena, opts.out_file
    ena_acc2data = ena_acc_data_dct(ena_table)
    #print "start"
    hmm_acc2data = hmm_its_parser(hmm_table)
    acc2taxon = index_reader(index)
    #hmm_acc2data, acc2taxon = hmm_acc_data_dct(hmm_table, index)
    #print "stop"
    all_acc = []
    all_acc.extend(ena_acc2data.keys())
    all_acc.extend(hmm_acc2data.keys())
    all_acc_set = set(all_acc)
    with open(_output, "w") as output:
        output.write("ACCESSION\tVERSION\tTAX_ID\tHAS_ENA\tITS1_F_ENA\tITS1_T_ENA\tHAS_HMM\tITS1_F_HMM\tITS1_T_HMM\tNOTE\tDAT_FILENAME\n")
        tt = 0
        tf = 0
        ft = 0
        for acc in all_acc_set:
            #IN ENTRAMBE
            ena = "False"
            try:
                version_e, tax_id_e, f_e, t_e, n_e = ena_acc2data[acc][0][1:]
                ena = "True"
                version_h, tax_id_h, f_h, t_h, n_h = hmm_acc2data[acc]
                note = ""
                if len(ena_acc2data[acc]) == 1:
                    if n_e == "0":
                        note = n_h
                    elif "0::" in n_e:
                        n_e = n_e.split(":")
                        note = "%s|%s" % (n_e[-1], n_h)
                else:
                    if n_e == "0":
                        note = "%s|MULTI" % n_h
                    elif "0::" in n_e:
                        n_e = n_e.split(":")
                        note = "%s|%s|MULTI" % (n_e[-1], n_h)
                output.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %
                             (acc, version_e, tax_id_h, "True", f_e, t_e,
                              "True", f_h, t_h, note, acc2taxon[acc][-1]))
                tt += 1
            except KeyError:
                note = ""
                if ena == "True":
                    version_e, tax_id_e, f_e, t_e, n_e = ena_acc2data[acc][0][1:]
                    if len(ena_acc2data[acc]) == 1:
                        if n_e == "0":
                            note = " "
                        elif "0::" in n_e:
                            n_e = n_e.split(":")
                            note = n_e[-1]
                    else:
                        if n_e == "0":
                            note = "MULTI"
                        elif "0::" in n_e:
                            n_e = n_e.split(":")
                            note = "%s|MULTI" % n_e[-1]
                    output.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %
                                 (acc, version_e, tax_id_e, "True", f_e, t_e,
                                  "False", " ", " ", note, acc2taxon[acc][-1]))
                    tf += 1
                else:
                    version_h, tax_id_h, f_h, t_h, n_h = hmm_acc2data[acc]
                    output.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %
                                 (acc, version_h, tax_id_h, "False", " ", " ",
                                  "True", f_h, t_h, n_h, acc2taxon[acc][-1]))
                    ft += 1

        return_time("""Resume of HMM and ENA procedure
        ACCESSION ENA-TRUE HMM-TRUE %i
        ACCESSION ENA-TRUE HMM-FALSE %i
        ACCESSION ENA-FALSE HMM-TRUE %i""" % (tt, tf, ft))

