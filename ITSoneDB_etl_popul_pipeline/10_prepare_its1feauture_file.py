#!/usr/bin/env python3

from librs.return_time import return_time
import os
import argparse
import argcomplete


def input_options():
    parser = argparse.ArgumentParser(description="split fasta and tsv files", prefix_chars="-")
    parser.add_argument("-i", "--mapping_file", type=str,
                        help="mapping file created by the 2_create_mapping_ENA_ITSoneDB_acc.py",
                        action="store", required=True)
    parser.add_argument("-l", "--its1_loc_file", type=str,
                        help="tsv containg ITS1 location and annotation obtained from ENA",
                        action="store", required=True)
    # parser.add_argument("-r", "--representative", type=str,
    #                     help="list of representative sequences",
    #                     action="store", required=True)
    parser.add_argument("-o", "--output_folder", type=str, help="output folder", action="store", required=False,
                        default="etl")
    argcomplete.autocomplete(parser)
    return parser.parse_args()


# noinspection PyPep8Naming
def dataITS1(its1_tsv):
    acc2data = {}
    with open(its1_tsv, 'rt') as a:
        a.readline()
        for line in a:
            s = list(map(str.strip, line.split("\t")))
            acc2data[s[0]] = s
    return acc2data


# def species_rep(rep_list):
#     rep = []
#     with open(rep_list) as a:
#         for line in a:
#             s = map(strip, line.split("|"))
#             rep.append(s[1].split("|")[0])
#     return rep


def acc2annotazione(loc_file):
    acc2ann = {}
    with open(loc_file, 'rt') as a:
        for line in a:
            s = list(map(str.strip, line.split("\t")))
            acc_num, annot = s[0], list(map(str.strip, s[-1].split(":")))
            acc2ann[acc_num] = [annot[0], annot[-1]]
    return acc2ann


if __name__ == "__main__":
    return_time("Start 10_prepare_its1feature_file")
    param = input_options()
    its1, tsv_data, out = param.mapping_file, param.its1_loc_file, param.output_folder
    acc2annotation = acc2annotazione(tsv_data)
    # representative_list = species_rep(representative)
    tmp = open(os.path.join(out, "its1feature.tsv"), "wt")
    tmp.write(
        "GBentry_Accession\tAccessionVersion\thasGBannotation\tGBITS1Accession\tGBkey\tGBQualifier\tGBValue\tGBstart\tGBend\tGBleftComplete\tGBrightComplete\thasHMM\tHMMITS1Accession\tHMMstart\tHMMend\tHMMleftComplete\tHMMrightComplete\trepresentativeForSpecie\n")
    ena2itsonedb_data = dataITS1(its1)
    for acc, data_list in ena2itsonedb_data.items():
        stringa = [acc, data_list[1]]
        if data_list[2] != "Null":
            stringa.append(1)
            stringa.append(data_list[2])
            stringa.append(acc2annotation[acc][0])
            stringa.append("product")
            stringa.append(acc2annotation[acc][1])
            stringa.append(data_list[3].replace("<","").replace(">",""))
            stringa.append(data_list[4].replace("<","").replace(">",""))
            if data_list[3].find("<") != -1 or data_list[3].find(">"):
                stringa.append(1)
            else:
                stringa.append(0)
            if data_list[4].find("<") != -1 or data_list[4].find(">"):
                stringa.append(1)
            else:
                stringa.append(0)
        else:
            stringa.append(0)
            stringa.append(0)
            stringa.append(0)
            stringa.append(0)
            stringa.append(0)
            stringa.append(0)
            stringa.append(0)
            stringa.append(0)
            stringa.append(0)
        if data_list[6] != "Null":
            stringa.append(1)
            stringa.append(data_list[6])
            stringa.append(data_list[7])
            stringa.append(data_list[8])
            stringa.append(1)
            stringa.append(1)
        else:
            stringa.append(0)
            stringa.append(0)
            stringa.append(0)
            stringa.append(0)
            stringa.append(0)
            stringa.append(0)
        stringa.append(0)
        stringa = map(str, stringa)
        tmp.write("\t".join(stringa) + "\n")
    tmp.close()
    return_time('DONE')