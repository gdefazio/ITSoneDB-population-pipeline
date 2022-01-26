#!/usr/bin/env python3
import os
import gzip
from Bio.SeqIO.FastaIO import SimpleFastaParser
import argparse
import argcomplete
import sys
from librs.return_time import return_time


def split_options():
    parser = argparse.ArgumentParser(description="split fasta and tsv files", prefix_chars="-")
    parser.add_argument("-t", "--tsv_folder", type=str,
                        help="folder containing the tsv files",
                        action="store", required=True)
    parser.add_argument("-f", "--fasta_folder", type=str,
                        help="folder containing the fasta files",
                        action="store", required=True)
    parser.add_argument("-i", "--its_file", type=str,
                        help="tsv containg ITS1 location per accession",
                        action="store", required=True)
    parser.add_argument("-o", "--output_folder", type=str, help="output folder", action="store", required=False,
                        default="etl")
    argcomplete.autocomplete(parser)
    return parser.parse_args()


def acc_list(input_its1_tsv):
    accession_list = []
    with open(input_its1_tsv, 'rt') as a:
        a.readline()
        for line in a:
            accession_list.append(list(map(str.strip, line.split("\t")))[0])
    return accession_list


def tsv_splitting(input_tsv, out_folder, accession):
    tsv_folder = os.path.join(out_folder, "splitted_tsv")
    relative = "etl%s" % tsv_folder.split("/etl")[-1]
    if os.path.exists(tsv_folder) is False:
        os.mkdir(tsv_folder)
    acc2file = {}
    if os.path.exists(input_tsv) is False:
        sys.exit("The selected TSV folder does not exist")
    for name in os.listdir(input_tsv):
        return_time("%s elaboration..." % name)
        data_diz = dict()
        if name.endswith(".tsv.gz"):
            # import gzip
            with gzip.open(os.path.join(input_tsv, name), 'rt') as a:
                stringa = list()
                for line in a:
                    if line.strip() != "--------------------------------------------------------":
                        stringa.append(line)
                    else:
                        stringa.append(line)
                        acc_num = stringa[1].strip().split(".")[0]
                        data_diz[acc_num] = stringa
                        stringa = list()
        elif name.endswith(".tsv"):
            with open(os.path.join(input_tsv, name), 'rt') as a:
                stringa = list()
                for line in a:
                    if line.strip() != "--------------------------------------------------------":
                        stringa.append(line)
                    else:
                        stringa.append(line)
                        # print count
                        acc_num = stringa[1].strip().split(".")[0]
                        data_diz[acc_num] = stringa
                        stringa = list()
        for acc_num in accession:
            if data_diz.get(acc_num):
                tmp_file = gzip.open(os.path.join(tsv_folder, "%s.gz" % acc_num), "wt")
                tmp_file.write("".join(data_diz[acc_num]))
                tmp_file.close()
                acc2file[acc_num] = os.path.join(relative, "%s.gz" % acc_num)
        data_diz.clear()
    return tsv_folder, acc2file


def split_fasta(input_fasta, out_folder, accession):
    fasta_folder = os.path.join(out_folder, "splitted_fasta")
    relative = "etl%s" % fasta_folder.split("/etl")[-1]
    if os.path.exists(fasta_folder) is False:
        os.mkdir(fasta_folder)
    acc2file = {}
    if os.path.exists(input_fasta) is False:
        sys.exit("The selected FASTA folder does not exist")
    for name in os.listdir(input_fasta):
        data_diz = {}
        return_time("%s elaboration..." % name)
        if name.endswith(".fasta.gz"):
            # import gzip
            with gzip.open(os.path.join(input_fasta, name), 'rt') as a:
                for title, seq in SimpleFastaParser(a):
                    acc_num = title.split("|")[0].split(".")[0]
                    data_diz[acc_num] = [title, seq]
        elif name.endswith(".fasta"):
            with open(os.path.join(input_fasta, name),'rt') as a:
                for title, seq in SimpleFastaParser(a):
                    acc_num = title.split("|")[0].split(".")[0]
                    data_diz[acc_num] = [title, seq]
        for acc_num in accession:
            if data_diz.get(acc_num):
                tmp_file = gzip.open(os.path.join(fasta_folder, "%s.gz" % acc_num), "wt")
                tmp_file.write(">%s\n%s\n" % (data_diz[acc_num][0], data_diz[acc_num][1]))
                tmp_file.close()
                acc2file[acc_num] = os.path.join(relative, "%s.gz" % acc_num)
        data_diz.clear()
    return fasta_folder, acc2file


if __name__ == "__main__":
    param = split_options()
    tsv, fasta, out, its1 = param.tsv_folder, param.fasta_folder, param.output_folder, param.its_file
    lista_accession = acc_list(its1)
    if os.path.exists(out) is False:
        os.mkdir(out)
    cartella_tsv, diz_tsv = tsv_splitting(tsv, out, lista_accession)
    return_time("Splitted tsv stored in %s" % cartella_tsv)
    with gzip.open(os.path.join(out, "acc2tsv_info.csv.gz"), "wt") as tmp:
        tmp.write("ACC\tTSV\n")
        for acc in lista_accession:
            tmp.write("%s\t%s\n" % (acc, diz_tsv[acc]))
    cartella_fasta, diz_fasta = split_fasta(fasta, out, lista_accession)
    return_time("Splitted fasta stored in %s" % cartella_fasta)
    with gzip.open(os.path.join(out, "acc2fasta_info.csv.gz"), "wt") as tmp:
        tmp.write("ACC\tFASTA\n")
        for acc in lista_accession:
            tmp.write("%s\t%s\n" % (acc, diz_fasta[acc]))
