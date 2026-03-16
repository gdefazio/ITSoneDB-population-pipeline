#!/usr/bin/env python3
from librs.return_time import return_time
import os
import argparse
import argcomplete
from Bio import SeqIO
from gzip import open as gzopen


def input_options():
    parser = argparse.ArgumentParser(description="split fasta and tsv files", prefix_chars="-")
    parser.add_argument("-r", "--release", type=str,
                        help="ENA release",
                        action="store", required=True)
    parser.add_argument("-c", "--representative_file", type=str,
                        help="TSV file contatining the representative ITS1 sequence per species taxid",
                        action="store", required=True)
    parser.add_argument("-o", "--output_folder", type=str, help="output folder",
                        action="store", required=False,
                        default=os.getcwd())
    parser.add_argument("-i", "--its1_fasta", type=str,
                        help="FASTA_file containing ITS1 sequences",
                        action="store", required=True)
    parser.add_argument("-f", "--its1_flanking", type=str,
                        help="FASTA_file containing ITS1 sequences with flanking regions",
                        action="store", required=True)
    argcomplete.autocomplete(parser)
    return parser.parse_args()


def rep_list(rep_tsv):
    lista = set()
    with open(rep_tsv, 'rt') as a:
        for line in a:
            s = list(map(str.strip, line.split("\t")))
            lista.add(s[1].split('|')[0])
    return lista


if __name__ == "__main__":
    return_time("Start 7_rep_fasta_generation")
    param = input_options()
    release, representative_file, its1_fasta, flanking_fasta, out = param.release, param.representative_file, param.its1_fasta, param.its1_flanking, param.output_folder
    representative_set = rep_list(representative_file)
    with gzopen(its1_fasta, 'rt') as fasta, \
            gzopen(os.path.join(out, "ITSoneDB_rep_seq_r%s.fasta.gz" % release), "wt") as tmp:
        for record in SeqIO.parse(fasta, "fasta"):
            if record.name.split('|')[0] in representative_set:
                tmp.write(record.format("fasta"))
    with gzopen(flanking_fasta, 'rt') as flanking_fasta, \
            gzopen(os.path.join(out, "ITSoneDB_rep_seq_and_flanking_r%s.fasta.gz" % release), "wt") as tmp:
        for record in SeqIO.parse(flanking_fasta, "fasta"):
            if record.name.split('|')[0] in representative_set:
                tmp.write(record.format("fasta"))
    return_time("Stop 7_rep_fasta_generation")