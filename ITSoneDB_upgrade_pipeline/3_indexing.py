#!/usr/bin/env python3
from librs.return_time import return_time
import os
import gzip
import argparse
import argcomplete


#lanciare il programma in release_138

def split_options():
    parser = argparse.ArgumentParser(
        description="The program return an index file of all downloaded entries",
        prefix_chars="-")
    parser.add_argument("-r", "--release", type=str,
                        help="the number of ena release",
                        action="store", required=True)
    parser.add_argument("-i", "--input", type=str,
                        help="folder of FASTA file",
                        action="store", required=True
                        )
    argcomplete.autocomplete(parser)
    return parser.parse_args()


if __name__ == '__main__':
    return_time("Index generation")
    opts = split_options()
    release, input = opts.release, opts.input
    if input.endswith('/'): input = input[:-1]
    base = "/".join(input.split('/')[:-1])
    # if os.getcwd().split("/")[-1] == release:
    try:
        os.mkdir(os.path.join(base, "index_%s" % release))
        return_time("%s created" % os.path.join(base, "index_%s" % release))
    except OSError:
        return_time("%s already exists" % os.path.join(base, "index_%s" % release))
    # else:
    #     sys.exit("Please start this program in the release %s directory!" % release)
    return_time("START fasta indexing elaboration")
    check = 0
    with gzip.open(os.path.join(base, "index_%s/index.csv.gz" % release), "wt") as index:
        index.write("ACCESSION\tVERSION\tTAX_ID\tLENGTH\tADDRESS\n")
        for fasta in os.listdir(input):
            if fasta.startswith("rel_") and fasta.endswith("_r%s.fasta.gz" % release):
                # return_time("START %s " % fasta)
                #path = os.path.join(input,fasta)
                # print(input, fasta)
                with gzip.open(os.path.join(input,fasta), 'rt') as fastoso:
                    line = fastoso.readline()
                    while line is not "":
                        if line.startswith(">"):
                            key = line[1:].strip().split(".")
                            acc = key[0]
                            version = key[1].split("|")[0]
                            tax_id = key[1].split("|")[1]
                            line = fastoso.readline()
                            val = len(line.strip())
                            index.write("%s\t%s\t%s\t%s\t%s\n" % (
                                acc, version, tax_id, val, fasta.split(".")[0]))
                            line = fastoso.readline()
                check += 1
    if len(os.listdir(input)) == check:
        return_time("""----------------------------------
        All of %s have just been elaborated""" % check)
    else:
        return_time("""
        Error: Only %s of %s have just been elaborated""" % (check, len(os.listdir(input))))

