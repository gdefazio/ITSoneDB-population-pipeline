import os
import gzip
from string import strip
import sys
from datetime import datetime
import argparse
import argcomplete

#lanciare il release_138

def input_options():
    parser = argparse.ArgumentParser(
        description="This program splits the alignment files only for the accession "
                    "contained in ITSoneDB_entries_generator_definitive.csv only if the "
                    "accession has an hmm-based inference",
        prefix_chars="-")
    parser.add_argument("-i", "--its_file", type=str,
                        help="tsv containing ITS1 location per accession",
                        action="store", required=True)
    parser.add_argument("-u", "--unsplit", type=str,
                        help="The unsplitted alignment directory",
                        action="store", required=True)
    parser.add_argument("-o", "--output_dir", type=str,
                        help="output folder in which store hmm_match_data directory",
                        action="store", required=False,
                        default="etl")
    argcomplete.autocomplete(parser)
    return parser.parse_args()


def verify_outpath_exist(outdir):
    try:
        os.mkdir(outdir)
    except OSError:
        print "+++WARNING+++\nThe output directory already exists:\nPlease " \
              "rename or remove ./hmm_match_data directory!"
        sys.exit(1)
    try:
        os.mkdir("%s/5_8S_alignment_file" % outdir)
    except OSError:
        pass
    try:
        os.mkdir("%s/18S_alignment_file" % outdir)
    except OSError:
        pass


def import_DB_acc(input):
    accession_list = []
    ITSoneDB_entries_file = input
    with open(ITSoneDB_entries_file, "r") as a:
        a.readline()
        for line in a:
            s = map(strip, line.split("\t"))
            if s[6] == "True":
                accession_list.append(s[0])
    return accession_list


def extract_align_5_8(accession_list, folder, outdir):
    p = os.path.join(os.path.join(unsplitted_dir, folder), "align_5_8S_data.txt.gz")
    print p
    print "5.8S ALIGNEMENT EXTRACTION", datetime.now().time()
    counter = []
    with gzip.open(p, "r") as f:
        line = f.readline()
        while line.startswith(">>") is False:
            if line.startswith("Internal pipeline statistics summary:") is False:
                line = f.readline()
            else:
                count = len(set(counter))
                print "5.8S accession found: %i" % count
                return
        else:
            while line.startswith("Internal pipeline statistics summary:") is False:
                if line.startswith(">>"):
                    s = map(strip, line.split(" "))
                    acc = s[1]
                    #print acc
                    if acc in accession_list:
                        counter.append(acc)
                        #print ">>>>", acc
                        with open("%s/5_8S_alignment_file/%s.align" % (outdir, acc), "a") as tmp:
                            tmp.write(line)
                            #print line
                            line = f.readline()
                            while line.startswith(">>") == 0: #or (line.startswith("Internal pipeline statistics summary:") == 0):
                                if line.startswith("Internal pipeline statistics summary:") is False:
                                    tmp.write(line)
                                    #print line
                                    line = f.readline()
                                    #print "last line", line
                                else:
                                    count = len(set(counter))
                                    print "5.8S accession found: %i" % count
                                    return
                    else:
                        line = f.readline()
                else:
                    line = f.readline()
    count = len(set(counter))
    print "5.8S accession found: %i" % count


def extract_align_18(accession_list, folder, outdir):
    print "18S ALIGNEMENT EXTRACTION", datetime.now().time()
    p = os.path.join(os.path.join(unsplitted_dir, folder), "align_18S_data.txt.gz")
    #print p
    counter = []
    with gzip.open(p, "r") as f:
        line = f.readline()
        while line.startswith(">>") is False:
            if line.startswith("Internal pipeline statistics summary:") is False:
                line = f.readline()
            else:
                count = len(set(counter))
                print "18S accession found: %i" % count
                return 
        else:
            while line.startswith("Internal pipeline statistics summary:") is False:
                if line.startswith(">>"):
                    s = map(strip, line.split(" "))
                    acc = s[1]
                    #print acc
                    if acc in accession_list:
                        counter.append(acc)
                        #print ">>>>", acc
                        with open("%s/18S_alignment_file/%s.align" % (outdir, acc), "a") as tmp:
                            tmp.write(line)
                            line = f.readline()
                            while line.startswith(">>") is False:
                                if line.startswith("Internal pipeline statistics summary:") is False:
                                    tmp.write(line)
                                    line = f.readline()
                                else:
                                    count = len(set(counter))
                                    print "18S accession found: %i" % count
                                    return
                    else:
                        line = f.readline()
                else:
                    line = f.readline()
    count = len(set(counter))
    print "18S accession found: %i" % count


if __name__ == '__main__':
    opts = input_options()
    final_table, unsplitted_dir, output_dir = opts.its_file, opts.unsplit, opts.output_dir
    outdir = os.path.join(output_dir, "hmm_match_data")
    accession_list = import_DB_acc(final_table)
    verify_outpath_exist(outdir)
    for folder in os.listdir(unsplitted_dir):
    #folder = "rel_gss_env_01_r138"
        extract_align_5_8(accession_list, folder, outdir)
        extract_align_18(accession_list, folder, outdir)
    print "DONE"
