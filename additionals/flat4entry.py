import os
import gzip
import datetime


def input_options():
    parser = argparse.ArgumentParser(
        description="This program splits the DAT files only for the accession "
                    "contained in ITSoneDB_entries_generator_definitive.csv only.",
        prefix_chars="-")
    parser.add_argument("-i", "--its_file", type=str,
                        help="tsv containing ITS1 location per accession",
                        action="store", required=True)
    parser.add_argument("-u", "--unsplit", type=str,
                        help="The unsplitted DAT directory",
                        action="store", required=True)
    parser.add_argument("-o", "--output_dir", type=str,
                        help="output folder in which store flat directory",
                        action="store", required=False,
                        default="etl")
    argcomplete.autocomplete(parser)
    return parser.parse_args()


def import_index(path_index):
    index = {}
    #path = "results_138/ITSoneDB_entries_generator_definitive_datname.csv"
    with open(path_index, "r") as findex:
        findex.readline()
        line = findex.readline()
        while line:
            s = line.strip().split("\t")
            index.setdefault(s[-1], [])
            index[s[-1]].append(s[0])
            line = findex.readline()
    return index


if __name__ == '__main__':
    opts = input_options()
    its1, dat_dir, outdir = opts.its1_file, opts.unsplit, opts.output_dir
    index = import_index(its1)
    for dat in index.keys():
        path = os.path.join(dat_dir, "%s.dat.gz" % dat)
        dat_index = {}
        print "start ", dat, " indexing", datetime.datetime.now().time()
        with gzip.open(path, "r") as source:
            c = 0
            line = source.readline()
            while line:
                if line.startswith("ID"):
                    acc = line.split(";")[0].split("   ")[1]
                    start = c
                    line = source.readline()
                    c += 1
                elif line.startswith("//"):
                    stop = c
                    dat_index[acc] = (start, stop)
                    line = source.readline()
                    c += 1
                else:
                    line = source.readline()
                    c += 1
        print "stop ", dat, " indexing", datetime.datetime.now().time()
        accession = index[dat]
        print "start ", dat, " flat generation", datetime.datetime.now().time()
        with gzip.open(path, "r") as source:
            lines = source.readlines()
        for acc in accession:
            with open(os.path.join(outdir, "flat/%s" % acc), "w") as flat:
                start, stop = dat_index[acc]
                for line in lines[start:stop + 1]:
                    flat.write(line)
        print "stop ", dat, " flat generation", datetime.datetime.now().time()
