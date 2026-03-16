import os
import argparse
import argcomplete
import datetime

def split_options():
    parser = argparse.ArgumentParser(description="Return an index table", prefix_chars="-")
    # parser.add_argument("-r", "--release", type=str,
    #                     help="the number of ena release",
    #                     action="store", required=True)
    parser.add_argument("-a", "--acc_db", type=str,
                        help="path of table with ITSoneDB entries",
                        action="store", required=True,
                        default="results_138/ITSoneDB_entries_generator_definitive.csv")
    parser.add_argument("-i", "--index", type=str,
                        help="path of index table with ENA entrie",
                        action="store", required=False, default="index/index.csv")
    parser.add_argument("-o", "--out_file", type=str, help="output file path name",
                        action="store", required=False,
                        default="index/final_index.csv")
    argcomplete.autocomplete(parser)
    return parser.parse_args()


def import_acc_DB(ITSoneDB_entries_file):
    acc_DB_dict = {}
    with open(ITSoneDB_entries_file, "r") as acc_DB:
        acc_DB.readline()
        line = acc_DB.readline()
        while line:
            s = line.strip().split("\t")
            acc_DB_dict.setdefault(s[-1], [])
            acc_DB_dict[s[-1]].append(s[0])
            line = acc_DB.readline()
    return acc_DB_dict


def final_index(ITSoneDB_entries_file, ena_index, final_index):
    accession_dict = import_acc_DB(ITSoneDB_entries_file)
    with open(ena_index, "r") as index:
        with open(final_index, "w") as final_index:
            final_index.write("ACCESSION\tVERSION\tITSoneDB\tTAX_ID\tSEQ_LEN\tADDRESS\n")
            index.readline()
            line = index.readline()
            while line:
                address, acc, version, taxid, seqlen = line.strip().split("\t")
                try:
                    accession_list = accession_dict[address]
                    if acc in accession_list:
                        final_index.write("%s\t%s\t1\t%s\t%s\t%s\n" % (
                            acc, version, taxid, seqlen, address))
                    else:
                        final_index.write("%s\t%s\t0\t%s\t%s\t%s\n" % (
                            acc, version, taxid, seqlen, address))
                    line = index.readline()
                except KeyError:
                    final_index.write("%s\t%s\t0\t%s\t%s\t%s\n" % (
                        acc, version, taxid, seqlen, address))
                    line = index.readline()



if __name__ == '__main__':
    print "Start final index generation..."
    print datetime.datetime.now().time()
    param = split_options()
    acc_db_table, index_ena, \
    outfolder = param.acc_db, param.index, param.out_file
    final_index(acc_db_table, index_ena, outfolder)
    print "DONE"
    print datetime.datetime.now().time()



