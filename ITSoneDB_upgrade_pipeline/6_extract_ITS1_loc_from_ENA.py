#!/usr/bin/env python3
from librs.return_time import return_time
import os
import gzip
import argparse
import argcomplete


def split_options():
    parser = argparse.ArgumentParser(
        description="""The TSV files are parsed out to extract the 
    annotation relative to ITS1 boundaries by means of a 
    commonly used ITS1 synonyms dictionary.
    Usage: python extract_ITS1_loc.py -r nr_ENA_rel -i tsv_folder""",
        prefix_chars="-")
    parser.add_argument("-r", "--release", type=str,
                        help="the number of ena release",
                        action="store", required=True)
    parser.add_argument("-i", "--input", type=str,
                        help="folder of TSV file",
                        action="store", required=True)
    argcomplete.autocomplete(parser)
    return parser.parse_args()


# def index()
#     acc2version = {}
#     with open("index_master", "r") as ind:
#         line = ind.readline()
#         while line:
#             s = line.split("\t")
#             acc2version[s[0]] = s[1]
#     return acc2version


def inner_splitter(open_file):
    with gzip.open(open_file, 'rt') as l:
        acc2data = {}
        data = l.readlines()
    i = 0
    while i < len(data):
        if data[i].startswith("NEW ACC"):
            acc = ""
            acc = data[i+1].strip()
            #print acc
            acc2data.setdefault(acc, [])
            while data[i].strip() != "--------------------------------------------------------":
                acc2data[acc].append(data[i])
                i += 1
            else:
                i += 1
    return acc2data


def aux(i, note):
    if note == " ":
        note = "%i" % i
    else:
        note = "%i::%s" % (i, note)
    return note


def data_parser(diz, ITS1_dictionary):
    #acc2data = index()
    mito_entries = []
    cplast_entries = []
    its1_counter = 0
    for acc, data in diz.items():
        # print('\n'.join(data))
        # break
        organism = ""
        version = data[2].strip()
        taxon = data[-1].strip()
        tax_path = list(map(str.strip, data[5].split(";")))
        #print acc,tax_path
        if "Eukaryota" in tax_path:
            its1_data = []
            i = 0
            note = " "
            #per ogni rigo associato ad un accession number
            for feature_info in data[7:]:
                # print(acc)
                # se non comincia con source vediamo che cominci con uno
                # degli elementi della lista
                if feature_info.split(":location:")[0] in ["misc_feature", "misc_RNA", "precursor_RNA", "rRNA"]:
                    g = list(map(str.strip, feature_info.split("\t")))
                    ps_splt = g[0].split(":location:")
                    p = g[1].split(" ### ")
                    for qual_value in p:
                        #print qual_value
                        #print acc, qual_value.split("=")[-1][1:-1]
                        #se tra gli attributti della riga abbiamo dei
                        #termini del dizionario controllato indicanti
                        #ITS1 allora estraiamo le coordinate
                        q_val = qual_value.split("=")[-1][1:-1]
                        if q_val in ITS1_dictionary:
                            #print acc, qual_value.split("=")[-1][1:-1]
                            if ps_splt[1].startswith("complement"):
                                psl = ps_splt[1].lstrip("complement(")[0:-1].split("..")
                                note = "COMPLEMENT::%s" % q_val
                                #print acc, psl
                            elif ps_splt[1].startswith("join"):
                                psl_t = ps_splt[1].lstrip("join(")[0:-1].split("..")
                                psl = [psl_t[0], psl_t[-1]]
                                note = "JOIN::%s" % q_val
                                #print acc, psl
                            else:
                                psl = ps_splt[1].split("..")
                                note = q_val
                            #print(psl)
                            its1_data.extend(psl)
                            # print(its1_data)
                            # per un rigo potremmo avere piu' attributi
                            # presenti nella lista ma riferiti alla stessa
                            # regione di sequenza, quindi dopo il primo fermo
                            # il ciclo for degli elementi splittati del rigo
                            break
            if len(its1_data) == 2:
                # print(its1_data)
                for feature_info in data[7:]:
                    # se il rigo comincia con source avro' delle info
                    # generali relative la sequenza
                    if feature_info.startswith("source"):
                        g = list(map(str.strip, feature_info.split("\t")))
                        p = g[1].strip().split(" ### ")
                        for line in p:
                            # print line
                            if line.strip("/").startswith("organelle"):
                                c = list(map(str.strip, line.split("=")))
                                # print "+++", c[1]
                                if c[1].strip('"') in mitoconlist:
                                    excl_file.write("%s\n" % " ".join([file_name, "organelle",
                                                    c[1].strip('"'), acc, taxon, *its1_data]))
                                    mito_entries.append(acc)
                                    # print(acc)
                                elif c[1].strip('"') in chloroplist:
                                    excl_file.write("%s\n" % " ".join([file_name,"organelle",
                                                    c[1].strip('"'), acc, taxon, *its1_data]))
                                    cplast_entries.append(acc)
                                # elif line.startswith("db_xref"):
                                #     t = map(strip, line.split("="))
                                #     taxon = t[1].strip('""')
                                #     if taxon.startswith("taxon:"):
                                #         taxon = taxon.strip("taxon:").split('"')[0]
                                break

                # print('yes')
                note = aux(i, note)
                # print(mito_entries)
                if acc not in mito_entries and acc not in cplast_entries:
                    #print acc, taxon, its1_data[0], its1_data[1], note
                    data_feature.write("%s\t%s\t%s\t%s\t%s\t%s\n" %
                                       (acc, version, taxon, its1_data[0],
                                        its1_data[1], note))
                    # print("%s\t%s\t%s\t%s\t%s\t%s\n" %
                    #                    (acc, version, taxon, its1_data[0],
                    #                     its1_data[1], note))
                    its1_counter += 1
                # elif :
                #     #print acc, taxon, its1_data[0], its1_data[1], note
                #     data_feature.write("%s\t%s\t%s\t%s\t%s\t%s\n" %
                #                        (acc, version, taxon, its1_data[0],
                #                         its1_data[1], note))
                # else:
                #     print("organelle", acc, taxon, its1_data)
                i += 1
                #else:
                #    print "lenght of its1_data is ", len(its1_data), "|", its1_data
    return mito_entries, cplast_entries, its1_counter


if __name__ == '__main__':
    return_time("Extraction of ITS1 locations from ENA feature tables START")
    opts = split_options()
    release, tsv = opts.release, opts.input
    if tsv.endswith('/'):
        base = tsv[:-1]
    else:
        base = tsv
    base = "/".join(base.split('/')[:-3])
    # print(base)
    # print(release, tsv)
    # if os.getcwd().split("/")[-1] == release:
    ITS1_dict_path = "/".join([base, "/aux/ITS1_dictionary"])
    # print(ITS1_dict_path)
    # else:
    #     sys.exit("Please launch this program in the %s release directory!" % release)
    mitoconlist = ["mitochondria", "mitochondrion", "mitochondrion:kinetoplast"]
    chloroplist = ["chloroplast", "chloroplasts", "plastid:chloroplast"]
    ITS1_dictionary = []
    with open(ITS1_dict_path) as l:
        for line in l:
            line = line.strip()
            line = line.replace("\'", "")
            line = line.replace("[", "")
            line = line.replace("]", "")
            ITS1_dictionary.append(line)
    output = os.path.join(base, 'releases', release, "ITS1_loc_ENA.csv")
    excluded = os.path.join(base, 'releases', release, "excluded_ITS1_ENA.csv")
    excl_file = open(excluded, "wt")
    with open(output, "wt") as data_feature:
        data_feature.write("ACCESSION\tVERSION\tTAX_ID\tITS1_FROM\tITS1_TO\tNOTE\n")
        for file_name in os.listdir(tsv):
            if file_name.endswith(".tsv.gz"):
                # file_name = 'rel_std_hum_22_r141.tsv.gz'
                # print("start %s" % file_name)
                mito_entries, cplast_entries, its1_counter = \
                    data_parser(inner_splitter(os.path.join(tsv, file_name)),
                                ITS1_dictionary)
                return_time("%s resume:\nITS1 found -> %s\nMitoc. entries -> %s\nCplast. entries -> %s\n" % (
                    file_name,its1_counter, len(mito_entries), len(cplast_entries)))
    excl_file.close()
    return_time("Extraction of ITS1 locations from ENA feature tables DONE")