import os
import gzip
import sys
from string import strip

mito_entries = []
mitoconlist = ["mitochondria", "mitochondrion","mitochondrion:kinetoplast"]
data_feature = open("ITS1_loc.txt", "w")
ITS1_dictionary =[]
with open("ITS1_dictionary") as l:
    for line in l:
        line = line.strip()
        line = line.replace("\'","")
        line = line.replace("[","")
        line = line.replace( "]", "" )
        ITS1_dictionary.append(line)


product_note = set()

folder = sys.argv[1]

def inner_splitter(open_file):
    l = gzip.open(open_file)
    acc2data = {}
    data = l.readlines()
    l.close()
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

def data_parser(diz):
    for acc, data in diz.items():
        organism = ""
        taxon = ""
        p = ""
        tax_path = map(strip,data[4].split(";"))
        print acc,tax_path
        its1_data = []
        if "Eukaryota" in tax_path:
            for feature_info in data[6:]:
                p = map(strip, feature_info.split("\t"))
                if p[0].startswith("source"):
                    s = map(strip, p[1].split("/"))
                    organism = (s[1].strip("###"))
                    for line in s:
                        if line.startswith("organelle"):
                            c = map(strip, line.split("="))
                            if c[1].strip("") in mitoconlist:
                                mito_entries.append( acc )
                        elif line.startswith("db_xref"):
                            t = map(strip, line.split("="))
                            taxon = t[1].strip('""')
                            if taxon.startswith("taxon"):
                                print taxon
                if p[0].split(":location:")[0] in ["misc_feature", "misc_RNA", "precursor_RNA", "rRNA"]:
                    c = map(strip, p[1].split("###"))
                    #print c
                    for qual_value in c:
                        if qual_value.split("=")[-1][1:-1] in ITS1_dictionary:
                            its1_data.append("%s:::%s" % (p[0],qual_value.split("=")[-1][1:-1]))
            if acc not in mito_entries and len(its1_data) > 0:
                data_feature.write("%s\t%s\t%s\t%s\n" % (acc, taxon, organism, "\t".join(its1_data)))
            else:
                print "%s\t%s\t%s\t%s\n" % (acc, taxon, organism, "\t".join(its1_data))



for file_name in os.listdir(folder):
    if file_name.endswith("gz"):
        data_parser(inner_splitter(os.path.join(folder,file_name)))

data_feature.close()


