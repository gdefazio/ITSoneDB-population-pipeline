import os
from string import strip
import sys
import gzip

tsv = sys.argv[1]
check_folder = sys.argv[2]


def acc2align_folder(tsv, release_n, align_f, lista_accession):
    """

    :param align_f: str
    :param tsv: str
    :param release_n: int
    :return :
    """
    folder2acc_list = {}
    counter = 0
    with open(tsv) as a:
        a.readline()
        for line in a:
            s = map(strip, line.split("\t"))
            if s[0].split("|")[0] in lista_accession:
                counter += 1
                folder2acc_list.setdefault("{}/rel_{}_r{}".format(align_f, s[-1], release_n), [])
                folder2acc_list["{}/rel_{}_r{}".format(align_f, s[-1], release_n)].append(s[0])
    print counter
    return folder2acc_list


def parse_file(file_obj, acclist):
    acc2_alignment_data = {}
    line = file_obj.readline()
    while line:
        if len(acc2_alignment_data) == len(acclist):
            print len(acc2_alignment_data), len(acclist)
            break
        elif line.startswith(">>"):
            s = map(strip, line.split(" "))
            acc_n = s[1]
            if acc_n in acclist:
                acc2_alignment_data.setdefault(acc_n, [])
                acc2_alignment_data[acc_n].append(line)
                line = file_obj.readline()
                while line.startswith(">>") is False and line.startswith("Internal pipeline") is False:
                    acc2_alignment_data[acc_n].append(line)
                    line = file_obj.readline()
            else:
                line = file_obj.readline()
        elif line.startswith("Internal pipeline"):
            break
        else:
            line = file_obj.readline()
    return acc2_alignment_data


acc2folder = {}
for folder in os.listdir(check_folder):
    if os.path.isdir(os.path.join(check_folder, folder)):
        for name in os.listdir(os.path.join(check_folder, folder)):
            acc2folder.setdefault(name.rstrip(".align"), [])
            acc2folder[name.rstrip(".align")].append(folder)

acc_list = []
with open(tsv) as a:
    a.readline()
    for line in a:
        s = map(strip, line.split("\t"))
        acc = s[0].split("|")[0]
        if not acc in acc2folder:
            print "{} missing ACC".format(acc)
            acc_list.append(acc)
        elif len(acc2folder[acc]) != 2:
            print "{} missing alignment".format(acc), acc2folder[acc]
            acc_list.append(acc)

hmm_its1_tsv, alignment_folder, release, output_folder = "results_file_141.tsv", "align_141", "141", "etl"
hmm_out_folder = os.path.join(output_folder, "hmm_match_data")
if not os.path.exists(hmm_out_folder):
    os.mkdir(hmm_out_folder)
folder_58S = os.path.join(hmm_out_folder, "alignment_5_8S")
folder_SSU = os.path.join(hmm_out_folder, "alignment_18S")
if not os.path.exists(folder_58S):
    os.mkdir(folder_58S)
if not os.path.exists(folder_SSU):
    os.mkdir(folder_SSU)
diz18S = {}
diz5_8S = {}
print len(acc_list)
alignment_folder2acc_list = acc2align_folder(hmm_its1_tsv, release, alignment_folder, acc_list)
print len(alignment_folder2acc_list)
for folder, acc_list in alignment_folder2acc_list.items():
    if os.path.exists(folder):
        for name in os.listdir(folder):
            print os.path.join(folder, name)
            if name.startswith("align_18S") and name.endswith(".gz"):
                with gzip.open(os.path.join(folder, name)) as f_obj:
                    diz18S.update(parse_file(f_obj, acc_list))
            elif name.startswith("align_18S") and name.endswith(".txt"):
                with open(os.path.join(folder, name)) as f_obj:
                    diz18S.update(parse_file(f_obj, acc_list))
            if name.startswith("align_5_8S") and name.endswith(".gz"):
                with gzip.open(os.path.join(folder, name)) as f_obj:
                    diz5_8S.update(parse_file(f_obj, acc_list))
            elif name.startswith("align_5_8S") and name.endswith(".txt"):
                with open(os.path.join(folder, name)) as f_obj:
                    diz5_8S.update(parse_file(f_obj, acc_list))
print len(diz5_8S), len(diz18S)
for acc in diz5_8S:
    with open(os.path.join(folder_58S, "{}.align".format(acc.split("|")[0])), "w") as tmp:
        tmp.write("".join(diz5_8S[acc]))
    with open(os.path.join(folder_SSU, "{}.align".format(acc.split("|")[0])), "w") as tmp:
        tmp.write("".join(diz18S[acc]))
