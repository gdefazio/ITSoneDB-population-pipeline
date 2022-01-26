#!/usr/bin/env python3
import os
import argparse
import argcomplete
import gzip
from librs.return_time import return_time


description_info = """This script maps extracts the HMM alignment for ITS1 annotated by using the HMM.\n"""


def input_options():
    parser = argparse.ArgumentParser(description=description_info, prefix_chars="-")
    parser.add_argument("-m", "--hmm_results_file", type=str,
                        help="hmm result file",
                        action="store", required=True)
    parser.add_argument("-a", "--alignment_folder", type=str, help="alignment folder", action="store", required=True)
    parser.add_argument("-r", "--release", type=int, help="release number", action="store", required=True)
    parser.add_argument("-o", "--output_folder", type=str, help="output folder", action="store", required=False,
                        default="etl")
    argcomplete.autocomplete(parser)
    return parser.parse_args()


def acc2align_folder(tsv, release_n):
    """

    :param align_f: str
    :param tsv: str
    :param release_n: int
    :return :
    """
    folder2acc_list = {}
    with open(tsv) as a:
        a.readline()
        for line in a:
            s = list(map(str.strip, line.split("\t")))
            folder2acc_list.setdefault("rel_%s_r%s" % (s[-1], release_n), [])
            folder2acc_list["rel_%s_r%s" % (s[-1], release_n)].append(s[0])
    return folder2acc_list


def parse_file(file_obj, acclist):
    # print(acclist)
    acc2_alignment_data = {}
    line = file_obj.readline()
    while line:
        # print(line)
        if len(acc2_alignment_data) == len(acclist):
            # print(len(acc2_alignment_data), len(acclist))
            break
        elif line.startswith(">>"):
            s = list(map(str.strip, line.split(" ")))
            acc_n = s[1].split('.')[0]
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
    # print(acc2_alignment_data)
    return acc2_alignment_data


if __name__ == "__main__":
    info = input_options()
    hmm_its1_tsv, alignment_folder, \
        release, output_folder = info.hmm_results_file, \
                                 info.alignment_folder, info.release, info.output_folder
    hmm_out_folder = os.path.join(output_folder, "hmm_match_data")
    try:
        os.mkdir(hmm_out_folder)
    except FileExistsError:
        pass
    folder_58S = os.path.join(hmm_out_folder, "alignment_5_8S")
    folder_SSU = os.path.join(hmm_out_folder, "alignment_18S")
    try:
        os.mkdir(folder_58S)
    except FileExistsError:
        pass
    try:
        os.mkdir(folder_SSU)
    except FileExistsError:
        pass
    diz18S = {}
    diz5_8S = {}
    alignment_folder2acc_list = acc2align_folder(hmm_its1_tsv, release)
    for folder, acc_list in alignment_folder2acc_list.items():
        path = os.path.join(alignment_folder, folder)
        if os.path.exists(path):
            for name in os.listdir(path):
                # print(os.path.join(path, name))
                if name.startswith("align_18S") and name.endswith(".gz"):
                    with gzip.open(os.path.join(path, name), 'rt') as f_obj:
                        diz18S.update(parse_file(f_obj, acc_list))
                elif name.startswith("align_18S") and name.endswith(".txt"):
                    with open(os.path.join(path, name), 'rt') as f_obj:
                        diz18S.update(parse_file(f_obj, acc_list))
                if name.startswith("align_5_8S") and name.endswith(".gz"):
                    with gzip.open(os.path.join(path, name), 'rt') as f_obj:
                        diz5_8S.update(parse_file(f_obj, acc_list))
                elif name.startswith("align_5_8S") and name.endswith(".txt"):
                    with open(os.path.join(path, name), 'rt') as f_obj:
                        diz5_8S.update(parse_file(f_obj, acc_list))
    return_time("""-------------
    %s 5.8S mapping reports selected 
    %s 18S mapping reports selected """ % (len(diz5_8S), len(diz18S)))
    for acc in diz5_8S:
        with gzip.open(os.path.join(folder_58S, "{}.align.gz".format(acc.split("|")[0])), "wt") as tmp:
            tmp.write("".join(diz5_8S[acc]))
        with gzip.open(os.path.join(folder_SSU, "{}.align.gz".format(acc.split("|")[0])), "wt") as tmp:
            tmp.write("".join(diz18S[acc]))
    return_time("Mapping reports located at %s" % hmm_out_folder)
