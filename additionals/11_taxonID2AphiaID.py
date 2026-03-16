#!/usr/bin/env python3
__author__ = 'Bruno Fosso'
__version__ = "1.0"

import os
import argparse
import sys
import argcomplete
import requests
# import gzip
# import urllib
from urllib import parse
from time import sleep


def input_options():
    parser = argparse.ArgumentParser(description="split fasta and tsv files", prefix_chars="-")
    parser.add_argument("-t", "--its1_tsv", type=str,
                        help="TSV file containing data per ENA entry",
                        action="store", required=True)
    parser.add_argument("-n", "--ncbi_taxonomy_folder", type=str,
                        help="path to ncbi taxonomy folder",
                        action="store", required=True)
    parser.add_argument("-o", "--output_folder", type=str, help="output folder", action="store", required=False,
                        default=os.getcwd())
    # parser.add_argument("-a", "--alredy_done", type=str,
    #                     help="if the process was stopped restart by using the previous out file",
    #                     action="store", required=False, default=None)
    parser.add_argument("-p", "--previous_release", type=str,
                        help="tsv data for marine section from the previous release",
                        action="store", required=False, default=None)
    argcomplete.autocomplete(parser)
    return parser.parse_args()


def generazione_dizionari_tassonomia(taxnomy_folder, tax_list):
    # if os.path.exists(os.path.join(taxnomy_folder, "nodes.dmp")):
    #     nodesfile = open(os.path.join(taxnomy_folder, "nodes.dmp"))
    # else:
    #     sys.exit("No NODESFILE")
    if os.path.exists(os.path.join(taxnomy_folder, "names.dmp")):
        namesfile = open(os.path.join(taxnomy_folder, "names.dmp"))
    else:
        sys.exit("no namesfile")

    # parent_dict = {}
    # order_dict = {}
    # for linea in nodesfile:
    #     linea = linea.strip()
    #     fields = map(strip, linea.split("|"))
    #     nodeid, parentid = fields[0], fields[1]
    #     parent_dict[nodeid] = parentid
    #     order_dict[fields[0]] = fields[1]

    name_dict = {}
    for linea in namesfile:
        linea = linea.strip()
        fields = list(map(str.strip, linea.split("|")))
        nodeid, nome = fields[0], fields[1]
        if nodeid in tax_list:
            name_dict.setdefault(nodeid, [])
            name_dict[nodeid].append(nome)
    # return parent_dict, name_dict, order_dict
    return name_dict


def previous_release_data(tsv):
    d_list = {}
    with open(tsv, 'rt') as aa:
        for linea in aa:
            field = list(map(str.strip, linea.split("\t")))
            if field[-1] != "Error":
                d_list.setdefault(field[0], [])
                d_list[field[0]].append(field)
    return d_list


def taxid_list(its1_file):
    tax_list = set()
    with open(its1_file, 'rt') as aa:
        aa.readline()
        for linea in aa:
            tax_list.add(list(map(str.strip, linea.split("\t")))[2])
    return tax_list


def api_query(url_link):
    """
    API querying system
    :param url_link: API link
    :return: request result or NONE
    """
    c = 0
    answer = "error"
    # while c < 3:
        # print(c)
    try:
        risposta = requests.get(url_link, timeout=10)
        if risposta.content == "":
            answer = "none"
            c = 5
        elif int(risposta.content):
            answer = int(risposta.content)
            c = 5
    except ValueError:
        pass
            # c += 1
    return answer


if __name__ == "__main__":
    param = input_options()
    its1_tsv, ncbi_taxonomy, out, previous_release = param.its1_tsv, param.ncbi_taxonomy_folder, param.output_folder, param.previous_release
    out = os.path.join(out, 'marine_species')
    uncomplete_data = os.path.join(out, 'taxid2aphiaID.tsv')
    taxonID = taxid_list(its1_tsv)
    taxid2name = generazione_dizionari_tassonomia(ncbi_taxonomy, taxonID)

    # l = open(os.path.join(out, "taxid2name_list"), "w")
    # for taxid, name_list in taxid2name.items():
    #     l.write("%s\t (%s)\n" % (taxid," ; ".join(name_list)))
    # l.close()

    if not os.path.exists(out):
        os.mkdir(out)

    verified = set()
    try:
        outfile = uncomplete_data
        with open(uncomplete_data, 'rt') as a:
            a.readline()
            for line in a:
                s = list(map(str.strip, line.split("\t")))
                verified.add(s[0])
        tmp = open(uncomplete_data, "at")
    except FileNotFoundError:
        tmp = open(os.path.join(out, "taxid2aphiaID.tsv"), "wt")

    to_process = set(list(taxid2name.keys())).difference(verified)
    if previous_release is not None:
        previous_release_info = previous_release_data(previous_release)
        to_process = to_process.difference(set(previous_release_info.keys()))
        for node, data_list in previous_release_info.items():
            for ll in data_list:
                tmp.write("%s\n" % "\t".join(ll))
    print(len(to_process))
    for taxid in to_process:
        species_list = taxid2name[taxid]
        for name in species_list:
            # print(name)
            search_name = parse.quote(name)
            # print(search_name)
            url = "http://www.marinespecies.org/rest/AphiaIDByName/%s?marine_only=true" % search_name
            response = api_query(url)
            print(url, response)
            try:
                i = int(response)
                print(taxid, name, str(response))
                tmp.write("%s\t%s\t%s\n" % (taxid, name, str(response)))
            except (ValueError, requests.exceptions.ReadTimeout):
                search_name = parse.quote(name.replace("'", ""))
                url = "http://www.marinespecies.org/rest/AphiaIDByName/%s?marine_only=true" % search_name
                response = api_query(url)
                if response == "error":
                    tmp.write("%s\t%s\t%s\n" % (taxid, name, "Error"))
                elif response == "none":
                    print(taxid, name, "None")
                    tmp.write("%s\t%s\t%s\n" % (taxid, name, "None"))
                else:
                    i = int(response)
                    print(taxid, name, str(response))
                    tmp.write("%s\t%s\t%s\n" % (taxid, name, str(response)))
        sleep(10)
    tmp.close()
