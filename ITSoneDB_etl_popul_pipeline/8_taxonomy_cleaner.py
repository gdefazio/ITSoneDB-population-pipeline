#!/usr/bin/env python3
import os
import argparse
import argcomplete
import shutil
from librs.return_time import return_time


def input_options():
    parser = argparse.ArgumentParser(description="split fasta and tsv files", prefix_chars="-")
    parser.add_argument("-n", "--NCBI_tax_folder", type=str,
                        help="path to the NCBI taxonomy folder",
                        action="store", required=True)
    parser.add_argument("-l", "--its1_file", type=str,
                        help="tsv containg ITS1 location and annotation",
                        action="store", required=True)
    parser.add_argument("-o", "--output_folder", type=str,
                        help="output folder", action="store", required=False,
                        default="./etl/cleaned_ncbi_taxonomy")
    argcomplete.autocomplete(parser)
    return parser.parse_args()


def generazione_dizionari(ncbi_folder):
    node2parent = {}
    node2merged = {}
    with open(os.path.join(ncbi_folder, "nodes.dmp"), 'rt') as node:
        for linea in node:
            field = list(map(str.strip, linea.split("|")))
            node2parent[field[0]] = field[1]

    with open(os.path.join(ncbi_folder, "merged.dmp"), 'rt') as node:
        for linea in node:
            field = list(map(str.strip, linea.split("|")))
            node2merged[field[0]] = field[1]
    return node2parent, node2merged


def node_list(tsv, tax_folder):
    lista = set()
    node2parent, node2merged = generazione_dizionari(tax_folder)
    with open(tsv, 'rt') as a:
        a.readline()
        for linea in a:
            field = list(map(str.strip, linea.split("\t")))
            lista.add(field[2])
    lista_nodi = set()
    for node in lista:
        try:
            parent = node2parent[node]
        except KeyError:
            parent = node2parent[node2merged[node]]
        while node != parent:
            lista_nodi.add(node)
            node = parent
            try:
                parent = node2parent[node]
            except KeyError:
                parent = node2parent[node2merged[node]]
    lista_nodi.add("1")
    return lista_nodi


if __name__ == "__main__":
    return_time('Start taxonomy cleaner')
    param = input_options()
    ncbi_taxonomy, its1, out = param.NCBI_tax_folder, param.its1_file, param.output_folder
    new_taxonomy_dir = os.path.join(out, 'cleaned_ncbi_taxonomy')
    if not os.path.exists(new_taxonomy_dir):
        os.mkdir(new_taxonomy_dir)

    for name in os.listdir(ncbi_taxonomy):
        if name not in ["names.dmp", "nodes.dmp"] and name.endswith('.dmp'):
            shutil.copy(os.path.join(ncbi_taxonomy, name), new_taxonomy_dir)

    all_node = node_list(its1, ncbi_taxonomy)

    tmp = open(os.path.join(new_taxonomy_dir, "nodes.dmp"), "wt")
    with open(os.path.join(ncbi_taxonomy, "nodes.dmp")) as nodes:
        for line in nodes:
            s = list(map(str.strip, line.split("|")))
            if s[0] in all_node:
                tmp.write(line)
    tmp.close()

    tmp = open(os.path.join(new_taxonomy_dir, "names.dmp"), "wt")
    with open(os.path.join(ncbi_taxonomy, "names.dmp")) as nodes:
        for line in nodes:
            s = list(map(str.strip, line.split("|")))
            if s[0] in all_node:
                tmp.write(line)
    tmp.close()
    return_time('Stop taxonomy cleaner')