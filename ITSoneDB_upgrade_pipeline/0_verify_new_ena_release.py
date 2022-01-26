#!/usr/bin/env python3
__author__ = "Giuseppe Defazio"
__version__ = 0
import os
import argparse
import argcomplete


def split_options():
    parser = argparse.ArgumentParser(
        description="It parses ENA release doc file and print the latest ENA release number",
        prefix_chars="-")
    parser.add_argument("-d", "--doc_release", type=str,
                        help="path to doc release file",
                        action="store", required=True)
    parser.add_argument("-r", "--releases_path", type=str,
                        help="path to releases directory",
                        action="store", required=True)
    argcomplete.autocomplete(parser)
    return parser.parse_args()


if __name__ == '__main__':
    opts = split_options()
    release_doc_path, releases_dir = opts.doc_release, opts.releases_path
    with open(release_doc_path, 'rt') as release_doc:
        line = release_doc.readline()
    s = line.split(" ")
    releases = []
    for n in os.listdir(releases_dir):
        try:
            rel = int(n)
            releases.append(rel)
        except ValueError:
            pass
    if int(s[1]) > max(releases):#valore della release precedente:
        print(s[1])
    else:
        print("false")
    # print("""
    # The latest version of ENA is %s.
    # You can try later for the upgrade!
    # See you later.
    # Bye.""")
    # sys.exit()
