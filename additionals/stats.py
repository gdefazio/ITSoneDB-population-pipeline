import argparse
import argcomplete
import os


def input_options():
    parser = argparse.ArgumentParser(description="split fasta and tsv files", prefix_chars="-")
    parser.add_argument("-r", "--result", type=str,
                        help="ITS1_loc_final_table_corrected.tsv in results dir of etl",
                        action="store", required=True)
    parser.add_argument("-t", "--taxonomy", type=str,
                        help="taxonomy location in etl",
                        action="store", required=True)
    # parser.add_argument("-n", "--ncbi_taxonomy_folder", type=str,
    #                     help="path to ncbi taxonomy folder",
    #                     action="store", required=True)
    # parser.add_argument("-o", "--output_folder", type=str, help="output folder", action="store", required=False,
    #                     default=os.getcwd())
    argcomplete.autocomplete(parser)
    return parser.parse_args()


def taxid_list4release(ITSoneDB_ITS1_summary_path):
    with open(ITSoneDB_ITS1_summary_path,'rt') as summary:
        summary.readline()
        acc2taxid = dict()
        for line in summary:
            s = line.strip().split('\t')
            acc = s[0]
            taxid = s[2]
            acc2taxid[acc] = taxid
    return acc2taxid


def import_lineage():
    lineage = dict()
    with open(os.path.join(taxonomy,'taxidlineage.dmp'), 'rt') as lin:
        for line in lin:
            s = list(map(str.strip, line.split('|')))
            #print(s)
            k = s[0]
            v = s[1].split(' ')
            #print(v)
            lineage[k] = v
    return lineage


def node2name():
    with open(os.path.join(taxonomy, "names.dmp"), "rt") as names:
        node2name = dict()
        # name2id = dict()
        line = names.readline()
        while line:
            s = list(map(str.strip, line.split("|")))
            line = names.readline()
            # ['100603', 'Peliosanthes sp. Tamura & Pooma 7031', '', 'equivalent name', '']
            if s[3] == 'scientific name':
                node2name[s[0]] = s[1]
        return node2name


def mergd():
    merd = dict()
    with open(os.path.join(taxonomy,'merged.dmp'), 'rt') as m:
        for line in m:
            s = list(map(str.strip, line.split('|')))
            #print(s)
            disappeared = s[0]
            newrepr = s[1]
            merd[disappeared] = newrepr
    return merd


if __name__ == '__main__':
    args = input_options()
    results, taxonomy = args.result, args.taxonomy
    counts = dict()
    merged = mergd()
    idlin = import_lineage()
    taxl = taxid_list4release(results)
    id2name = node2name()
    not_in_taxonomy = 0
    print('nr total species:', len(set([v for v in taxl.values()])), '\n')
    acc2lin = dict()
    complessivo = 0
    for acc, id in taxl.items():
        try:
            lineage = idlin[id]
            acc2lin[acc] = lineage
        except KeyError:
            try:
                newid = merged[id]
                lineage = idlin[newid]
                acc2lin[acc] = lineage
            except KeyError:
                print(acc, id, "not in taxonomy")
                #print(id)
                acc2lin[acc] = list()
                not_in_taxonomy += 1
    for lin in acc2lin.values():
        if len(lin) == 0:
            counts.setdefault('not_in_tax', 0)
            counts['not_in_tax'] += 1
        elif len(lin) == 1:
            #print(lineage)
            rank = id2name[lin[0]]
            counts.setdefault(rank, 0)
            counts[rank] += 1
        # se la lunghezza del lineage è maggiore di uno
        else: # if len(lin) > 1:
            # rank = id2name[lin[1]]
            if lin[1] == '2759': # eukaryota
                if lin[2] != '33154': # opistiokonta
                    rank = id2name[lin[2]]
                else:
                    rank = id2name[lin[3]]
                #print(rank)
                counts.setdefault(rank, 0)
                counts[rank] += 1
            else:
                rank = id2name[lin[1]]
                counts.setdefault(rank, 0)
                counts[rank] += 1

        complessivo += 1
    tot = 0
    print('\nGroups:')
    for k,v in counts.items():
        print(k, v)
        tot += v

    print('tot', tot)
    print('COMPLESSIVO',complessivo)