#!/usr/bin/env python3

# This script parses the dat files and create fasta and tsv dirs.
# Then stores fasta and tsv files
from librs.return_time import return_time
import os
import sys
#print(sys.version_info)
from string import ascii_letters, whitespace
import gzip
import multiprocessing as mp
import argparse
#import argcomplete

entry_sel_dataclass = ("gss", "htc", "htg", "mga", "wgs", "tsa", "sts", "std")
entry_sel_taxonomy = ("env", "fun", "hum", "inv", "mam", "vrt", "mus", "pln", "rod", "unc")
file_selection = []

feature_selection = ("misc_feature", "misc_RNA", "precursor_RNA", "rRNA", "source")
letters = list( ascii_letters )
whitespace = list( whitespace )


def split_options():
    parser = argparse.ArgumentParser(
        description="For each DAT file return a FASTA and TSV file",
        prefix_chars="-")
    parser.add_argument("-r", "--release", type=str,
                        help="the number of ena release",
                        action="store", required=True)
    parser.add_argument("-i", "--input", type=str,
                        help="folder of DAT file",
                        action="store", required=True
                        )
    parser.add_argument("-p", "--process", type=int,
                        help="number of max processes running",
                        action="store", required=True)
    # parser.add_argument("-o", "--out_file", type=str, help="output file path name",
    #                     action="store", required=False,
    #                     default="index/final_index.csv")
    #argcomplete.autocomplete(parser)
    return parser.parse_args()


def inner_splitter(open_file):
    """ determina l'accession
    :param open_file:
    :return:
    """
    acc2data = {}
    stringa = ""
    acc = ""
    for line in open_file:
        if line.startswith( "ID" ):
            linea = line.lstrip( "ID" )
            acc = list(map( str.strip, linea.split( ";" ) ))[0]
            acc2data.setdefault( acc, "" )
            stringa += line
        elif line.startswith( "//" ):
            acc2data[acc] = stringa
            stringa = ""
            acc = ""
        else:
            if len( line.strip( "\n" ).split() ) > 0:
                stringa += line
    return acc2data


def define_section_line(data_list):
    section2line = {}
    excl_line = []
    for section in ["ID", "AC", "AH", "AS", "DE", "KW", "DT", "OS", "OC", "RT", "RP", "RN", "RF", "RX", "RA", "RL", "DT", "DR", "FT", "FH", "SQ"]:
        # print section
        i = 0
        while i < len( data_list ):
            linea = data_list[i]
            if linea.startswith( section ):
                # print section, data_list[i]
                section2line[section] = i
                break
            i += 1
        i = 0
        while i < len( data_list ):
            if data_list[i].startswith( "XX" ):
                excl_line.append( i )
            i += 1
    # print section2line
    return section2line, excl_line


# per ogni file, scrivi accession, version, taxonomy, lunghezza seq estratta, feature table

def define_key_indent(linea):
    # definisce l'indentazione
    indent = 2
    linea = linea.strip( "\n" )
    for i in linea[2:]:
        if i in letters:
            break
        else:
            indent += 1
    return indent


def table_parser(item_list):
    """
    parse a ENA entries to find elements of interest
    :param item_list:
    :return: list of results
    """
    if len( item_list ) != 2:
        sys.exit( "the items list is longer than 2!!!" )
    else:
        acc = item_list[0]
        #print acc
        # data_file = item_list[1].split( "\n" )
        data_file = [i for i in item_list[1].split( "\n" ) if len( i ) >= 1]
        # return acc,len(data_file)
        #define_section_line toglie le linee che cominciano con XX
        section, exclusion_lines = define_section_line( data_file )
        #print data_file[0]
        lun = data_file[0].split()[-2]
        feature_data = {}
        species = ""
        taxonomy = []
        seq = ""
        version = ""
        description = []
        #####aggiunta g.d #######
        # mi porto dietro la sezione date
        date = []
        if "DT"  in section:
            for line in data_file[section["DT"]:section["DE"]]:
                s = list(map( str.strip, line.split() ))
                if s[0] == "DT":
                    date.append( " ".join(s[1:]) )
        date = "-->".join( date )
        #########################
        if "DE" in section:
            for line in data_file[section["DE"]:section["KW"]]:
                s = list(map( str.strip, line.split() ))
                if s[0] == "DE":
                    description.extend( s[1:] )
        description = " ".join( description )
        #print description
        # verifichiamo che la entry contenga il campo AC
        if "ID" in section: #nella versione precedente era AC e prendeva l'accession i
                                    #invece che il version
            version = list(map( str.strip, data_file[section["ID"]].split(";") ))[1].split()[1]
        for line in data_file[section["OS"]:section["RN"] - 1]:
            s = list(map( str.strip, line.split() ))
            if s[0] == "OS":
                species = " ".join( s[1:] )
            elif s[0] == "OC":
                taxonomy.extend( [a.strip( ";" ) for a in s[1:]] )
        else:
            pass
        # print taxonomy
        # verifichiamo che la entry contenga il campo taxonomy
        taxonomy[-1] = taxonomy[-1].strip( "." )
        taxonomy.append( species )
        taxonomy = "; ".join( taxonomy )
        indentazione = define_key_indent( data_file[section["FT"]] )
        # print indentazione
        key_placement = []
        i = section["FT"]
        while i < (section["SQ"] - 1):
            # print data_file[i]
            if i not in exclusion_lines:
                if data_file[i][indentazione] in letters:
                    key_placement.append( i )
            i += 1
        key_placement.append( section["SQ"] )
        # aggiungiamo il campo sequenza
        # print key_placement
        i = 0
        while i < len( key_placement ) - 1:
            # print i , i +2
            for linea in data_file[key_placement[i:i + 2][0]:key_placement[i:i + 2][1]]:
                if linea.startswith( "XX" ) is False:
                    # print line
                    if linea[indentazione] in letters:
                        # print linea
                        s = list( map( str.strip, linea.split() ) )
                        # print s
                        # individua la localizzazione
                        if s[1] in feature_selection:
                            chiave = "%s:location:%s" % (s[1], s[2])
                            #print chiave
                            feature_data.setdefault( chiave, [] )
                        else:
                            chiave = "pass"
                            # print chiave
                    else:
                        if chiave is not "pass":
                            s = list( map( str.strip, linea.split() ) )
                            if s[1].startswith( "/" ):
                                feature_data[chiave].append( " ".join( s[1:] ) )
                            else:
                                if len( feature_data[chiave] ) != 0:
                                    feature_data[chiave][-1] = feature_data[chiave][-1] + " ".join( s[1:] )
                                else:
                                    del feature_data[chiave]
                                    chiave += " ".join( s[1:] )
                                    feature_data.setdefault( chiave, [] )
            i += 1
        if data_file[section["SQ"]] == data_file[-1]:  # individua la sequenza
            seq = "NONE"
        else:
            i = section["SQ"] + 1
            while i < len( data_file ):
                s = list( map( str.strip, data_file[i].split() ) )
                if float( s[-1] ):
                    seq += "".join( s[:-1] )
                    i += 1
                else:
                    i += 1
        # print feature_data
        if seq is "NONE":
            return acc, version, date, description, taxonomy, "NONE", feature_data, seq
        else:
            return acc, version, date, description, taxonomy, len( seq ), feature_data, seq


def file_splitter(file_path):
    """
    Split GenBank multi-flatfile
    :param file_path: flat file GenBank
    :return: oggetto entry
    """
    if os.path.exists(file_path) is False:
        sys.exit( "the indicated flat file doesn't exists" )
    else:
        if file_path.endswith(".gz"):  # se il file termina con "gz", apri il file, e applica la funzione inner split
            input_data = gzip.open(file_path, "rt")
            acc2data = inner_splitter( input_data )
            input_data.close()
        else:  # altrimenti apri il file, chiamandolo input_data e poi inner split
            with open(file_path) as input_data:
                acc2data = inner_splitter(input_data)
        name = file_path.split("/")[-1].split('.')[0]
        with gzip.open(os.path.join(base, "tsv_%s" % release, "%s.tsv.gz" % name), "wt") as tmp, \
                gzip.open(os.path.join(base, "fasta_%s" % release, "%s.fasta.gz" % name), "wt") as tmp2:
            for items in acc2data.items():
                accession, version, date, description, taxonomy, \
                lenght_seq_extracted, features_dic, \
                sequence = table_parser( items )
                # print accession, length, version, taxonomy, locus_line, lenght_seq_extracted
                stringa = ["NEW ACC", accession, version, date,
                           description, taxonomy,
                           str( lenght_seq_extracted )]
                tax_id = ""
                for item in features_dic.items():
                    chiave, qualifier = item[0], " ### ".join( item[1] )
                    stringa.append( "%s\t%s" % (chiave, qualifier) )
                    #estrae il tax_id
                    for ft in item[1]:
                        if "taxon:" in ft:
                            tax_id = ft.split("taxon:")[1].strip('"')
                stringa.append(tax_id)
                stringa.append("--------------------------------------------------------")
                tmp.write("%s\n" % "\n".join(stringa))
                tmp2.write(">%s.%s|%s\n%s\n" % (
                    accession.strip(), version.strip(), tax_id, sequence))
    return 1


if __name__ == '__main__':
    return_time("""----------------------------------
    It splits each ENA dat file in Fasta and Tsv file""")
    param = split_options()
    processes, release, input_dir = param.process, param.release, param.input
    if input_dir.endswith('/'): input_dir = input_dir[:-1]
    base = '/'.join(input_dir.split("/")[:-1])
    return_time("Base release path %s" % base)
    # if input_dir.split("/")[-2] == str(release):
    try:
        os.mkdir(os.path.join(base, "fasta_%s" % release))
    except OSError:
        pass
    try:
        os.mkdir(os.path.join(base, "tsv_%s" % release))
    except OSError:
        pass
    # else:
    #     sys.exit("Please run this script in the %s release directory!" % release)
    list_file = []
    for file_name in os.listdir(input_dir):
        if file_name.startswith("rel_") and file_name.endswith(".dat.gz"):
            list_file.append(os.path.join(input_dir, file_name))
    # print(list_file)
    len_list_file = len(list_file)
    if len(list_file) < processes:
        processes = len(list_file)
        chunksize=1
    else:

        chunksize=len_list_file//processes
    with mp.Pool(processes=processes, initargs=[release]) as pool:
        check_list = pool.map(func=file_splitter, iterable=list_file,
                              chunksize=chunksize)
    len_check_list = check_list.count(1)
    if len_list_file == len_check_list:
        return_time('All the %s files have been processed' % len_check_list)
    else:
        raise Exception('ERROR: Only %s out of %s files have been processed' % (len_list_file, len_check_list))