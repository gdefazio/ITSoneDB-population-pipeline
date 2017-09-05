import os
import sys
from string import strip, letters, whitespace
import gzip
import time

entry_sel_dataclass = ("gss", "htc", "htg", "mga", "wgs", "tsa", "sts", "std")
entry_sel_taxonomy = ("env", "fun", "hum", "inv", "mam", "vrt", "mus", "pln", "rod", "unc")
file_selection = []

feature_selection = ("misc_feature", "misc_RNA", "precursor_RNA", "rRNA", "source")
letters = list( letters )
whitespace = list( whitespace )


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
            acc = map( strip, linea.split( ";" ) )[0]
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


def file_splitter(l):
    """
    Split GenBank multi-flatfile
    :param l: flat file GenBank
    :return: oggetto entry
    """
    if os.path.exists( l ) is False:
        sys.exit( "the indicated flat file doesn't exists" )
    else:
        if l.endswith( "gz" ):  # se il file termina con "gz", apri il file, e applica la funzione inner split
            input_data = gzip.open( l )
            acc2data = inner_splitter( input_data )
            input_data.close()
        else:  # altrimenti apri il file, chiamandolo input_data e poi inner split
            with open( l ) as input_data:
                acc2data = inner_splitter( input_data )
        tmp = open( "%s_output.tsv" % l.split( "." )[0], "w" )
        tmp2 = open("%s.fasta" % l.split(".")[0], "w")
        for items in acc2data.items():
            accession, version, description, taxonomy, lenght_seq_extracted, features_dic, sequence = table_parser( items )

            # print accession, length, version, taxonomy, locus_line, lenght_seq_extracted
            stringa = ["NEW ACC", accession, version, description, taxonomy, str( lenght_seq_extracted )]
            for item in features_dic.items():
                chiave, qualifier = item[0], " ### ".join( item[1] )
                stringa.append( "%s\t%s" % (chiave, qualifier) )
            stringa.append( "--------------------------------------------------------" )
            tmp.write( "%s\n" % "\n".join( stringa ) )
            tmp2.write(">%s\n%s" % (accession.strip(), sequence))
        tmp.close()
        tmp2.close()


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
    parse a GenBank entries to find elements of interest
    :param item_list:
    :return: list of results
    """
    if len( item_list ) != 2:
        sys.exit( "the items list is longer than 2!!!" )
    else:
        acc = item_list[0]
        print acc
        # data_file = item_list[1].split( "\n" )
        data_file = [i for i in item_list[1].split( "\n" ) if len( i ) >= 1]
        # return acc,len(data_file)
        section, exclusion_lines = define_section_line( data_file )
        print data_file[0]
        lun = data_file[0].split()[-2]
        feature_data = {}
        species = ""
        taxonomy = []
        seq = ""
        version = ""
        description = []
        if section.has_key( "DE" ):
            for line in data_file[section["DE"]:section["KW"]]:
                s = map( strip, line.split() )
                if s[0] == "DE":
                    description.extend( s[1:] )
        description = " ".join( description )
        print description
        # verifichiamo che la entry contenga il campo AC
        if section.has_key( "AC" ):
            version = map( strip, data_file[section["AC"]].split() )[1].strip( ";" )
        for line in data_file[section["OS"]:section["RN"] - 1]:
            s = map( strip, line.split() )
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
                        s = map( strip, linea.split() )
                        # print s
                        # individua la localizzazione
                        if s[1] in feature_selection:
                            chiave = "%s:location:%s" % (s[1], s[2])
                            print chiave
                            feature_data.setdefault( chiave, [] )
                        else:
                            chiave = "pass"
                            # print chiave
                    else:
                        if chiave is not "pass":
                            s = map( strip, linea.split() )
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
                s = map( strip, data_file[i].split() )
                if float( s[-1] ):
                    seq += "".join( s[:-1] )
                    i += 1
                else:
                    i += 1
        # print feature_data
        if seq is "NONE":
            return acc, version, description, taxonomy, "NONE", feature_data, seq
        else:
            return acc, version, description, taxonomy, len( seq ), feature_data, seq


start = time.localtime()

for file_name in os.listdir( "." ):
    print file_name
    if file_name.endswith( "dat.gz" ) and os.path.exists( "%s_output.tsv" % file_name.split( "." )[0] ) is False:
        s = map( strip, file_name.split( "_" ) )
        if s[1] in entry_sel_dataclass and s[2] in entry_sel_taxonomy:
            s = map( strip, file_name.split( "_" ) )
            file_splitter( file_name )

end = time.localtime()
print start
print end
