__author__ = 'Bruno Fosso'
import getopt
import sys
import os
from string import find, strip
from Bio import SeqIO
import math


def usage():
    print """
    18S and 5.8S HMM profile mapping results are parsed in order to infer ITS1 boundaries.\n
    Usage:\n
         python estrazione_localizzazione_from_hmm.py  -f fasta_input -p path\n
    option:\n
         -h print this help\n
         -f FASTA sequences used as input for the HMM profiles mapping [MANDATORY]\n
         -p path for the working directory\n
    """


tab_18S = ""
tab_58S = ""
fasta_input = ""
directory = ""
try:
    opts, args = getopt.getopt( sys.argv[1:], "hf:p:" )
except getopt.GetoptError, err:
    print str( err )
    usage()
    sys.exit()
for o, a in opts:
    if o == "-h":
        usage()
        sys.exit()
    elif o == "-f":
        fasta_input = a
    elif o == "-p":
        directory = a
    elif o == "":
        usage()
        sys.exit()
    else:
        print "Unhandled option."
        usage()
        sys.exit()

if fasta_input != "":
    if os.path.exists( fasta_input ):
        pass
    else:
        print "The indicated FASTA file is inexistent"
        sys.exit()
else:
    print "The FASTA file option is missing!!!"
    sys.exit()

print """
---------------------------------------------------------------------------------------------
FASTA ACCESSION NUMBER EXTRACTION
\n
"""

lista_file = [fasta_input]
FASTA = SeqIO.index_db( "fungal_collection.idx", lista_file, "fasta" )
acc_list = set()
for acc in FASTA.keys():
    acc = acc.strip()
    acc_list.add( acc )
print

#################################################################################
#                                                                               #
# FASE DI ANALISI DEGLI OUTPUT PRODOTTI PER IL PROFILO RELATIVO ALL'rRNA 5.8S   #
#                                                                               #
#################################################################################
print """
THE SCRIPT CONSIDERS:\n
- THE E-VALUE OBTAINED FOR EACH ALIGNMENT\n
- COMPARISON OF THE ALGINED PORTION OF THE PROFILE COMPARED TO THE ALIGNED REGION OF THE SEQUENCE\n 
\n
\n
"""

# print "5.8S"
rna58 = {}

if os.path.exists( directory + "hmm_match_data/5_8S_alignment_file" ):
    pass
else:
    print "The alignment data for 5.8S HMM profile are missing"
    exit()

for alignment in os.listdir( directory + "hmm_match_data/5_8S_alignment_file" ):
    l = open( directory + "hmm_match_data/5_8S_alignment_file/" + alignment )
    lines = l.readlines()
    l.close()
    acc = alignment.rstrip( ".align" )
    tlen = len( FASTA[acc].seq )
    stop = len( lines )
    inizio = []
    fine = 0
    # print lines[0]
    index = 0
    print acc
    while index < stop:
        line = lines[index].strip()
        if find( line, "#    score" ) != -1:
            inizio.append( index + 1 )
        elif find( line, "Alignments for each domain:" ) != -1:
            fine = index
        index += 1
    # if acc == "EU343825.1":
    #    print inizio,fine
    #    print lines[min(inizio):fine]
    if len( inizio ) != 0 and fine != 0:
        # print lines[0]
        indice = min( inizio )
        while indice < fine:
            line = lines[indice].strip()
            s = map( strip, line.split( " " ) )
            data = []
            if len( s ) > 1:
                for i in s:
                    if i:
                        data.append( i )
            if len( data ) == 16:
                match = data[1]
                e_value = float( data[4] )
                pp = float( data[15] )
                if e_value <= 0.001 and tlen < 1000000:
                    rna58.setdefault( acc, [] )
                    rna58[acc].append( data )
                # elif e_value > 0.001 and pp >= 0.85 and tlen < 1000000:
                #     end_ali = int( data[10] )
                #     hmm_end = int( data[7] )
                #     if hmm_end <= 50 and end_ali == len( FASTA[acc].seq ):
                #         rna58.setdefault( acc, [] )
                #         rna58[acc].append( data )
            indice += 1

################################################################################
#                                                                              #
# FASE DI ANALISI DEGLI OUTPUT PRODOTTI PER IL PROFILO RELATIVO ALL'rRNA 18S   #
#                                                                              #
################################################################################

if os.path.exists( directory + "hmm_match_data/18S_alignment_file" ):
    pass
else:
    print "The alignment data for 18S HMM profile are missing"
    exit()

rna18 = {}
for alignment in os.listdir( directory + "hmm_match_data/18S_alignment_file" ):
    l = open( directory + "hmm_match_data/18S_alignment_file/" + alignment )
    lines = l.readlines()
    l.close()
    acc = alignment.rstrip( ".align" )
    tlen = len( FASTA[acc].seq )
    stop = len( lines )
    index = 0
    inizio = []
    fine = 0
    # print lines[0]
    while index < stop:
        line = lines[index].strip()
        if find( line, "#    score" ) != -1:
            inizio.append( index + 1 )
        elif find( line, "Alignments for each domain:" ) != -1:
            fine = index
        index += 1
    if len( inizio ) != 0 and fine != 0:
        # print lines[0]
        indice = min( inizio )
        while indice < fine:
            line = lines[indice].strip()
            s = map( strip, line.split( " " ) )
            data = []
            if len( s ) > 1:
                for i in s:
                    if i:
                        data.append( i )
            if len( data ) == 16:
                match = data[1]
                e_value = float( data[4] )
                pp = float( data[15] )
                if e_value <= 0.001 and tlen < 1000000:
                    rna18.setdefault( acc, [] )
                    rna18[acc].append( data )
                elif e_value > 0.001 and pp >= 0.85 and tlen < 1000000:
                    end_ali = int( data[10] )
                    hmm_end = int( data[7] )
                    hmm_start = int( data[6] )
                    if math.fabs( (hmm_end - hmm_start + 1) - end_ali ) <= 10:
                        rna18.setdefault( acc, [] )
                        rna18[acc].append( data )
                        # print acc
            indice += 1

print """
___________________________________________________________________________________________
___________________________________________________________________________________________
___________________________________________________________________________________________
ITS1 boundaries inferring
"""

common_acc = set( rna18.keys() ).intersection( set( rna58.keys() ) )

validated_rna58 = {}
for acc in common_acc:
    if len( rna58[acc] ) == 1:
        validated_rna58[acc] = rna58[acc]
    else:
        if len( rna58[acc] ) == 2:
            evalue = []
            mapped_hmm = []
            if (int( rna58[acc][0][10] ) - int( rna58[acc][1][9] )) <= 200:
                evalue.append( float( rna58[acc][0][4] ) )
                evalue.append( float( rna58[acc][1][4] ) )
                if float( rna58[acc][0][4] ) == min( evalue ) and float( rna58[acc][1][4] ) != min( evalue ):
                    validated_rna58.setdefault( acc, [] )
                    validated_rna58[acc].append( rna58[acc][0] )
                elif float( rna58[acc][0][4] ) != min( evalue ) and float( rna58[acc][1][4] ) == min( evalue ):
                    validated_rna58.setdefault( acc, [] )
                    validated_rna58[acc].append( rna58[acc][1] )
                elif float( rna58[acc][0][4] ) == min( evalue ) and float( rna58[acc][1][4] ) == min( evalue ):
                    mapped_hmm.append( 153 - (float( rna58[acc][0][7] ) - float( rna58[acc][0][6] ) + 1) )
                    mapped_hmm.append( 153 - (float( rna58[acc][1][7] ) - float( rna58[acc][1][6] ) + 1) )
                    if max( mapped_hmm ) == mapped_hmm[0]:
                        validated_rna58.setdefault( acc, [] )
                        validated_rna58[acc].append( rna58[acc][0] )
                    elif max( mapped_hmm ) == mapped_hmm[1]:
                        validated_rna58.setdefault( acc, [] )
                        validated_rna58[acc].append( rna58[acc][1] )
            else:
                len_match_1 = ((float( rna58[acc][0][10] ) - float( rna58[acc][0][9] ) + 1) / 153) * 100
                len_match_2 = ((float( rna58[acc][1][10] ) - float( rna58[acc][1][9] ) + 1) / 153) * 100
                if len_match_1 > 90 and len_match_2 > 90:
                    print >> tmp, "\t".join( rna58[acc][0] )
                    print >> tmp, "\t".join( rna58[acc][1] )
                elif len_match_1 > 90 and len_match_2 < 90:
                    validated_rna58.setdefault( acc, [] )
                    validated_rna58[acc].append( rna58[acc][0] )
                elif len_match_1 < 90 and len_match_2 > 90:
                    validated_rna58.setdefault( acc, [] )
                    validated_rna58[acc].append( rna58[acc][1] )

validated_rna18 = {}
for acc in common_acc:
    if len( rna18[acc] ) == 1:
        validated_rna18.setdefault( acc, [] )
        lista = [rna18[acc][0][9], rna18[acc][0][10], rna18[acc][0][4], rna18[acc][0][15]]
        validated_rna18[acc].append( lista )
    # elif len(rna18[acc]) == 2:
    #     #controlliamo che non sia un match spezzato
    #     if int(rna18[acc][1][7]) - int(rna18[acc][0][6]) <= 20:
    #         evalues = []
    #         evalues.append(rna18[acc][0][4])
    #         evalues.append(rna18[acc][1][4])
    #         pp_values = []
    #         pp_values.append(rna18[acc][0][15])
    #         pp_values.append(rna18[acc][1][15])
    #         validated_rna18.setdefault(acc,[])
    #         lista = [rna18[acc][0][6],rna18[acc][1][7],min(evalues),max(pp_values)]
    #         validated_rna18[acc].append(lista)
    elif len( rna18[acc] ) <= 3:
        evalues = []
        pp_values = []
        index = 0
        ali_start = rna18[acc][0][9]
        ali_end = rna18[acc][0][10]
        hmm_start = rna18[acc][0][6]
        hmm_end = rna18[acc][0][7]
        while index < len( rna18[acc] ) - 1:
            pp_values.append( rna18[acc][index][15] )
            evalues.append( rna18[acc][index][4] )
            if (int( rna18[acc][index + 1][6] ) - int( hmm_end )) <= 20:
                hmm_end = rna18[acc][index + 1][7]
                ali_end = rna18[acc][index + 1][10]
            index += 1
        pp_values.append( rna18[acc][index][15] )
        evalues.append( rna18[acc][index][4] )
        if hmm_end != rna18[acc][0][7]:
            validated_rna18.setdefault( acc, [] )
            lista = [ali_start, ali_end, min( evalues ), max( pp_values )]
            validated_rna18[acc].append( lista )
            # print acc
        else:
            parag = {}
            for dom in rna18[acc]:
                n_match = dom[0]
                hmm_start = dom[6]
                hmm_end = dom[7]
                ali_start = dom[9]
                ali_end = dom[10]
                c_value = dom[4]
                pp = dom[15]
                parag[n_match] = hmm_start + "\t" + hmm_end + "\t" + c_value + "\t" + pp + "\t" + ali_start + "\t" + ali_end
            ignora = []
            for dom in parag.keys():
                c = parag[dom].split( "\t" )
                hmm_start = int( c[0] )
                hmm_end = int( c[1] )
                c_value = float( c[2] )
                for counter in parag.keys():
                    if counter != dom:
                        s = parag[counter].split( "\t" )
                        hmm_start_conf = int( s[0] )
                        hmm_end_conf = int( s[1] )
                        c_value_conf = float( s[2] )
                        if math.fabs( hmm_start - hmm_start_conf ) <= 10 and math.fabs( hmm_end - hmm_end_conf ) <= 10:
                            if c_value_conf > c_value:
                                ignora.append( counter )
                            elif c_value_conf < c_value:
                                ignora.append( dom )
                            elif c_value_conf == c_value:
                                ignora.append( dom )
                                ignora.append( counter )
            sizes = []
            evalues = []
            pp_values = []
            for dom in parag.keys():
                if dom in ignora:
                    pass
                else:
                    s = parag[dom].split( "\t" )
                    sizes.append( int( s[-1] ) )
                    sizes.append( int( s[-2] ) )
                    evalues.append( float( s[2] ) )
                    pp_values.append( float( s[3] ) )
            if len( sizes ) != 0:
                start = min( sizes )
                end = max( sizes )
                e_value = min( evalues )
                pp = max( pp_values )
                validated_rna18.setdefault( acc, [] )
                lista = [str( start ), str( end ), str( e_value ), str( pp )]
                validated_rna18[acc].append( lista )
    else:
        parag = {}
        for dom in rna18[acc]:
            n_match = dom[0]
            hmm_start = dom[6]
            hmm_end = dom[7]
            ali_start = dom[9]
            ali_end = dom[10]
            c_value = dom[4]
            pp = dom[15]
            parag[n_match] = hmm_start + "\t" + hmm_end + "\t" + c_value + "\t" + pp + "\t" + ali_start + "\t" + ali_end
        ignora = []
        for dom in parag.keys():
            c = parag[dom].split( "\t" )
            hmm_start = int( c[0] )
            hmm_end = int( c[1] )
            c_value = float( c[2] )
            for counter in parag.keys():
                if counter != dom:
                    s = parag[counter].split( "\t" )
                    hmm_start_conf = int( s[0] )
                    hmm_end_conf = int( s[1] )
                    c_value_conf = float( s[2] )
                    if math.fabs( hmm_start - hmm_start_conf ) <= 10 and math.fabs( hmm_end - hmm_end_conf ) <= 10:
                        if c_value_conf > c_value:
                            ignora.append( counter )
                        elif c_value_conf < c_value:
                            ignora.append( dom )
                        elif c_value_conf == c_value:
                            ignora.append( dom )
                            ignora.append( counter )
        sizes = []
        evalues = []
        pp_values = []
        for dom in parag.keys():
            if dom in ignora:
                pass
            else:
                s = parag[dom].split( "\t" )
                sizes.append( int( s[-1] ) )
                sizes.append( int( s[-2] ) )
                evalues.append( float( s[2] ) )
                pp_values.append( float( s[3] ) )
        if len( sizes ) != 0:
            start = min( sizes )
            end = max( sizes )
            e_value = min( evalues )
            pp = max( pp_values )
            validated_rna18.setdefault( acc, [] )
            lista = [str( start ), str( end ), str( e_value ), str( pp )]
            validated_rna18[acc].append( lista )
            # print acc

result = open( "ITS1_boundaries.csv", "w" )
result.write(
    "accession \tlun \t start_18S \t end_18S \t e_value \t pp \t star_5.8S \t end_5.8S \t e_value \t pp\n" )

common_acc = set( validated_rna18.keys() ).intersection( set( validated_rna58.keys() ) )

for acc in common_acc:
    # print acc
    if len( validated_rna18[acc] ) > 1 and len( validated_rna58[acc] ) > 1:
        pass
    elif len( validated_rna18[acc] ) == 1 and len( validated_rna58[acc] ) == 1:
        # print acc
        s = (validated_rna18[acc][0])
        c = (validated_rna58[acc][0])
        tlen = str( len( FASTA[acc].seq ) )
        print >> result, acc + "\t" + tlen + "\t" + s[0] + "\t" + s[1] + "\t" + s[2] + "\t" + s[3] + "\t" + c[9] + "\t" + c[10] + "\t" + c[4] + "\t" + c[15]

result.close()
