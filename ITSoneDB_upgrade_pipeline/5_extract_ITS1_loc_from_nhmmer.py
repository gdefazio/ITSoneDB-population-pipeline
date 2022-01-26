#!/usr/bin/env python3
__author__ = "Giuseppe Defazio"
__version__ = 2

# versione utilizzata per analizzare gli autput ottenuti da nhmmer

from librs.return_time import return_time
import os
import sys
import gzip
from math import fabs
import getopt
# import shlex
# import Bio
# from Bio import SeqIO
import multiprocessing as mp
from shutil import rmtree
from numpy import mean


def usage():
    "Usage"
    print("""
    This program has a folder path in input.\n
    In this folder there are the folders containing $nhmmer output files for\n
    18S and 5.8S rRNA alignment results and the related FASTA file.\n
    The program infers the ITS1 boundaries with related metadata stored in a .csv file\n
    USAGE:\n
    $ nohup estraz_localiz_nhmmer.py -f align_path -s fasta_path -p 10 -r 141 -o results_file.csv > log_file.log &\n
    \n
    \n
    """)


def pars_arg():
    """
    The pars_arg() function parses the indicated arguments.
    The arguments are the folder that contains the folders in which there
    are gzipped fasta file with the output files of $nhmmer and the name of ITS1 output
    csv file.
    :return:the folder path string, the name of output file.
    """
    path_folder_of_folders = ""
    path_fasta_folder = ""
    outname = "ITS1_boundaries.csv"
    processes = ""
    release = ""
    try:
        getopt.getopt(sys.argv[1:], "hf:s:o:p:r:")
    except getopt.GetoptError:
        usage()
        sys.exit()
    opt, arg = getopt.getopt(sys.argv[1:], "hf:s:o:p:r:")
    if len(opt) == 0:
        usage()
        sys.exit()
    for o, a in opt:
        if o == "-h":
            usage()
            sys.exit()
        elif o == "-f":
            path_folder_of_folders = a
        elif o == "-s":
            path_fasta_folder = a
        elif o == "-o":
            outname = a
        elif o == "-p":
            processes = a
        elif o == "-r":
            release = a
        else:
            usage()
            sys.exit()
    return path_folder_of_folders, path_fasta_folder, outname, processes, release


def data_extr_valid_18(acc_curs_dic, lines_list): #read_fasta):
    """For each common accession, this function extracts data accession related
    about the alignment between the object sequence and 18S-HMM profile.
    If there is only one alignment, this function verifies that the alignment
    passes through the filter and adds the data of alignment in a tupla that
    represent the value object in which the key is the accession number.
    If there are two alignments, this function tries to reconstruct a
    single alignment and store the data into the dictionary. If there are more than
    2 alignments for a single sequence, the sequence is discarded.
    than 2 domains are discarded.
    :param acc_curs_dic:dictionary with the accession numbers as keys and the associated
    row number in the file, to retrieve the information.
    :param lines_list:list containing the rows of parsed document.
    :return:pre_result: a dictionary with the accession numbers as keys and a
    sub-dictionary as value. In the sub-dictionary is/are stored the validated alignment/s.
    """
    #estrazione
    #sorto le chiavi indicanti il numero di riga in cui trovo l'header >> per ogni accession
    #comune presente nel documento. Se ad un accession corrispondono piu' di due domini
    #cio' che e' stato acquisito viene cancellato dalla memoria. I dati acquisiti vengono
    #messi in un dizionario temporaneo che al n di dominio associa una lista di dati.
    #Dopo di che doms4acc viene analizzato.
    #print "start data_estr_valid_18()"
    pre_results = {}
    for acc in acc_curs_dic.keys():
        annote = ""
        doms4acc = {}
        cursors_list = acc_curs_dic[acc]
        for n in range(0, len(cursors_list)):
            line = lines_list[cursors_list[n] + 3]
            data = []
            s = list(map(str.strip, line.split(" ")))
            for val in s:
                if val:
                    try:
                        data.append(float(val))
                    except ValueError:
                        data.append(val)
            if len(data) == 15:
                doms4acc[n] = data
        numbers_of_domains = len(doms4acc)
        if numbers_of_domains == 1:
            #se l'accession presenta un solo dominio di allineamento hmm si impone la condizione che
            #l'evalue sia inferiore a 1e-3 oppure superiore ma con una posterior probability maggiore di
            #uguale di 0.85.
            evalue = doms4acc[0][3]
            hmmfrom = doms4acc[0][4]
            hmmto = doms4acc[0][5]
            alignfrom = doms4acc[0][7]
            alignto = doms4acc[0][8]
            seqlen = doms4acc[0][13]
            pstprob = doms4acc[0][14]
            rigth_end = (1800 - hmmto) + 1
            hmm_coverage = round(((fabs((hmmto - hmmfrom) + 1)) / 1800.0), 2)
            #annote = ""
            if (evalue <= 0.001) or ((evalue > 0.001) and (pstprob >= 0.85)):
                #verifico che la distanza del dominio dalla
                #ultima posizione dell'hmm sia <= a 20 nt.
                if rigth_end <= 20:
                    #annoto la direzione dell'allineamento
                    if (alignto - alignfrom) > 0:
                        annote = "+"
                    else:
                        annote = "-"
                    pre_results[acc] = {}
                    pre_results[acc][0] = (evalue, hmm_coverage, pstprob,
                        alignfrom, alignto, seqlen, annote)
        if numbers_of_domains == 2:
            #se i domini sono due provo ad unirli, ma per farlo verifico che l'orientamento
            #dell'allineamento dei due domini sia concorde. Se discorde vengono scartati
            #####ordino i domini in ordine posizionale#######
            m_dict = {}
            for k, v in doms4acc.items():
                tpl = (v[7], v[8])
                m = mean(tpl)
                m_dict[m] = v
            doms4acc[0] = m_dict[min(m_dict.keys())]
            doms4acc[1] = m_dict[max(m_dict.keys())]
            #################################################
            evalue_0 = doms4acc[0][3]
            hmmfrom_0 = doms4acc[0][4]
            hmmto_0 = doms4acc[0][5]
            hmm_coverage_0 = round((fabs((hmmto_0 - hmmfrom_0) + 1) / 1800.0), 2)
            alignfrom_0 = doms4acc[0][7]
            alignto_0 = doms4acc[0][8]
            seqlen_0 = doms4acc[0][13]
            pstprob_0 = doms4acc[0][14]
            rigth_end_0 = (1800 - hmmto_0) + 1
            evalue_1 = doms4acc[1][3]
            hmmfrom_1 = doms4acc[1][4]
            hmmto_1 = doms4acc[1][5]
            hmm_coverage_1 = round((fabs((hmmto_1 - hmmfrom_1) + 1) / 1800.0), 2)
            alignfrom_1 = doms4acc[1][7]
            alignto_1 = doms4acc[1][8]
            seqlen_1 = doms4acc[1][13]
            pstprob_1 = doms4acc[1][14]
            rigth_end_1 = (1800 - hmmto_1) + 1
            orient_0 = alignto_0 - alignfrom_0
            orient_1 = alignto_1 - alignfrom_1
            #l'orientamento dell'allineamento pregiudica l'ordine con cui si procede alla
            # validazione. In caso di orientamento positivo, procedo alla validazione dal
            #secondo al primo dominio.
            pre_results[acc] = {}
            if (orient_0 > 0) and (orient_1 > 0):
                #orientazione concorde e positiva
                annote = "+"
                #########################################
                #VERIFICO SE ENTRAMBE I DOMINI SONO VALIDI E MARGINI
                #verifico che il secondo dominio abbia il margine destro validato
                #e che il suo coverage non sia inferiore a 0.85
                #validazione di evalue di entrambe gli allineamenti
                if ((evalue_1 <= 0.001) or ((evalue_1 > 0.001) and (pstprob_1 >= 0.85))) and \
                        ((evalue_0 <= 0.001) or ((evalue_0 > 0.001) and (pstprob_0 >= 0.85))):
                    #validazione di margine di entrambe gli allineamenti
                    if (rigth_end_1 <= 20) and (rigth_end_0 <= 20):
                        #se entrambe gli allineamenti hanno un coverage di 0.85 li valido
                        if (hmm_coverage_0 >= 0.85) and (hmm_coverage_1 >= 0.85):
                            pre_results[acc][1] = (
                                    evalue_1, hmm_coverage_1, pstprob_1,
                                    alignfrom_1, alignto_1, seqlen_1, annote)
                            pre_results[acc][0] = (
                                evalue_0, hmm_coverage_0, pstprob_0,
                                alignfrom_0, alignto_0, seqlen_0, annote)
                        if (hmm_coverage_0 < 0.85) and (hmm_coverage_1 >= 0.85):
                            #se il primo allineamento ha un coverage minore di 0.85
                            #e il secondo maggiore, valido il secondo.
                            pre_results[acc][1] = (
                                evalue_1, hmm_coverage_1, pstprob_1,
                                alignfrom_1, alignto_1, seqlen_1, annote)
                            coverage_ponderato = \
                                ((hmmto_0 - hmmfrom_0) + 1)/float(alignto_0)
                            if coverage_ponderato >= 0.85:
                                #per il primo verifico che il coverage della regione iniziale
                                #di sequenza coperta dall'allineamento sia maggiore uguale
                                #di 0.85
                                pre_results[acc][0] = (
                                    evalue_0, hmm_coverage_0, pstprob_0,
                                    alignfrom_0, alignto_0, seqlen_0, annote)
                        if (hmm_coverage_0 < 0.85) and (hmm_coverage_1 < 0.85):
                            #se per entrambe gli allineamenti il coverage e' inferiore a 0.85
                            #scelgo tra i due domini quello piu' coerente biologicamente ovvero
                            #quello posto ad inizio sequenza.
                            annote = "+other_align(%i..%i)" % (alignfrom_1, alignto_1)
                            pre_results[acc][0] = (
                                evalue_0, hmm_coverage_0, pstprob_0,
                                alignfrom_0, alignto_0, seqlen_0, annote)
                            #print "18s", acc, pre_results[acc]
                    #se il primo dominio non e' margine provo a fondererlo con il secondo
                    if (rigth_end_1 <= 20) and (rigth_end_0 > 20):
                        hmmdistance = ((hmmfrom_1 - hmmto_0) + 1)
                        #verifico le distanze dei domini su hmm
                        aligndistance = ((alignfrom_1 - alignto_0) + 1)
                        if hmmdistance <= 0:
                            global_coverage = hmm_coverage_0 + hmm_coverage_1 + (hmmdistance/1800.0)
                        elif hmmdistance > 0:
                            global_coverage = hmm_coverage_0 + hmm_coverage_1
                        #coverage_ponderato = (global_coverage * 1800) / float(alignto_1)
                        #print coverage_ponderato, acc
                        if (global_coverage >= 0.85) and (global_coverage <= 1) and \
                                (mean([int(hmmfrom_0), int(hmmto_0)]) <
                                 mean([int(hmmfrom_1), int(hmmto_1)])):
                            #il coverage dei due allineamenti fusi deve essere compreso tra 0.85 e 1.
                            annote = "+join(%i...%i,%i...%i)" % (
                                alignfrom_0, alignto_0, alignfrom_1, alignto_1)
                            #print annote, acc
                            evalue = max(evalue_0, evalue_1)
                            pstprob = min(pstprob_0, pstprob_1)
                            pre_results[acc][1] = (
                                evalue, global_coverage, pstprob,
                                alignfrom_0, alignto_1, seqlen_0, annote)
                    if (rigth_end_1 > 20) and (rigth_end_0 <= 20):
                        #print acc, "uno e mezzo0"
                        #se il secondo dominio non e' margine sinistro provo a vedere
                        #se  margine destro. L'info potra' essermi utile in fase di valizazione
                        if (hmmfrom_1 <= 20) and (hmmfrom_0 <= 20):
                            #print acc, "uno e mezzo1"
                            pre_results[acc][1] = [alignfrom_1]
                            pre_results[acc][0] = (
                                evalue_0, hmm_coverage_0, pstprob_0,
                                alignfrom_0, alignto_0, seqlen_0, annote)

                #lo stesso algoritmo decisionale viene utilizzato ed adattato al caso in cui
                #l'orientamento degli allineamenti fosse negativo.
            ########################################################################
            ###########            N E G A T I V O                 #################
            ########################################################################
            if (orient_0 < 0) and (orient_1 < 0):#orientazione domini hmm negativa
                #print acc_curs_dic[cursor], "YES"
                annote = "-"
                if ((evalue_1 <= 0.001) or ((evalue_1 > 0.001) and (pstprob_1 >= 0.85))) and \
                        ((evalue_0 <= 0.001) or ((evalue_0 > 0.001) and (pstprob_0 >= 0.85))):
                    #sono entrambe margine destro
                    if (rigth_end_1 <= 20) and (rigth_end_0 <= 20):#il coverage di entrambe gli allieamenti e' superiore all'85%
                        #il coverage di entrambe gli allieamenti e' superiore all'85%
                        if (hmm_coverage_1 >= 0.85) and (hmm_coverage_0 >= 0.85):
                            pre_results[acc][1] = (
                                evalue_1, hmm_coverage_1, pstprob_1,
                                alignfrom_1, alignto_1, seqlen_1, annote)
                            pre_results[acc][0] = (
                                evalue_0, hmm_coverage_0, pstprob_0,
                                alignfrom_0, alignto_0, seqlen_0, annote)
                        # il coverage di entrambe gli allieamenti e' inferiore all'85%
                        if (hmm_coverage_1 < 0.85) and (hmm_coverage_0 < 0.85):
                            annote = "-other_align(%i..%i)" % (alignfrom_0, alignto_0)
                            pre_results[acc][1] = (
                                evalue_1, hmm_coverage_1, pstprob_1,
                                alignfrom_1, alignto_1, seqlen_1, annote)
                        if (hmm_coverage_1 < 0.85) and (hmm_coverage_0 >= 0.85):
                            pre_results[acc][0] = (
                                evalue_0, hmm_coverage_0, pstprob_0,
                                alignfrom_0, alignto_0, seqlen_0, annote)
                            coverage_ponderato = \
                                ((hmmto_1 - hmmfrom_1) + 1) / float(alignto_1)
                            if coverage_ponderato >= 0.85:
                                pre_results[acc][1] = (
                                    evalue_1, hmm_coverage_1, pstprob_1,
                                    alignfrom_1, alignto_1, seqlen_1, annote)
                    if (rigth_end_1 > 20) and (rigth_end_0 <= 20):
                        hmmdistance = ((hmmfrom_0 - hmmto_1) + 1)
                        #verifico le distanze dei domini su hmm
                        aligndistance = ((alignfrom_1 - alignto_0) + 1)
                        if hmmdistance <= 0:
                            global_coverage = hmm_coverage_0 + hmm_coverage_1 + (hmmdistance/1800.0)
                        if hmmdistance > 0:
                            global_coverage = hmm_coverage_0 + hmm_coverage_1
                        #coverage_ponderato = global_coverage * 1800 / float(alignfrom_1 - alignto_0)
                        #print hmm_coverage_0, hmm_coverage_1, coverage_ponderato, global_coverage, acc, hmmdistance/1800.0, " REV"
                        if (global_coverage >= 0.85) and (global_coverage <= 1) and \
                                (mean([int(hmmfrom_0), int(hmmto_0)]) > mean([int(hmmfrom_1), int(hmmto_1)])):
                            #aggiungi un passaggio in cui verifichi le N nell'introne
                            annote = "-join(%i...%i,%i...%i)" % (
                                alignfrom_1, alignto_1, alignfrom_0, alignto_0)
                            evalue = max(evalue_0, evalue_1)
                            pstprob = min(pstprob_0, pstprob_1)
                            pre_results[acc][1] = (
                                evalue, global_coverage, pstprob,
                                alignfrom_1, alignto_0, seqlen_0, annote)
                            #print acc, pre_results[acc][1], "REV"
                    if (rigth_end_1 <= 20) and (rigth_end_0 > 20):
                        # se il secondo dominio non e' margine sinistro provo a vedere
                        # se  margine destro. L'info potra' essermi utile in fase di valizazione
                        if (hmmfrom_1 <= 20) and (hmmfrom_0 <= 20):
                            pre_results[acc][0] = [alignfrom_0]
                            pre_results[acc][1] = (
                                evalue_1, hmm_coverage_1, pstprob_1,
                                alignfrom_1, alignto_1, seqlen_1, annote)
    #print "stop data_estr_valid_18()"
    return pre_results


def data_extr_valid_58(acc_curs_dic, lines_list):
    """For each common accession, this function extracts the data accession related
    about the alignment between query sequence and 5.8S HMM profile. If there is only
    one alignment domain, this function verifies that the aligment passes through the
    filter and adds the data of validated alignment in a tupla. The tupla is stored as
    value, in a dictionary in which the key is the accession number. The alignments with
    more than 1 domain are discarded.
    :param acc_curs_dic: dictionary with the accession numbers as keys and the associated
    row number in the file, to retrieve the information.
    :param lines_list:list containing the rows of parsed document.
    :return:pre_result: a dictionary with the accession numbers as keys and a
    tupla as value. In the tupla are stored the information about the validated alignments.
    """
    #estrazione
    # sorto le chiavi indicanti il numero di riga in cui trovo l'header >> per ogni accession
    # comune presente nel documento. Se ad un accession corrispondono piu' di un dominio
    # cio' che e' stato acquisito viene cancellato dalla memoria. I dati acquisiti vengono messi
    # in un dizionario remporaneo che al n di dominio associa una lista di dati. Dopo di che
    # doms4acc viene analizzato.
    #print "start data_estr_valid_58()"
    pre_results = {}
    ###################################################
    for acc in acc_curs_dic.keys():                   #
        doms4acc = {}                                 #
        cursors_list = acc_curs_dic[acc]              #
        for n in range(0, len(cursors_list)):        #
            line = lines_list[cursors_list[n] + 3]    #
            data = []                                 #
            s = list(map(str.strip, line.split(" ")))           #Estrazione informazioni
            for val in s:                             #
                if val:                               #
                    try:                              #
                        data.append(float(val))       #
                    except ValueError:
                        data.append(val)
            if len(data) == 15:
                doms4acc[n] = data
    ####################################################
        numbers_of_domains = len(doms4acc)
        #print doms4acc
        if numbers_of_domains == 1:
            #print numbers_of_domains
            #Analizzo solo gli allineamenti con un dominio
            evalue = doms4acc[0][3]
            hmmfrom = doms4acc[0][4]
            hmmto = doms4acc[0][5]
            alignfrom = doms4acc[0][7]
            alignto = doms4acc[0][8]
            seqlen = doms4acc[0][13]
            pstprob = doms4acc[0][14]
            hmm_coverage = round(((fabs((hmmto - hmmfrom) + 1)) / 154.0), 2)
            #verifico che il dominio passi i filtri di evalue e posterior probability
            if evalue <= 0.001 or ((evalue > 0.001) and (pstprob >= 0.85)):
                #verifico che la distanza del dominio dall'inizio
                #dell'HMM 5.8 sia non superione a 20nt
                if hmmfrom <= 20:
                    pre_results[acc] = (
                            evalue, hmm_coverage, pstprob,
                            alignfrom, alignto)
    #print "stop data_estr_valid_58()"
    return pre_results


def assign(coordinate_18, coordinate_58):
    """
    For each sequence, from information about 18S and 5.8S alignments coordinates, this
    function assigns the ITS1 coordinates. ITS1 coordinates are then stored into a tupla.
    :param coordinate_18:
    :param coordinate_58:
    :return: its1_from_to: a tupla that contains in order: 18S e-value, 18S HMM coverage,
    18S posterior probability, 18S align from, 18S align to, 5.8S e-value,
    5.8S HMM coverage, 5.8S posterior probability, 5.8S align from, 5.8S align to,
    length of sequence, annotations(e.g. strand polarity)
    """
    its1_from_to = (coordinate_18[0], coordinate_18[1], coordinate_18[2],
                    int(coordinate_18[3]), int(coordinate_18[4]),
                    coordinate_58[0], coordinate_58[1], coordinate_58[2],
                    int(coordinate_58[3]), int(coordinate_58[4]),
                    (coordinate_18[-2]), coordinate_18[-1])
    return its1_from_to


def validation_concordance(rslts18, rslts58, sample_name):
    """
    This function try to find a biological coherence between the data of 18S and 5.8S
    alignments for each sequence. Then the ITS1 are validated.
    :param rslts18:a dictionary with the accession numbers as keys and a
    sub-dictionary as value. In the sub-dictionary is/are stored the validated alignment/s.
    :param rslts58:a dictionary with the accession numbers as keys and a
    tupla as value. In the tupla are stored the information about the validated alignments
    :param sample_name:
    :return:validated_its_boundaries: a dictionary with the accession numbers as keys
    and the tuplas with ITS1 boundaries as value.
    """
    #print "start validation_concordance()"
    validated_its_boundaries = {}
    #its1_from_to = ()
    #verifico tra gli accession che hanno passato la fase di
    # validazione quali sono quelli comuni
    common_accession = set(rslts18.keys()).intersection(set(rslts58.keys()))
    return_time("""- - - -
    [%s]>> ACCESSION VALIDATION
    [%s]accession_18: %s accession_58: %s common_accession: %s\n""" % (
       sample_name, sample_name, len(rslts18.keys()), len(rslts58.keys()),
       len(common_accession)))
    #verifico che gli allineamenti di 18S e 5.8S validati siano concordi tra loro,
    #in caso negativo vengono scartati. In caso positivo i dati vengono immagazzinati
    # in una tupla che viene associata all'accession che rappresenta la chiave
    # del dizionario validated_its_boundaries.
    for acc in common_accession:
        orient58 = rslts58[acc][4] - rslts58[acc][3]
        coordinate_58 = rslts58[acc]
        #print rslts18[acc]
        if len(rslts18[acc].items()) == 1:
        #se l'allineamento validato e' solo uno
            coordinate_18 = list(rslts18[acc].values())[0]
            #print rslts18[acc], acc
            #print coordinate_18, acc
            orient18 = coordinate_18[4] - coordinate_18[3]
            if orient18 > 0 and orient58 > 0:
            #se l'orientamento e' concorde e positivo tra i due geni fiancheggianti
                if coordinate_18[4] < coordinate_58[3]:
                #se 18S to e' minore di 5.8S from
                    #validiamo
                    validated_its_boundaries[acc] = assign(coordinate_18, coordinate_58)
                    #print "1 dom 18s, 1 dom 5.8 concordanza +", acc
            if orient58 < 0 and orient18 < 0:
            #se l'orientamento e' concorde e negativo tra i due geni fiancheggianti
                if coordinate_18[4] > coordinate_58[3]:
                #se 18s to e' maggiore di 5.8S from
                    #validiamo
                    validated_its_boundaries[acc] = assign(coordinate_18, coordinate_58)
                    #print "1 dom 18s, 1 dom 5.8 concordanza -", acc
        elif len(rslts18[acc].items()) == 2:
        #se i geni 18S individuati sono 2
            if len(rslts18[acc][0]) == 1 and orient58 < 0:
            #se abbiamo la coordinata iniziale del primo allineamento di 18S
            #se l'orientazione di 5.8S e' negativa
                coordinate_18 = rslts18[acc][1]
                orient18 = coordinate_18[4] - coordinate_18[3]
                if (orient18 < 0) and (rslts18[acc][0][0] < coordinate_58[4]):
                    #se orientamento concorde
                    # se alifrom58s maggiore di alito18s
                    if coordinate_18[4] > coordinate_58[3]:
                        validated_its_boundaries[acc] = assign(coordinate_18, coordinate_58)
                        #print "solo margine sinistro primo dominio -", acc
            elif len(rslts18[acc][1]) == 1 and orient58 > 0:
            #se alifrom18s dell'operone seguente e' maggiore di
            #alito58s dell'operone precedente
                coordinate_18 = rslts18[acc][0]
                orient18 = coordinate_18[4] - coordinate_18[3]
                if orient18 >0 and (rslts18[acc][1][0] > coordinate_58[4]):
                # orientazione positiva
                #se alifrom58s maggiore di alito18s
                    if coordinate_18[4] < coordinate_58[3]:
                        validated_its_boundaries[acc] = assign(coordinate_18, coordinate_58)
                        #print "solo margine sinistro secondo dominio +", acc
            else:
                #entrambe i domini sono validi
                coordinate_18_0 = rslts18[acc][0]
                coordinate_18_1 = rslts18[acc][1]
                orient18 = coordinate_18_0[4] - coordinate_18_0[3]
                if orient18 > 0 and orient58 > 0:
                    #verifico a quale dei due 18s appartiene il 58s
                    if (coordinate_18_0[4] < coordinate_58[3]) and \
                            (coordinate_58[4] < coordinate_18_1[3]):
                        #il 5.8 si trova tra i due 18s
                        coordinate_18 = coordinate_18_0
                        validated_its_boundaries[acc] = assign(coordinate_18, coordinate_58)
                        #print "2 18s validati, 5.8 appartiene al 1 dom +", acc
                    if (coordinate_18_0[4] < coordinate_18_1[4]) and \
                            (coordinate_18_1[4] < coordinate_58[3]):
                        #il 5.8 appartiene al secondo 18s
                        coordinate_18 = coordinate_18_1
                        validated_its_boundaries[acc] = assign(coordinate_18, coordinate_58)
                        #print "2 18s validati, 5.8 appartiene al 2 dom +", acc
                if orient58 < 0 and orient18 < 0:
                    # verifico a quale dei due 18s appartiene il 58s
                    if (coordinate_18_0[4] < coordinate_58[4]) and \
                            (coordinate_58[4] < coordinate_18_1[4]):
                        #print acc
                        #print "+++", coordinate_18_0[4], rslts58[acc][4], coordinate_18_1[4]
                        # il 5.8 si trova tra i due 18s
                        coordinate_18 = coordinate_18_1
                        validated_its_boundaries[acc] = assign(coordinate_18, coordinate_58)
                        #print "2 18s validati, 5.8 appartiene al 2 dom -", acc
                    #print "test:", acc
                    #print coordinate_18_1[4], coordinate_18_0[4], rslts58[acc][4]
                    if (coordinate_58[4] < coordinate_18_0[4]) and \
                            (coordinate_18_0[4] < coordinate_18_1[4]):
                        #print "test:", acc
                        # il 5.8 appartiene al primo 18s
                        coordinate_18 = coordinate_18_0
                        validated_its_boundaries[acc] = assign(coordinate_18, coordinate_58)
                        #print "2 18s validati, 5.8 appartiene al 1 dom -", acc
    #print "stop validation_concordance()"
    return validated_its_boundaries


def last_check(validated_its_boundaries, sample_name, path_fasta):
    """
    For each validated ITS1, this function verifies that the % of uncalled bases into
    the ITS1 sequence, is lower than 3%. The same threshold is applied for the possible
    introns that are in the 18S genes. Also, the sequences with a distance between the
    beginning of 18S alignment and the beginning of sequence upper than 25 are discarded.
    :param validated_its_boundaries:a dictionary with the accession numbers as keys
    and the tuplas with ITS1 boundaries as value.
    :param sample_name:
    :param path_fasta:
    :return:validated_its_boundaries:a dictionary with the accession numbers as keys
    and the tuplas with ITS1 boundaries as value, cleaned of the ITS1 that don't pass
    the filters.
    """
    if len(validated_its_boundaries.keys()) != 0:
        #print "start parsing fasta %s %s" % (sample_name, str(datetime.datetime.now().time()))
        read_fasta = {}
        #parsing del file fasta
        with gzip.open(path_fasta, 'rt') as a:
            line = a.readline()
            while line != "":
                if line.startswith(">"):
                    key = line[1:].strip()
                    line = a.readline()
                    val = line.strip()
                    read_fasta[key] = val
                    line = a.readline()
        #print "stop parsing fasta", datetime.datetime.now().time()
        for k, v in validated_its_boundaries.copy().items():
            seq = read_fasta[k]
            if "+" in v[-1]:
            #se l'orientamento del cluster e' positivo
            #per coverage minori di 0.85 scarto le sequenze con 18s-from maggiore di 25
                if ((int(v[3]) + 1) > 25) and (float(v[1]) <= 0.85):
                    excluded_ITS.write("[%s] %s removed: 18S-from distance > 20.\n(%s)\n" % (sample_name, k, str(v)))
                    try:
                        validated_its_boundaries.pop(k)
                    except KeyError:
                        pass
                    pass
                # conta la percent di N in INTRONE. se inferiore a 3 percent ok
                if "join" in v[-1]:
                    inti = v[-1].split("...")
                    int_i = inti[1].split(",")
                    intron = seq[int(int_i[0]):int(int_i[1])]
                    if len(intron) != 0:
                        #print k, v
                        inf = intron.count("n")/float(len(intron))
                        if inf > 0.03:
                            excluded_ITS.write("[%s] %s removed: N percentage in the intron %f\n(%s)\n" % (sample_name, k, inf, str(v)))
                            try:
                                validated_its_boundaries.pop(k)
                            except KeyError:
                                pass
                            pass
                #conta la percent di N in ITS1. se inferiore a 3 percent ok
                its = seq[int(v[4]):int(v[8])]
                its_nf = its.count("n")/float(len(its))
                if its_nf > 0.03:
                    excluded_ITS.write("[%s] %s removed: N percentage in the ITS1 %f\n(%s)\n" % (sample_name, k, its_nf, str(v)))
                    try:
                        validated_its_boundaries.pop(k)
                    except KeyError:
                        pass
                    pass
            if "-" in v[-1]:
                if ((int(v[-2]) - int(v[3]) + 1) > 25) and (float(v[1]) <= 0.85):
                    excluded_ITS.write("[%s] %s removed: 18S-from distance > 20.\n(%s)\n" % (sample_name, k, str(v)))
                    try:
                        validated_its_boundaries.pop(k)
                    except KeyError:
                        pass
                    pass
                if "join" in v[-1]:
                    inti = v[-1].split("...")
                    int_i = inti[1].split(",")
                    intron = seq[int(int_i[1]):int(int_i[0])]
                    if len(intron) != 0:
                        inf = intron.count("n")/float(len(intron))
                        if inf > 0.03:
                            excluded_ITS.write("[%s] %s removed: N percentage in the intron %f\n(%s)\n" % (sample_name, k, inf, str(v)))
                            try:
                                validated_its_boundaries.pop(k)
                            except KeyError:
                                pass
                            pass
                its = seq[int(v[8]):int(v[4])]
                its_nf = its.count("n")/float(len(its))
                if its_nf > 0.03:
                    excluded_ITS.write("[%s] %s removed: N percentage in the ITS1 %f\n(%s)\n" % (sample_name, k, its_nf, str(v)))
                    try:
                        validated_its_boundaries.pop(k)
                    except KeyError:
                        pass
                    pass
    return validated_its_boundaries


def extraction_ITS1(path18, path58, sample_name, path_fasta):
    """In order to analyze the only query sequences that contemporary align with
    18s and 5.8s hmm profiles, the $nhmmer output files are parsed and the common
    accessions are included into a list object. The 18s and 5.8s dictionary objects
    are created. The dictionary keys are the accession numbers, the values are the
    integer numbers that represent the line positional indexes of accession numbers
    into the readlines() list of strings,
    The uncommon accession are popped by dictionaries.
    Then data_extra_valid_18() and data_extra_valid_58() are recalled, followed by
    validation_concordance() and last_check().
    This function return the ITS1 boundaries dictionary.
    :param path18: the path string of align_18S_data.txt file.
    :param path58: the path string of align_5_8S_data.txt file.
    :param sample_name: the name of ENA section of the sequences.
    :param path_fasta: the path string of related fasta file.
    """
    #print "start accession_extraction()"
    accession_18, accession_58, read_fasta, common_acc = {}, {}, {}, []
    #indicizza, sia per 18S che per 5.8S, i file output di nhmmer per l'allineamento
    # estraendo gli accession delle sequenze object allineate ed associando
    # all'accession il numero o i numeri di rigo su cui sono riportati i dati
    # di allineamento.
    with gzip.open(path18, "rt") as read_18:
        lines_18 = read_18.readlines()
    cursor = 0
    for line in lines_18:
        if line[0:2] == ">>":
            s = line.split(" ")
            #print cursor, s[1]
            accession_18.setdefault(s[1], [])
            accession_18[s[1]].append(cursor)
            cursor += 1
        else:
            cursor += 1
    with gzip.open(path58, "rt") as read_58:
        lines_58 = read_58.readlines()
    cursor = 0
    for line in lines_58:
        if line[0:2] == ">>":
            s = line.split(" ")
            accession_58.setdefault(s[1], [])
            accession_58[s[1]].append(cursor)
            cursor += 1
        else:
            cursor += 1
    #individua gli accession number comuni tra i due allineamenti e cancella
    #quelli non comuni
    common_acc = set(accession_18.keys()).intersection(set(accession_58.keys()))
    return_time("""- - - -
    [%s]>> ACCESSION EXTRACTION
    [%s]accession_18: %s accession_58: %s common_accession: %s\n""" % (
        sample_name, sample_name,
        len(accession_18.values()), len(accession_58.keys()), len(common_acc)))
    for k, v in accession_18.copy().items():
        if k not in common_acc:
            accession_18.pop(k)
    for k1, v1 in accession_58.copy().items():
        if k1 not in common_acc:
            accession_58.pop(k1)
    #richiamo le due funzioni data_extr_valid_18 data_extr_valid_58,
    # validation_concordance e last_check
    results18 = data_extr_valid_18(accession_18, lines_18)
    results58 = data_extr_valid_58(accession_58, lines_58)
    validated_its_boundaries = validation_concordance(results18, results58, sample_name)
    #print "stop accession_extraction()"
    return last_check(validated_its_boundaries, sample_name, path_fasta)


def import_tax():
    lineage = dict()
    with open(os.path.join(basepath,
                           'taxonomy_%s' % release,
                           'taxidlineage.dmp'), 'rt') as taxy:
        for line in taxy:
            s = list(map(str.strip, line.split('|')))
            s1 = s[1].split(' ')
            lineage[s[0]] = s1

    with open(os.path.join(basepath,
                           'taxonomy_%s' % release,
                           "delnodes.dmp"), "rt") as delnodes:
        delnodes_list = dict()
        line = delnodes.readline()
        while line:
            delnodes_list.setdefault(line.split(" ")[0], True)
            line = delnodes.readline()

    merged = dict()
    with open(os.path.join(basepath,
                           'taxonomy_%s' % release,
                           'merged.dmp'), 'rt') as taxy:
        for line in taxy:
            s = list(map(str.strip, line.split('|')))
            merged[s[0]] = s[1]

    return lineage, delnodes_list, merged


def print_result(boundaries_dict, sample_name, tmp_path):
    """
    This funtion writes the results in the temporary output file.
    :param boundaries_dict: the dictionary of validated ITS1 boundaries.
    :param sample_name: The folder name of the analyzed file-sequences type.
    :param tmp_path: the folder in which are put the temporary output files
    :return:
    """
    #print "start print_results()"
    #per ogni cartella vengono printati i risultati in un file temporaneo corrisponde
    lineage, deleted, merged = import_tax()
    if len(boundaries_dict.keys()) != 0:
        tmp_name = "%s.tmp" % sample_name
        path_tmp_file = os.path.join(tmp_path, tmp_name)
        # print(path_tmp_file)
        with open(path_tmp_file, "wt") as results:
            for k, v in boundaries_dict.items():
                #print v
                if k.split('|')[1] in list(lineage.keys()):
                    taxtry = k.split('|')[1]
                    # print(taxtry)
                else:
                    if k.split('|')[1] in list(deleted.keys()):
                        taxtry = '77133'
                    elif k.split('|')[1] in list(merged.keys()):
                        taxtry = merged[k.split('|')[1]]
                        # print(taxtry)
                    else:
                        taxtry = '77133'
                if '2157' not in lineage[taxtry] and '2' not in lineage[taxtry]:
                    if len(v) == 12:
                        results.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
                            str(k), str(int(v[-2])), str(v[0]), str(v[1]), str(v[2]), str(v[3]), str(v[4]),
                            str(v[5]), str(v[6]), str(v[7]), str(v[8]), str(v[9]), v[-1], sample_name))
                    else:
                        return_time("Len tupla %s. Nessun risultato viene printato per l'accession %s\n"%
                              (str(len(v)), k))
                else:
                    excluded_ITS.write('%s is not an Eukaryotes\n' % k)
                    excluded_ITS.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n\n" % (
                            str(k), str(int(v[-2])), str(v[0]), str(v[1]), str(v[2]), str(v[3]), str(v[4]),
                            str(v[5]), str(v[6]), str(v[7]), str(v[8]), str(v[9]), v[-1], sample_name))
    #print "stop print_results()"
    return


def worker(arg):#folder, path_folder_of_folders, fasta_folder, tmp_path):
    """
    This function analyzes the single folder in which are stored the output files of
    HMM 18S e 5.8S for a single fasta file.
    :param folder: are stored the output files of
    HMM 18S e 5.8S for a single fasta file
    :param path_folder_of_folders: path of up folder of the analyzed folder
    :param tmp_path: path of temporary files' folder.
    :return:
    """
    folder = arg[0]
    path_folder_of_folders = arg[1]
    fasta_folder = arg[2]
    tmp_path = arg[3]
    if os.path.isdir(os.path.join(path_folder_of_folders, folder)):
        path_18, path_58 = "", ""
        return_time("The '%s' dir elaboration started\n" % folder)
        name_fasta = "%s.fasta.gz" % folder
        path_fasta = os.path.join(fasta_folder, name_fasta)#path_folder_of_folders
        sample_split = folder.split("_")
        name = "_".join(sample_split[1:4])
        chpath = os.path.join(path_folder_of_folders, folder)
        l = os.listdir(chpath)
        if "align_18S_data.txt.gz" in l:
            path_18 = os.path.join(chpath, "align_18S_data.txt.gz")
            # print path_18
        else:
            return_time("The align_18S_data.txt file there isn't in %s folder\n" % folder)
        if "align_5_8S_data.txt.gz" in l:
            path_58 = os.path.join(chpath, "align_5_8S_data.txt.gz")
        else:
            return_time("The align_5_8S_data.txt file there isn't in %s folder\n" % folder)
        if os.path.isfile(path_18) and os.path.isfile(path_58):
            print_result(extraction_ITS1(path_18, path_58, name, path_fasta),
                         name, tmp_path)
        return_time("The '%s' finished\n" % folder)
    return 1


def main(path_folder_of_folders, path_fasta, outname, processes, release):
    """
    The function menages the other functions. The analysis for
    each folder are launched recursively in multi-process mode.
    :param path_folder_of_folders:the folder path string.
    :param outname: the name of the output file
    :return:
    """
    #print "start giostra()"
    if os.path.exists(path_folder_of_folders):
        out_name = outname
        tmp_path = os.path.join(path_folder_of_folders, "tmp")
        try:
            os.mkdir(tmp_path)
        except OSError:
            rmtree(tmp_path)
            os.mkdir(tmp_path)
        # print('+++', tmp_path)
        list_folders = [folder for folder in os.listdir(path_folder_of_folders)
                        if (folder.startswith("rel_") and folder.endswith("_r%s" % release))]
        ########debug mode single process
        # for folder in list_folders:
        # worker('rel_std_fun_02_r141', path_folder_of_folders, path_fasta, tmp_path)
        # def aux(folder):
        #     worker(folder, path_folder_of_folders, path_fasta, tmp_path)
        ####### multiprocess
        if len(list_folders) < int(processes):
            processes = len(list_folders)
            chunksize = 1
        else:
            processes = int(processes)
            chunksize= len(list_folders)//int(processes)
        with mp.Pool(processes=processes) as po:
            args = [(folder, path_folder_of_folders, path_fasta, tmp_path) for
                    folder in list_folders]
            future = po.map(func=worker, iterable=args,
                            chunksize=chunksize)
        # po.shutdown(wait=False)
        ######################################## OLD MULTIPROCESS PROCEDURE
        # list_folders = [folder for folder in os.listdir(path_folder_of_folders)
        #                 if (folder.startswith("rel_") and folder.endswith("_r%s" % release))]
        # shuffle(list_folders)
        # number_folders = len(list_folders)
        # list_tuplae = [(n, n+4) for n in xrange(0, number_folders, 4) if (n+4) < number_folders]
        # list_tuplae.append((list_tuplae[-1][-1], None))
        # for t in list_tuplae:
        #     for folder in list_folders[t[0]:t[1]]:
        #         Prs = mp.Process(target=worker, args=(folder, path_folder_of_folders, path_fasta, tmp_path,))
        #         Prs.start()
        #     print "The nr of running processes is %s" % len(mp.active_children())
        #     while len(mp.active_children()) > int(processes):
        #         pass
        # while len(mp.active_children()) > 0:
        #     pass
        # else:
        #######################################################################
        #if future.done():
        with open(out_name, "wt") as results:
            results.write(
                "ACCESSION\tVERSION\tTAXID\tLENGTH\tE-VALUE_18\tHMM-COVERAGE_18\tPP_18\tSTART_18\tEND_18\tE-VALUE_5.8\tHMM-COVERAGE_5.8\tPP_5.8\tSTART_5.8\tEND_5.8\tNOTE\tSAMPLE_NAME\n")
            for tmp in os.listdir(tmp_path):
                with open(os.path.join(tmp_path, tmp), "rt") as handle:
                    lines = handle.readlines()
                    for line in lines:
                        s = line.split('\t')
                        acc = s[0].split('.')[0]
                        vers, taxid = s[0].split('.')[1].split('|')
                        results.write('\t'.join([acc, vers, taxid, *s[1:]]))
        if future.count(1) == len(list_folders):
            pass
        else:
            print('\n##################### ERRORS OCCURRED #####################\n')
            print(format(list(future)))
    return


# def kill_child():
#     print('+++', os.geteuid())
#     # pool.shutdown(wait=False)


if __name__ == '__main__':
    path_align, path_fasta, nameout, processes, release = pars_arg()
    if path_align.endswith('/'):
        path = path_align[:-1]
    else:
        path = path_align
    basepath = path.split('/')
    basepath = "/".join(basepath[:-1])
    excluded_path = os.path.join(basepath, 'hmmer_excluded_ITS1.log')
    excluded_ITS = open(excluded_path, 'wt')
    main(path_align, path_fasta, nameout, processes, release)
    excluded_ITS.close()
