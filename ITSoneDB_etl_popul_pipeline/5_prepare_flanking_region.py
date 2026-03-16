#!/usr/bin/env python3
# coding=utf-8
from librs.return_time import return_time
import os
import argparse
import argcomplete
import sys
import gzip


def input_options():
    parser = argparse.ArgumentParser(description="split fasta and tsv files", prefix_chars="-")
    parser.add_argument("-m", "--mapping_file", type=str,
                        help="mapping file created by the 2_create_mapping_ENA_ITSoneDB_acc.py",
                        action="store", required=True)
    parser.add_argument("-b", "--hmm_results_file", type=str,
                        help="tsv containing ITS1 location inferred by HMM",
                        action="store", required=True)
    parser.add_argument("-f", "--r18S_align_res", type=str,
                        help="folder containing the 18S align results",
                        action="store", required=True)
    parser.add_argument("-r", "--r5_8S_align_res", type=str,
                        help="folder containing the 5.8S align results",
                        action="store", required=True)
    parser.add_argument("-a", "--fasta_tsv", type=str,
                        help="tsv file containing the fasta file location",
                        action="store", required=True)
    parser.add_argument("-o", "--output_folder", type=str, help="output folder", action="store", required=False,
                        default="etl")
    argcomplete.autocomplete(parser)
    return parser.parse_args()


def hmm_data(hmm_tsv):
    """
    This function associate ENA accession number to hmm alignment data
    :param hmm_tsv: string
    :return data_align: dict
    """
    data_aling = {}
    with open(hmm_tsv, 'rt') as a:
        a.readline()
        for linea in a:
            s = list(map(str.strip, linea.split("\t")))
            data_aling[s[0]] = s
    return data_aling


def def_len(tsv):
    """

    :param tsv: str
    :return acc_len: dict
    """
    from Bio import SeqIO
    acc_len = {}
    with gzip.open(tsv, 'rt') as a:
        a.readline()
        for linea in a:
            acc, path = list(map(str.strip, linea.split("\t")))
            with gzip.open(os.path.join(out.split('/etl')[0],path), 'rt') as b:
                record = SeqIO.read(b, "fasta")
                acc_len[acc] = len(record)
    return acc_len


if __name__ == "__main__":
    param = input_options()
    mapping, hmm_res, folder18S_align, folder58S_align, out, fasta_tsv = param.mapping_file, param.hmm_results_file, param.r18S_align_res, param.r5_8S_align_res, param.output_folder, param.fasta_tsv
    if folder18S_align.startswith('/home') and folder58S_align.startswith('/home'):
        folder18S_align = 'etl%s' % folder18S_align.split('/etl')[-1]
        folder58S_align = 'etl%s' % folder58S_align.split('/etl')[-1]
    hmm_info = hmm_data(hmm_res)
    return_time('HMM table parsing DONE')
    acc2len = def_len(fasta_tsv)
    return_time('acc2fasta_info parsing DONE')
    tmp = open(os.path.join(out, "flankingregion.tsv"), "wt")
    tmp.write(
        "GBentry_Accession\tAccessionVersion\thasGBannotation\tGBITS1Accession\thasHMM\tHMMITS1Accession\tGBfeatureKey\tDescription\tQualifier\tValue\tStart\tEnd\tLefComplete\tRightcomplete\tStart2\tEnd2\tLeftComplete\tRightComplete\teValue\tscore\tposteriorProbability\talignmentInfo\n")
    count = 1
    hmm18_count = 1
    hmm58_count = 1
    with open(mapping, 'rt') as l:
        l.readline()
        for line in l:
            field = list(map(str.strip, line.split("\t")))
            # return_time(' '.join([field[0], str(count)]))
            ssu_stringa = [field[0], field[1]]  # GBentry_Accession,AccessionVersion
            stringa_58s = [field[0], field[1]]  # GBentry_Accession,AccessionVersion
            if field[2] != "Null":
                ssu_stringa.append("1")  # hasGBannotation
                ssu_stringa.append(field[2])  # GBITS1Accession
                stringa_58s.append("1")  # hasGBannotation
                stringa_58s.append(field[2])  # GBITS1Accession
            else:
                ssu_stringa.append("0")  # hasGBannotation
                ssu_stringa.append("0")  # GBITS1Accession
                stringa_58s.append("0")  # hasGBannotation
                stringa_58s.append("0")  # GBITS1Accession
            if field[6] != "Null":
                ssu_stringa.append("1")  # hasHMM
                ssu_stringa.append(field[6])  # HMMITS1Accession
                stringa_58s.append("1")  # hasHMM
                stringa_58s.append(field[6])  # HMMITS1Accession
            else:
                ssu_stringa.append("0")  # hasHMM
                ssu_stringa.append("0")  # HMMITS1Accession
                stringa_58s.append("0")  # hasHMM
                stringa_58s.append("0")  # HMMITS1Accession
            if ssu_stringa[2] != "0":
                start_ITS1 = int(field[3].replace("<", "").replace(">", ""))
                if start_ITS1 == 1:
                    # non c'e' la regione fiancheggiante di sinistra
                    start_left_flank = 'NULL'
                    end_left_flank = 'NULL'
                else:
                    # c'e' la regione fiancheggiante di sinistra
                    end_left_flank = start_ITS1 - 1
                    # se la fine della regione fiancheggiante di sinistra
                    # coincide con l'inizio della sequenza
                    if end_left_flank == 1:
                        # la regione fiancheggiante è di un singolo nt quindi start ed end coincidono
                        start_left_flank = 1
                    else:
                        #calcolo l'inizio della flank sinistra
                        start_left_flank = end_left_flank - 150
                        if start_left_flank < 1:
                            start_left_flank = 1
                end_ITS1 = int(field[4].replace("<", "").replace(">", ""))
                # se la fine della its1 coincide con la fine della sequenza
                # non c'è una flanking region
                if end_ITS1 == acc2len[field[0]]:
                    start_right_flank = 'NULL'
                    end_right_flank = 'NULL'
                elif end_ITS1 < acc2len[field[0]]:
                    start_right_flank = end_ITS1 + 1
                    if start_right_flank == acc2len[field[0]]:
                        # la regione fiancheggiante è di un singolo nt quindi start ed end coincidono
                        end_right_flank = acc2len[field[0]]
                    else:
                        end_right_flank = start_right_flank + 150
                        if end_right_flank > acc2len[field[0]]:
                            end_right_flank = acc2len[field[0]]
                else:
                    sys.exit('%s\t%s' % (end_ITS1, acc2len[field[0]]))
                ssu_stringa.append("rRNA")  # GBfeatureKey
                ssu_stringa.append("18S")  # Description
                ssu_stringa.append("product")  # Qualifier
                ssu_stringa.append("18S ribosomal RNA")  # Value
                ssu_stringa.append(str(start_left_flank))  # Start
                ssu_stringa.append(str(end_left_flank)) # End #int(field[3].replace("<", "").replace(">", "")) - 1))  # End
                ssu_stringa.append("0")  # LefComplete
                ssu_stringa.append("0")  # Rightcomplete
                stringa_58s.append("rRNA")  # GBfeatureKey
                stringa_58s.append("5.8S")  # Description
                stringa_58s.append("product")  # Qualifier
                stringa_58s.append("5.8S ribosomal RNA")  # Value
                stringa_58s.append(str(start_right_flank)) #int(field[4].replace("<", "").replace(">", "")) + 1))  # Start
                stringa_58s.append(str(end_right_flank))  # End
                stringa_58s.append("0")  # LefComplete
                stringa_58s.append("0")  # Rightcomplete
            else:
                ssu_stringa.append("NULL")  # GBfeatureKey
                ssu_stringa.append("18S")  # Description
                ssu_stringa.append("NULL")  # Qualifier
                ssu_stringa.append("NULL")  # Value
                ssu_stringa.append("NULL")  # start
                ssu_stringa.append("NULL")  # end
                ssu_stringa.append("NULL")  # leftcomplete
                ssu_stringa.append("NULL")  # rightcomplete
                stringa_58s.append("NULL")  # GBfeatureKey
                stringa_58s.append("5.8S")  # Description
                stringa_58s.append("NULL")  # Qualifier
                stringa_58s.append("NULL")  # Value
                stringa_58s.append("NULL")  # start
                stringa_58s.append("NULL")  # end
                stringa_58s.append("NULL")  # leftcomplete
                stringa_58s.append("NULL")  # rightcomplete
            if field[6] != "Null":
                ssu_stringa.append(hmm_info[field[0]][7])  # Start --> questo è da verificare
                ssu_stringa.append(hmm_info[field[0]][8])  # End --> questo è da verificare
                ssu_stringa.append("1")  # LeftComplete
                ssu_stringa.append("1")  # RightComplete
                ssu_stringa.append(hmm_info[field[0]][4])  # eValue --> questo è da verificare
                ssu_stringa.append("NULL")  # score
                ssu_stringa.append(hmm_info[field[0]][6])  # posteriorProbability --> questo è da verificare
                ssu_stringa.append(os.path.join(folder18S_align, "%s.align" % field[0]))  # alignmentInfo
                if not os.path.exists(os.path.join(folder18S_align, "%s.align" % field[0])):
                    #return_time(os.path.join(folder18S_align, "%s.align" % field[0]))
                    hmm18_count += 1
                stringa_58s.append(hmm_info[field[0]][12])  # Start --> questo è da verificare
                stringa_58s.append(hmm_info[field[0]][13])  # End  --> questo è da verificare
                stringa_58s.append("1")  # LeftComplete
                stringa_58s.append("1")  # RightComplete
                stringa_58s.append(hmm_info[field[0]][9])  # eValue --> questo è da verificare
                stringa_58s.append("NULL")  # score
                stringa_58s.append(hmm_info[field[0]][11])  # posteriorProbability --> questo è da verificare
                stringa_58s.append(os.path.join(folder58S_align, "%s.align" % field[0]))  # alignmentInfo
                if not os.path.exists(os.path.join(folder58S_align, "%s.align" % field[0])):
                    # return_time(os.path.join(folder58S_align, "%s.align" % field[0]))
                    hmm58_count += 1
            else:
                ssu_stringa.append("NULL")  # 16
                ssu_stringa.append("NULL")  # 17
                ssu_stringa.append("NULL")  # 18
                ssu_stringa.append("NULL")  # 19
                ssu_stringa.append("NULL")  # 20
                ssu_stringa.append("NULL")  # 21
                ssu_stringa.append("NULL")  # 22
                ssu_stringa.append("NULL")  # 23
                stringa_58s.append("NULL")  # 16
                stringa_58s.append("NULL")  # 17
                stringa_58s.append("NULL")  # 18
                stringa_58s.append("NULL")  # 19
                stringa_58s.append("NULL")  # 20
                stringa_58s.append("NULL")  # 21
                stringa_58s.append("NULL")  # 22
                stringa_58s.append("NULL")  # 23
            tmp.write("%s\n" % "\t".join(ssu_stringa))
            tmp.write("%s\n" % "\t".join(stringa_58s))
            count += 1
    tmp.close()
    return_time("%s Accession elaborated for %s" % (count, os.path.join(out, "flankingregion.tsv")))
    return_time("%s 18S mapping reported for %s" % (hmm18_count, os.path.join(out, "flankingregion.tsv")))
    return_time("%s 58S mapping reported for %s" % (hmm58_count, os.path.join(out, "flankingregion.tsv")))

