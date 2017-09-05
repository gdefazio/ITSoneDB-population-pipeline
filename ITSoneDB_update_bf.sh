#!/bin/bash
# This script computes the 18S and 5.8S HMM profile mapping against the candidate sequences

# Add here the path to script folder
script=/path/to/script/folder/

main_folder=`pwd`
for compressed_fasta in `cat sample_list`
    do
    echo "STEP1 - FASTA extraction"
    echo "Starting at"
    echo `date`
    echo
    gzip -d ${compressed_fasta}
    fasta=$(basename ${compressed_fasta} .gz)

    if [[ -e ${fasta} && -s ${fasta} ]]
        then
        echo "STEP1 - Ending at"
        echo `date`
    else
        echo "Errors during the FASTA estraction"
        exit 1
    fi

    echo "STEP2 - HMMER execution"
    echo "Starting at"
    echo `date`
    echo

    folder=$(basename ${fasta} .fasta)
    mkdir ${folder}
    cd ${folder}
    mkdir hmm_match_data
    hmmsearch -o hmm_match_data/align_5_8S_data.txt --cpu 8 ''${main_folder}''/hmm_motifs/solo_fungi_5_8S.hmm ../''${fasta}''
    hmmsearch -o hmm_match_data/align_18S_data.txt  --cpu 8 ''${main_folder}''/hmm_motifs/total_rfam_18S.hmm ../''${fasta}''

    if [[ -e hmm_match_data/align_5_8S_data.txt && -s hmm_match_data/align_5_8S_data.txt ]] && [[ -e hmm_match_data/align_18S_data.txt && -s hmm_match_data/align_18S_data.txt ]]
        then
        echo "STEP2 - Ending at"
        echo `date`
        echo "STEP3 - Extraction of alignment file"
        echo "Starting at"
        echo `date`
        echo
        python ''${main_folder}''/${script}/hmmer_txt_parser.py

        echo "STEP3 - Ending at"
        echo `date`
        echo

        echo "STEP4 - Inference of ITS1 location"
        echo "Starting at"
        echo `date`
        echo
        python ''${main_folder}''/${script}/estrazione_localizzazione_from_hmm.py -f ../''${fasta}''
        echo "STEP4 - Ending at"
        echo `date`
        echo
    else
        echo "Errors during the HMMER execution"
        echo ${fasta}
    fi

    cd ${main_folder}
    gzip ${fasta}
    done