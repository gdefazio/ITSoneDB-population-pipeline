#!/usr/bin/env bash

help="
                     #####################################################
                     #                                                   #
                     #              ITSoneDB upgrade program             #
                     #               ETL_population_pipeline             #
                     #                                                   #
                     #####################################################


This script menages the procedure to populate etl directory after
ITSoneDB_upgrade_pipeline procedure completion.
USAGE:
$./ETL_population_pipeline.sh
        -s directory of scripts menaged by this .sh
        -n release number
        -r releases directory
        -p previous release number
"

if [[ $# -eq 0 ]]   			# se lo script è eseguito senza argomenti...
	then
  		echo -e "$help"		# ... stampo l'help...
  		exit 1         		# ... ed esco.
fi
while getopts ":h:s:n:r:p:" opzione
do
	case "${opzione}" in
	s) scripts=${OPTARG};;
	#x) aux=${OPTARG};;
	n) release=${OPTARG};;
	r) releases=${OPTARG};;
	p) prevrel=${OPTARG};;
	h) echo -e "\n$help";exit 1;;
	#*) echo -e "\nÈ stata inserita un'opzione non valida.\n\n $help";exit 1;;
	esac
done

function exit_verify()
{
    if [ ${?} -eq 0 ]
    then
        true
    else
        exit 1
    fi
}

echo "                     #####################################################"
echo "                     #                                                   #"
echo "                     #              ITSoneDB upgrade program             #"
echo "                     #               ETL_population_pipeline             #"
echo "                     #                                                   #"
echo "                     #####################################################"
echo ""
echo ""
echo ""
scripts=`realpath -s {scripts}`
releases=`realpath -s {releases}`
echo $scripts
#x) aux=${OPTARG};;
echo $release
echo $releases
echo $prevrel
echo ""
echo ""
echo "############################################################################################"
echo "${scripts}/0_splitter.py
                -t ${releases}/${release}/tsv_${release}/
                -f ${releases}/${release}/fasta_${release}/
                -i ${releases}/${release}/ITS1_loc_final_table.csv
                -o ${releases}/${release}/etl"
${scripts}/0_splitter.py\
                -t ${releases}/${release}/tsv_${release}/\
                -f ${releases}/${release}/fasta_${release}/\
                -i ${releases}/${release}/ITS1_loc_final_table.csv\
                -o ${releases}/${release}/etl
exit_verify

echo ""
echo ""
echo "############################################################################################"
echo "${scripts}/1_correct_final_table_taxid.py
                -i ${releases}/${release}/ITS1_loc_final_table.csv
                -t ${releases}/${release}/taxonomy_${release}
                -m ${releases}/${release}/excluded_ITS1_ENA.csv
                -p ${releases}/${prevrel}/taxonomy_${prevrel}
                -o ${releases}/${release}/etl"
${scripts}/1_correct_final_table_taxid.py\
                -i ${releases}/${release}/ITS1_loc_final_table.csv\
                -t ${releases}/${release}/taxonomy_${release}\
                -m ${releases}/${release}/excluded_ITS1_ENA.csv\
                -p ${releases}/${prevrel}/taxonomy_${prevrel}\
                -o ${releases}/${release}/etl
exit_verify
echo ""
echo ""
echo "############################################################################################"
echo "${scripts}/2_create_mapping_ENA_ITSoneDB_acc.py
                -i ${releases}/${release}/etl/ITS1_loc_final_table_corrected.tsv
                -o ${releases}/${release}/etl
                -r ${release}
                -f ${releases}/${release}/etl/acc2fasta_info.csv.gz
                -m ${releases}/${prevrel}/etl/mapping_data_ENA_release_${prevrel}.tsv
                -n ${releases}/${release}/taxonomy_${release}"
${scripts}/2_create_mapping_ENA_ITSoneDB_acc.py\
                -i ${releases}/${release}/etl/ITS1_loc_final_table_corrected.tsv\
                -o ${releases}/${release}/etl\
                -r ${release}\
                -f ${releases}/${release}/etl/acc2fasta_info.csv.gz\
                -m ${releases}/${prevrel}/etl/mapping_data_ENA_release_${prevrel}.tsv\
                -n ${releases}/${release}/taxonomy_${release}
exit_verify
echo ""
echo ""
echo "############################################################################################"
echo "${scripts}/3_prepare_gbentry.py
                -i ${releases}/${release}/etl/ITS1_loc_final_table_corrected.tsv
                -f ${releases}/${release}/etl/acc2fasta_info.csv.gz
                -t ${releases}/${release}/etl/acc2tsv_info.csv.gz
                -o ${releases}/${release}/etl"
${scripts}/3_prepare_gbentry.py\
                -i ${releases}/${release}/etl/ITS1_loc_final_table_corrected.tsv\
                -f ${releases}/${release}/etl/acc2fasta_info.csv.gz\
                -t ${releases}/${release}/etl/acc2tsv_info.csv.gz\
                -o ${releases}/${release}/etl
exit_verify
echo ""
echo ""
echo "############################################################################################"
echo "${scripts}/4_Alignment_splitter.py
                -m ${releases}/${release}/ITS1_loc_NHMMER.csv
                -a ${releases}/${release}/align_${release}
                -r ${release}
                -o ${releases}/${release}/etl"
${scripts}/4_Alignment_splitter.py\
                -m ${releases}/${release}/ITS1_loc_NHMMER.csv\
                -a ${releases}/${release}/align_${release}\
                -r ${release}\
                -o ${releases}/${release}/etl
exit_verify
echo ""
echo ""
echo "############################################################################################"
echo "${scripts}/5_prepare_flanking_region.py
                -m ${releases}/${release}/etl/mapping_data_ENA_release_${release}.tsv
                -b ${releases}/${release}/ITS1_loc_NHMMER.csv
                -f ${releases}/${release}/etl/hmm_match_data/alignment_18S
                -r ${releases}/${release}/etl/hmm_match_data/alignment_5_8S
                -a ${releases}/${release}/etl/acc2fasta_info.csv.gz
                -o ${releases}/${release}/etl"
${scripts}/5_prepare_flanking_region.py\
                -m ${releases}/${release}/etl/mapping_data_ENA_release_${release}.tsv\
                -b ${releases}/${release}/ITS1_loc_NHMMER.csv\
                -f ${releases}/${release}/etl/hmm_match_data/alignment_18S\
                -r ${releases}/${release}/etl/hmm_match_data/alignment_5_8S\
                -a ${releases}/${release}/etl/acc2fasta_info.csv.gz\
                -o ${releases}/${release}/etl
exit_verify
echo ""
echo ""
echo "############################################################################################"
echo "${scripts}/6_find_representative.py
                -f ${releases}/${release}/etl/FASTA_STATICI/ITSoneDB_total_fasta_rel${release}.fa.gz
                -n ${releases}/${release}/taxonomy_${release}
                -o ${releases}/${release}/etl
                -p 20"
${scripts}/6_find_representative.py\
                -f ${releases}/${release}/etl/FASTA_STATICI/ITSoneDB_total_fasta_rel${release}.fa.gz\
                -n ${releases}/${release}/taxonomy_${release}\
                -o ${releases}/${release}/etl\
                -p 20
exit_verify
echo ""
echo ""
echo "############################################################################################"
echo "${scripts}/7_rep_fasta_generation.py
                -r ${release}
                -c ${releases}/${release}/etl/Representative_computation/rep_sequences.csv
                -i ${releases}/${release}/etl/FASTA_STATICI/ITSoneDB_total_fasta_rel${release}.fa.gz
                -f ${releases}/${release}/etl/FASTA_STATICI/ITSoneDB_total_fasta_rel${release}_with_flanking.fa
                -o ${releases}/${release}/etl/Representative_computation"
${scripts}/7_rep_fasta_generation.py\
                -r ${release}\
                -c ${releases}/${release}/etl/Representative_computation/rep_sequences.csv\
                -i ${releases}/${release}/etl/FASTA_STATICI/ITSoneDB_total_fasta_rel${release}.fa.gz\
                -f ${releases}/${release}/etl/FASTA_STATICI/ITSoneDB_total_fasta_rel${release}_with_flanking.fa.gz\
                -o ${releases}/${release}/etl/Representative_computation
exit_verify
echo ""
echo ""
echo "############################################################################################"
echo "${scripts}/8_taxonomy_cleaner.py
                -n ${releases}/${release}/taxonomy_${release}
                -l ${releases}/${release}/etl/ITS1_loc_final_table_corrected.tsv
                -o ${releases}/${release}/etl"
${scripts}/8_taxonomy_cleaner.py\
                -n ${releases}/${release}/taxonomy_${release}\
                -l ${releases}/${release}/etl/ITS1_loc_final_table_corrected.tsv\
                -o ${releases}/${release}/etl

echo ""
echo ""
echo "############################################################################################"
echo "${scripts}/9_prepare_fasta_input.py
                -m ${releases}/${release}/etl/mapping_data_ENA_release_${release}.tsv
                -f ${releases}/${release}/etl/FASTA_STATICI/ITSoneDB_total_fasta_rel${release}.fa.gz
                -n ${releases}/${release}/etl/cleaned_ncbi_taxonomy
                -o ${releases}/${release}/etl/FASTA_STATICI"
${scripts}/9_prepare_fasta_input.py\
                -m ${releases}/${release}/etl/mapping_data_ENA_release_${release}.tsv\
                -f ${releases}/${release}/etl/FASTA_STATICI/ITSoneDB_total_fasta_rel${release}.fa.gz\
                -n ${releases}/${release}/etl/cleaned_ncbi_taxonomy\
                -o ${releases}/${release}/etl/FASTA_STATICI
exit_verify

echo ""
echo ""
echo "############################################################################################"
echo "${scripts}/10_prepare_its1feauture_file.py
                -i ${releases}/${release}/etl/mapping_data_ENA_release_${release}.tsv
                -l ${releases}/${release}/etl/ITS1_loc_final_table_corrected.tsv
                -o ${releases}/${release}/etl"
${scripts}/10_prepare_its1feauture_file.py\
                -i ${releases}/${release}/etl/mapping_data_ENA_release_${release}.tsv\
                -l ${releases}/${release}/etl/ITS1_loc_final_table_corrected.tsv\
                -o ${releases}/${release}/etl
exit_verify
