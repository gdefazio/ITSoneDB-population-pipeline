#!/usr/bin/env bash

help="

This script menages the procedure to extract ITS1 location from ENA notations
and from NHMMER 18S and 5.8 mapping procedure.
USAGE:
$./ITSoneDB_upgrade_pipeline.sh
        -s full path directory containing scripts executed here
        -x full path aux directory
        -r full path releases directory
        -c cpus number
"

if [[ $# -eq 0 ]]   			# se lo script è eseguito senza argomenti...
	then
  		echo -e "$help"		# ... stampo l'help...
  		exit 1         		# ... ed esco.
fi
while getopts ":h:s:x:r:c:" opzione
do
	case "${opzione}" in
	s) scripts=${OPTARG};;
	x) aux=${OPTARG};;
	r) releases=${OPTARG};;
	c) cpus=${OPTARG};;
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
echo "                     #      Welcome to the ITSoneDB upgrade program      #"
echo "                     #                                                   #"
echo "                     #####################################################"
echo ""
echo ""
echo ""
wget -O ./release_doc  ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/doc/Release_* -q
exit_verify
release=$(${scripts}/0_verify_new_ena_release.py -d ./release_doc -r ${releases})
exit_verify
# echo ${release}
rm ./release_doc
echo "This program try to find the new release of ENA"

if [[ ${release} == "false" ]]
then
echo "You are not lucky today! There is not a new release of ENA.
    You can try later for the upgrade!
    See you later.
    Bye."
exit 0
else
releases=$(realpath -s ${releases})
scripts=$(realpath -s ${scripts})
aux=$(realpath -s ${aux})

echo ""
echo "New release found: ${release}"
echo "START THE UPGRADE OF ITSoneDB BASED ON ENA RELEASE ${release}"

echo ""
echo ""
echo "############################################################################################"
echo "${scripts}/1_download_ena_sequences.sh
                               -r ${release}
                               -d ${releases}/${release}/dat_${release}"
${scripts}/1_download_ena_sequences.sh\
                            -r ${release}\
                            -d ${releases}/${release}/dat_${release}
exit_verify
echo ""
echo ""
echo "############################################################################################"
echo "${scripts}/2_ENA_parser.py
                            -r ${release}
                            -i ${releases}/${release}/dat_${release}
                            -p ${cpus}"
${scripts}/2_ENA_parser.py\
                    -r ${release}\
                    -i ${releases}/${release}/dat_${release}\
                    -p ${cpus}
exit_verify
echo ""
echo ""
echo "############################################################################################"
echo "${scripts}/3_indexing.py
                            -r ${release}
                            -i ${releases}/${release}/fasta_${release}"
${scripts}/3_indexing.py\
                -r ${release}\
                -i ${releases}/${release}/fasta_${release}
exit_verify
echo ""
echo ""
echo "############################################################################################"
echo "${scripts}/4_exec_nhmmer.sh
                             -r ${release}
                             -t ${cpus}
                             -i ${releases}/${release}/fasta_${release}
                             -g ${aux}"
${scripts}/4_exec_nhmmer.sh -r ${release}\
                             -t ${cpus}\
                             -i ${releases}/${release}/fasta_${release}\
                             -g ${aux}
exit_verify
echo ""
echo ""
echo "############################################################################################"
echo "${scripts}/5_extract_ITS1_loc_from_nhmmer.py
                                             -f ${releases}/${release}/align_${release}
                                             -s ${releases}/${release}/fasta_${release}
                                             -p ${cpus}
                                             -r ${release}
                                             -o ${releases}/${release}/ITS1_loc_NHMMER.csv"
${scripts}/5_extract_ITS1_loc_from_nhmmer.py -f ${releases}/${release}/align_${release}\
                                             -s ${releases}/${release}/fasta_${release}\
                                             -p ${cpus}\
                                             -r ${release}\
                                             -o ${releases}/${release}/ITS1_loc_NHMMER.csv
exit_verify
echo ""
echo ""
echo "############################################################################################"
echo "${scripts}/6_extract_ITS1_loc_from_ENA.py
                                       -r ${release}
                                       -i ${releases}/${release}/tsv_${release}"
${scripts}/6_extract_ITS1_loc_from_ENA.py\
                            -r ${release}\
                            -i ${releases}/${release}/tsv_${release}
exit_verify
echo ""
echo ""
echo "############################################################################################"
echo "${scripts}/7_ITSoneDB_ITS1_summary.py
                                      -r ${release}
                                      -m ${releases}/${release}/ITS1_loc_NHMMER.csv
                                      -e ${releases}/${release}/ITS1_loc_ENA.csv
                                      -i ${releases}/${release}/index_${release}/index.csv.gz
                                      -o ${releases}/${release}/ITS1_loc_final_table.csv"
${scripts}/7_ITSoneDB_ITS1_summary.py -r ${release}\
                                      -m ${releases}/${release}/ITS1_loc_NHMMER.csv\
                                      -e ${releases}/${release}/ITS1_loc_ENA.csv\
                                      -i ${releases}/${release}/index_${release}/index.csv.gz\
                                      -o ${releases}/${release}/ITS1_loc_final_table.csv
exit_verify

fi