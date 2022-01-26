#!/usr/bin/env bash

help="questo script serve per lanciare il download della release desiderata di ENA
in modo sequenziale senza l'intervento dell'operatore per tutte le categorie di sequenze
utili all'uprgade del database ITSoneDB. Il download viene effettuato nella cartella
in cui viene lanciato lo script.
usage ./download_ena_sequences.sh -r release di ENA -f cartella in cui scaricare
Options:
	-r    release: numero della release di ENA dalla quale si vuole avviare il download
	-d    folder: cartella nella quale scaricare
	-h    stampa questo help
"

if [[ $# -eq 0 ]]   			# se lo script è eseguito senza argomenti...
	then
  		echo -e "$help"		# ... stampo l'help...
  		exit 1         		# ... ed esco.
fi
while getopts ":h:r:d:" opzione
do
	case "${opzione}" in
	r) release=${OPTARG};;
	d) folder=${OPTARG};;
	h) echo -e "\n$help";exit 1;;
	#*) echo -e "\nÈ stata inserita un'opzione non valida.\n\n $help";exit 1;;
	esac
done

#DOWNLOAD TAXONOMY
# echo ${folder}
taxdir=`dirname ${folder}`/taxonomy_${release}
# echo ${taxdir}
# mkdir ${taxdir}
echo ""
echo "Download NCBI taxonomy at ${taxdir}/new_taxdump.tar.gz"
if [[ -a ${taxdir}/new_taxdump.tar.gz ]]
then
echo ''
else
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz -q -P ${taxdir}
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz -q -P ${taxdir}
tar -xf ${taxdir}/new_taxdump.tar.gz -C ${taxdir}
fi


echo "Checking if exists, create otherwise: ${folder}"
#wd=$(pwd)
if [[ -d ${folder} ]]
then
	echo "${folder} exists"
else
	mkdir -p ${folder}
	echo "The path ${folder} has just been totally or partially created"
fi

echo $(pwd)
for dv in std gss htc htg mga wgs tsa sts
do
for group in fun env hum inv mam vrt mus pln rod unc
do
echo "DOWNLOAD rel_${dv}_${group}_*_r${release}.dat.gz"
#     ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/std/rel_est_inv_23_r138.dat.gz
wget  ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/std/rel_${dv}_${group}_*_r${release}.dat.gz -q -P ${folder}
done
done