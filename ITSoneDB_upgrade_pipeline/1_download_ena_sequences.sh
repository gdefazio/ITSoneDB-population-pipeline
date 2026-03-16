#!/bin/bash


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



echo $(pwd)
INCLUDE=""
for dv in STD GSS HTC HTG WGS TSA STS CON ; do
for group in FUN ENV HUM INV MAM VRT MUS PLN ROD UNC ; do

INCLUDE="--include='${dv}_${group}_*.dat.gz' ${INCLUDE}"
done
done

ecode=1

while [[ $ecode -ne 0 ]]; do
  bash -c "rsync -av ${INCLUDE} --exclude='*' rsync://ftp.ebi.ac.uk/pub/databases/ena/sequence/snapshot_latest/std/ /home3/gdefazio/ITSoneDB/ena_repo"
ecode=$?
done

#bash -c "rsync -av ${INCLUDE} --exclude='*' rsync://ftp.ebi.ac.uk/pub/databases/ena/sequence/snapshot_latest/expanded_con/ /home3/gdefazio/ITSoneDB/ena_repo"
