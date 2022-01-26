#!/bin/bash
# Master script per l'aggiornamento periodico di ITSoneDB
# L'aggiornamento prende in considerazione le divisioni pln e inv per intero
# versione modificata da Bruno Fosso
#lancia nella cartella in cui sono i fasta

help="Script per l'avvio del programma nhmmer in maniera ricursiva, sulle sequenze
fasta. Il programma ha bisogno di conoscere il numero della release a cui sta
lavorando.
    -r nr della release
    -t numero delle cpu da utilizzare
    -i input directory
    -g global aux dir path
    "

if [ $# -eq 0 ]   			# se lo script è eseguito senza argomenti...
	then
  		echo -e "$help"		# ... stampo l'help...
  		exit 1         		# ... ed esco.
fi
while getopts ":h:r:t:i:g:" opzione
do
	case "${opzione}" in
	r) release=${OPTARG};;
	t) thread=${OPTARG};;
	i) inputdir=${OPTARG};;
	g) gauxdp=${OPTARG};;
	h) echo -e "\n$help";exit 1;;
	#*) echo -e "\nÈ stata inserita un'opzione non valida.\n\n $help";exit 1;;
	esac
done

main_folder=`pwd`
#siccome viene lanciato in fasta la mainfolder sarà fasta
out_dir=`dirname ${inputdir}`/align_${release}
mkdir ${out_dir}
#if [[ ${thread} == "" ]]
#then
#thread=40
#fi


for compressed_fasta in $(ls ${inputdir})
    do
    fasta=$(basename ${compressed_fasta} .gz) #toglie il suffisso .gz dal nome del fasta
    folder=$(basename ${fasta} .fasta)        #toglie il suffisso .fasta dal nome del fasta
    mkdir ${out_dir}/${folder}
    echo ${gauxdp}/solo_fungi_5_8S.hmm
    nhmmer -o ${out_dir}/${folder}/align_5_8S_data.txt --cpu ${thread} -T 5 ${gauxdp}/solo_fungi_5_8S.hmm ${inputdir}/${compressed_fasta}
    nhmmer -o ${out_dir}/${folder}/align_18S_data.txt  --cpu ${thread} -T 5 ${gauxdp}/total_rfam_18S.hmm ${inputdir}/${compressed_fasta}
    #cd $main_folder
    gzip ${out_dir}/${folder}/align_5_8S_data.txt
    gzip ${out_dir}/${folder}/align_18S_data.txt
    done
