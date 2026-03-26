#!/bin/bash

echo -e "#################################################################################################" "\n"
echo -e ===== Ejecutando SpeciesFinder sobre ensambles para la identificación taxonómica de virus  ===== "\n"
echo -e "\t" ===== Inicio: $(date) ===== "\n"
echo -e "#################################################################################################" "\n"

#-------------------------------------
dirfa="/home/admcenasa/Analisis_corridas/SPAdes/virus"
dirout="/home/admcenasa/Analisis_corridas/speciesfinder/virus"
dirgen="/home/admcenasa/Analisis_corridas/Resultados_all_virus/Ensambles"
dirfq="/home/admcenasa/Analisis_corridas/Archivos_postrim/virus/bowtie_filter"
diroutfq="/home/admcenasa/Analisis_corridas/Resultados_all_virus/Archivos_trimming"
#-------------------------------------

cd ${dirfa}

for assembly in *.fa; do
    ID=$(basename ${assembly} | cut -d '-' -f '1')

echo -e "########## ${ID} ##########"

speciesfinder -i ${assembly} \
              -o ${dirout}/SpF_${ID} \
              -db ${SPF_DB_PATH}/virus/virus_14_AT_110126 \
              -x \
              -tax ${SPF_DB_PATH}/virus/virus_14_AT_110126.tax

mv ${dirout}/SpF_${ID}/results.res ${dirout}/${ID}_results.res
rm -R ${dirout}/SpF_${ID}

	done

cd ${dirout}

for file in *.res; do
    tax=$(cat ${file} | sed -n '2p' | cut -d ' ' -f '2,3' | tr ' ' '_' )
    ID=$(basename ${file} | cut -d '_' -f '1')

for trimm in ${dirfq}/*.fastq.gz; do
    IDt=$(basename ${trimm} | cut -d '_' -f '1')

if [[ ${ID} != ${IDt} ]]; then

	continue

   else

mkdir -p ${diroutfq}/${tax}
echo -e "Moviendo ${trimm} a ${tax}"
mv -n ${trimm} ${diroutfq}/${tax}

	fi

   done
 done

for res in *.res; do
    ID=$(basename ${res} | cut -d '_' -f '1')

echo -e "\n ##### ${ID} ##### \n$(cat ${res})"

	done >> SpeciesFinder_viral_results.tsv

rm ${dirfa}/*.log
