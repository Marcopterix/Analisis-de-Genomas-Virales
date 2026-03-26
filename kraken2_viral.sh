#!/bin/bash

echo -e "######################################################################################################################" "\n"
echo -e  ========= Determinación de abundancia de especies con Bracken en secuencias virales ========== "\n"
echo -e "\t" ===== Inicio: $(date) ===== "\n"
echo -e "######################################################################################################################" "\n"

#---------------------------------------------------------
dirfq="$HOME/Analisis_corridas/Archivos_postrim/virus"
dirout="$HOME/Analisis_corridas/kraken2/virus"
dirbrackout="$HOME/Analisis_corridas/Bracken/virus"
#---------------------------------------------------------

cd ${dirfq}

for R1 in *_R1_*; do
    R2=${R1/_R1_/_R2_}
    ID=$(basename ${R1} | cut -d '_' -f '1')

# ------------------------------------------------------------------------

echo -e "########## ${ID} ##########"

kraken2 --paired ${R1} ${R2} \
        --gzip-compressed \
        --db $K2_DB_PATH \
        --use-names \
        --threads 25 \
        --report ${dirout}/${ID}_K2report.txt > ${dirout}/${ID}_kraken2_out.txt

	done

rm ${dirout}/*kraken2_out*

# ----------------------------------------------------

cd ${dirout}

for K2 in *K2report*; do
    ID=$(basename ${K2} | cut -d '_' -f '1')

bracken -d $K2_DB_PATH \
        -i ${K2} \
        -o ${dirbrackout}/${ID}_bracken_report.tsv \
        -t 15

	done

#-------------------------------------------------------

for k2 in *K2report*; do
    ID=$(basename ${k2} | cut -d '_' -f '1')

ImportTaxonomy.pl -t 5 -s 1 -m 3 ${k2} -o ${dirbrackout}/${ID}_kraken2_krona.html

	done

rm ${dirout}/*K2report.txt
rm ${dirout}/*K2report_bracken_species.txt

echo -e "###########################################################################################" "\n"
echo -e  ========== Análisis de abundancia de especies terminado ========== "\n"
echo -e "\n" ===== Fin: $(date) ===== "\n"
echo -e "###########################################################################################"
