#!/bin/bash

echo -e "###########################################################################################" "\n"

echo -e ===== Identificación taxonómica de virus con BLASTn en ensambles obtenidos con SPAdes ===== "\n"

echo -e          "\t"                 ===== Inicio: $(date) ===== "\n"

echo -e "###########################################################################################" "\n"

#Crear la base de datos de BLASTn: makeblastdb -in archivo.fa -dbtype nucl -out ./virus_db
#La base de datos de virus se puede descargar en: https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=Influenza%20A%20virus,%20taxid:11320&HostLineage_ss=NOT%20Homo%20sapiens%20(human),%20taxid:HostId_i:*%20NOT%20HostId_i:9606

#---------------------------------------------------------
dirfa="/home/admcenasa/Analisis_corridas/SPAdes/virus"
dirout="/home/admcenasa/Analisis_corridas/SPAdes/virus/BLASTn_results"
dirblas="/home/admcenasa/Analisis_corridas/SPAdes/virus/BLAST_assembly"
dirkmer="/home/admcenasa/Analisis_corridas/kmerfinder/virus"
diroutens="/home/admcenasa/Analisis_corridas/Resultados_all_virus/Ensambles"
#---------------------------------------------------------

cd ${dirfa}

db="virus_db"

echo -e "db = ${db}"

for ensamble in *.fa; do
    ID="$(basename ${ensamble} | cut -d '-' -f '1')"

# ------------------------------------------------------------------------

blastn -query ${ensamble} -db $BnV_DB_PATH/${db} \
       -outfmt "6 qseqid salltitles sstrand pident qcovs bitscore evalue" \
       -max_hsps 1 -culling_limit 1 \
       -perc_identity 80 -evalue 1e-10 \
       -out ${dirout}/${ID}_results.tsv

# --------------------------------

cat ${dirout}/${ID}_results.tsv | tr " " "_" > ${dirout}/${ID}_results_2.tsv
cat ${dirout}/${ID}_results_2.tsv | awk '{print $1}' > ${dirout}/${ID}_nodos.txt
cat ${dirout}/${ID}_results_2.tsv | awk '{print $2}' > ${dirout}/${ID}_gen.txt
cat ${dirout}/${ID}_results_2.tsv | awk '{print $3}' > ${dirout}/${ID}_sentido.txt
cat ${dirout}/${ID}_results_2.tsv | awk '{print $4}' > ${dirout}/${ID}_ident.txt
cat ${dirout}/${ID}_results_2.tsv | awk '{print $5}' > ${dirout}/${ID}_cov.txt
cat ${dirout}/${ID}_results_2.tsv | awk '{print $6}' > ${dirout}/${ID}_bitsc.txt
cat ${dirout}/${ID}_results_2.tsv | awk '{print $7}' > ${dirout}/${ID}_eval.txt
cat ${dirout}/${ID}_gen.txt | cut -d ',' -f '1' > ${dirout}/${ID}_gen1.txt
cat ${dirout}/${ID}_gen1.txt | tr "_" " " > ${dirout}/${ID}_gen2.txt
cat ${dirout}/${ID}_gen2.txt | tr "( )" " | " > ${dirout}/${ID}_gen3.txt
paste ${dirout}/${ID}_nodos.txt ${dirout}/${ID}_sentido.txt ${dirout}/${ID}_ident.txt ${dirout}/${ID}_cov.txt ${dirout}/${ID}_bitsc.txt ${dirout}/${ID}_eval.txt  ${dirout}/${ID}_gen3.txt > ${dirout}/${ID}_BLASTn_results_tmp.tsv
sed -i '1i Contig\tSentido\t%Identidad\t%Cobertura\tBitscore\te_value\tSec_ref' ${dirout}/${ID}_BLASTn_results_tmp.tsv
cat ${dirout}/${ID}_BLASTn_results_tmp.tsv | uniq > ${dirout}/${ID}_BLASTn_results.tsv

# -------------------------

rm ${dirout}/${ID}_results.tsv
rm ${dirout}/${ID}_results_2.tsv
rm ${dirout}/*.txt*
rm ${dirout}/*tmp*

	done

# -----------------------------------------------------

for BLAST in ${dirout}/*tsv; do
    ID=$(basename ${BLAST} | cut -d '_' -f '1')

for assembly in ${dirfa}/*fa; do
    IDa=$(basename ${assembly} | cut -d '-' -f '1')


awk '$2 == "minus" {print $1}' ${BLAST} > ${dirblas}/${IDa}_minus_contigs.txt
awk '$2 == "plus" {print $1}' ${BLAST} > ${dirblas}/${IDa}_plus_contigs.txt

	if [[ ${ID} == ${IDa} ]]; then

echo -e "If control: ${ID} ${IDa}"

#-----

seqtk subseq ${assembly} ${dirblas}/${IDa}_minus_contigs.txt | seqtk  seq -r > ${IDa}_metaSPAdes_plus.fa
seqtk subseq ${assembly} ${dirblas}/${IDa}_plus_contigs.txt > ${IDa}_SPAdes_plus_contigs.fa

#-----

cat ${IDa}_SPAdes_plus_contigs.fa ${IDa}_metaSPAdes_plus.fa > ${dirblas}/${ID}-metaSPAdes-assembly-plus_tmp.fasta

#-----

seqkit sort -l \
            -r ${dirblas}/${ID}-metaSPAdes-assembly-plus_tmp.fasta > ${dirblas}/${ID}-metaSPAdes-assembly-plus.fasta


		else
         continue
       fi
    done
done

rm *_metaSPAdes_plus.fa
rm *_SPAdes_plus_contigs.fa
rm ${dirblas}/*_plus_contigs*
rm ${dirblas}/*_minus_contigs.txt
rm ${dirout}/*_BLASTn_results.tsv
rm ${dirblas}/*-assembly-plus_tmp.fasta
#find /home/admcenasa/Analisis_corridas/SPAdes/virus/BLAST_assembly -type f -size 0 -exec rm -f {} \;

# ----------------------------------------------

echo -e "db = ${db}"

for ens in ${dirblas}/*fasta; do
    ID=$(basename ${ens} | cut -d '-' -f '1')

blastn -query ${ens} -db $BnV_DB_PATH/${db} \
       -outfmt "6 qseqid salltitles sstrand pident qcovs bitscore evalue" \
       -max_hsps 1 -culling_limit 1 \
       -perc_identity 80 -evalue 1e-10 \
       -out ${dirblas}/${ID}_results.tsv

# --------------------------------

cat ${dirblas}/${ID}_results.tsv | tr " " "_" > ${dirblas}/${ID}_results_2.tsv
cat ${dirblas}/${ID}_results_2.tsv | awk '{print $1}' > ${dirblas}/${ID}_nodos.txt
cat ${dirblas}/${ID}_results_2.tsv | awk '{print $2}' > ${dirblas}/${ID}_gen.txt
cat ${dirblas}/${ID}_results_2.tsv | awk '{print $3}' > ${dirblas}/${ID}_sentido.txt
cat ${dirblas}/${ID}_results_2.tsv | awk '{print $4}' > ${dirblas}/${ID}_ident.txt
cat ${dirblas}/${ID}_results_2.tsv | awk '{print $5}' > ${dirblas}/${ID}_cov.txt
cat ${dirblas}/${ID}_results_2.tsv | awk '{print $6}' > ${dirblas}/${ID}_bitsc.txt
cat ${dirblas}/${ID}_results_2.tsv | awk '{print $7}' > ${dirblas}/${ID}_eval.txt
cat ${dirblas}/${ID}_gen.txt | cut -d ',' -f '1' > ${dirblas}/${ID}_gen1.txt
cat ${dirblas}/${ID}_gen1.txt | tr "_" " " > ${dirblas}/${ID}_gen2.txt
cat ${dirblas}/${ID}_gen2.txt | tr "( )" " | " > ${dirblas}/${ID}_gen3.txt
paste ${dirblas}/${ID}_nodos.txt ${dirblas}/${ID}_sentido.txt ${dirblas}/${ID}_ident.txt ${dirblas}/${ID}_cov.txt ${dirblas}/${ID}_bitsc.txt ${dirblas}/${ID}_eval.txt ${dirblas}/${ID}_gen3.txt > ${dirblas}/${ID}_BLASTn_results_tmp.tsv
sed -i '1i Contig\tSentido\t%Identidad\t%Cobertura\tBitscore\te_value\tSec_ref' ${dirblas}/${ID}_BLASTn_results_tmp.tsv
cat ${dirblas}/${ID}_BLASTn_results_tmp.tsv | uniq > ${dirblas}/${ID}_BLASTn_results.tsv

# -------------------------

rm ${dirblas}/${ID}_results.tsv
rm ${dirblas}/${ID}_results_2.tsv
rm ${dirblas}/*.txt*
rm ${dirblas}/*tmp*

done

#------------------------------------------------------

for file in ${dirblas}/*BLASTn_results*; do
    ename=$(basename ${file} | cut -d '_' -f '1')
echo -e "\n########## ${ename} ########## \n$(cat ${file})"

   done >> ${dirblas}/BLASTn_all_rev.tsv
rm ${dirblas}/*_BLASTn_results.tsv

# -----------------------------------------------------------------------------

cd ${dirblas}

for file in ${dirkmer}/*.spa; do
    genero=$(cat ${file} | sed -n '2p' | cut -d ' ' -f '2,3,4' | cut -d ',' -f '1'| tr ' ' '_')
    ID=$(basename ${file} | cut -d '_' -f '1')
echo -e "${genero}"

for assembly in *.fasta; do
    assembly_ID=$(basename ${assembly} | cut -d '-' -f '1')

if [[ ${ID} == ${assembly_ID} ]]; then
        echo -e "If control: ${ID} ${assembly_ID}"
if [[ ${genero} != "Influenza_A_virus" ]]; then
	echo -e "If control: ${genero}"

seqtk seq -L 100 ${assembly} > ${ID}-metaSPAdes-assembly-plus.fa

	else

echo -e "else control: ${genero}"
mv ${assembly} ${ID}-metaSPAdes-assembly-plus.fa
echo -e "mv ${assembly} a ${ID}-metaSPAdes-assembly-plus.fa"

	continue

echo -e "Else control: ${genero}"

	fi
      fi
  done
done

rm *.fasta

# ------------------------------------------------------------------------------------

cd ${dirkmer}

for file in *.spa; do
    genero=$(cat ${file} | sed -n '2p' | cut -d ' ' -f '2,3,4' | cut -d ',' -f '1'| tr ' ' '_')
    ID=$(basename ${file} | cut -d '_' -f '1')

for assembly in ${dirblas}/*.fa; do
    assembly_ID=$(basename ${assembly} | cut -d '-' -f '1')

if [[ ${ID} != ${assembly_ID} ]]; then

       continue

 else

mkdir -p ${diroutens}/${genero}

echo -e "Moviendo ${assembly} a ${genero}"
     mv ${assembly} ${diroutens}/${genero}

        fi
    done
done

rm ${dirkmer}/*.spa

echo -e "############################################" "\n"
echo -e   "\t" ===== Fin: $(date) =====  "\n"
echo -e "############################################"  "\n"
