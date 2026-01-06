#!/bin/bash

echo -e "###########################################################################################" "\n"

echo -e ===== Identificación taxonómica de virus con BLASTn en ensambles obtenidos con SPAdes ===== "\n"

echo -e          "\t"                 ===== Inicio: $(date) ===== "\n"

echo -e "###########################################################################################" "\n"

#Crear la base de datos de BLASTn: makeblastdb -in archivo.fa -dbtype nucl -out ./virus_db
#La base de datos de virus se puede descargar en: https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=Influenza%20A%20virus,%20taxid:11320&HostLineage_ss=NOT%20Homo%20sapiens%20(human),%20taxid:HostId_i:*%20NOT%20HostId_i:9606

#---------------------------------------------------------
# Definir rutas de directorios de entrada y salida
dirfa="/home/user/Analisis_corridas/SPAdes/virus"
dirout="/home/user/Analisis_corridas/SPAdes/virus/BLASTn_results"
dirblas="/home/user/Analisis_corridas/SPAdes/virus/BLAST_assembly"
dirkmer="/home/user/Analisis_corridas/kmerfinder/virus"
diroutens="/home/user/Analisis_corridas/Resultados_all_virus/Ensambles"
db="virus_db"
#---------------------------------------------------------

cd ${dirfa}

echo -e "db = ${db}"

for ensamble in *.fa; do
    ID="$(basename ${ensamble} | cut -d '-' -f '1')"

# -------------------------------------------------------------------------
# Ejecutar BLASTn sobre los ensambles para identificar los contigs de virus
# -------------------------------------------------------------------------

blastn -query ${ensamble} -db $Bn_DB_PATH/${db} -outfmt "6 qseqid salltitles sstrand pident qcovs bitscore evalue" \
       -max_target_seqs 1 -max_hsps 1 -culling_limit 1 -perc_identity 80 -evalue 1e-10 \
       -out ${dirout}/${ID}_results.tsv

#Para conocer el % de identidad: -outfmt "6 pident"
#ID de secuencia de consulta: -outfmt "6 qseqid"
#Sentido de la secuencia: "sstrand"
#Imprimir e-value:  "evalue"
#Imprimir el bitscore: "bitscore"

# --------------------------------
# Modificar los archivos de salida
# --------------------------------

cat ./BLASTn_results/${ID}_results.tsv | tr " " "_" > ./BLASTn_results/${ID}_results_2.tsv
cat ./BLASTn_results/${ID}_results_2.tsv | awk '{print $1}' > ./BLASTn_results/${ID}_nodos.txt
cat ./BLASTn_results/${ID}_results_2.tsv | awk '{print $2}' > ./BLASTn_results/${ID}_gen.txt
cat ./BLASTn_results/${ID}_results_2.tsv | awk '{print $3}' > ./BLASTn_results/${ID}_sentido.txt
cat ./BLASTn_results/${ID}_results_2.tsv | awk '{print $4}' > ./BLASTn_results/${ID}_ident.txt
cat ./BLASTn_results/${ID}_results_2.tsv | awk '{print $5}' > ./BLASTn_results/${ID}_cov.txt
cat ./BLASTn_results/${ID}_results_2.tsv | awk '{print $6}' > ./BLASTn_results/${ID}_bitsc.txt
cat ./BLASTn_results/${ID}_results_2.tsv | awk '{print $7}' > ./BLASTn_results/${ID}_eval.txt
cat ./BLASTn_results/${ID}_gen.txt | cut -d ',' -f '1' > ./BLASTn_results/${ID}_gen1.txt
cat ./BLASTn_results/${ID}_gen1.txt | tr "_" " " > ./BLASTn_results/${ID}_gen2.txt
cat ./BLASTn_results/${ID}_gen2.txt | tr "( )" " | " > ./BLASTn_results/${ID}_gen3.txt
paste ./BLASTn_results/${ID}_nodos.txt ./BLASTn_results/${ID}_sentido.txt ./BLASTn_results/${ID}_ident.txt ./BLASTn_results/${ID}_cov.txt ./BLASTn_results/${ID}_bitsc.txt ./BLASTn_results/${ID}_eval.txt  ./BLASTn_results/${ID}_gen3.txt > ./BLASTn_results/${ID}_BLASTn_results_tmp.tsv
sed -i '1i Contig\tSentido\t%Identidad\t%Cobertura\tBitscore\te_value\tSec_ref' ./BLASTn_results/${ID}_BLASTn_results_tmp.tsv
cat ./BLASTn_results/${ID}_BLASTn_results_tmp.tsv | uniq > ./BLASTn_results/${ID}_BLASTn_results.tsv

# -------------------------
# Remover archivos sin uso
# -------------------------

rm ./BLASTn_results/${ID}_results.tsv
rm ./BLASTn_results/${ID}_results_2.tsv
rm ./BLASTn_results/*.txt*
rm ./BLASTn_results/*tmp*

done

# -----------------------------------------------------
# Invertir secuencias en minus a su complemento inverso
# -----------------------------------------------------

for BLAST in ${dirout}/*tsv; do
    ID=$(basename ${BLAST} | cut -d '_' -f '1')

for assembly in ${dirfa}/*fa; do
    IDa=$(basename ${assembly} | cut -d '-' -f '1')


awk '$2 == "minus" {print $1}' ${BLAST} > ./BLAST_assembly/${IDa}_minus_contigs.txt
awk '$2 == "plus" {print $1}' ${BLAST} > ./BLAST_assembly/${IDa}_plus_contigs.txt
if [[ ${ID} == ${IDa} ]]; then
echo -e "If control: ${ID} ${IDa}"
seqtk subseq ${assembly} ./BLAST_assembly/${IDa}_minus_contigs.txt | seqtk  seq -r > ${IDa}_metaSPAdes_plus.fa
seqtk subseq ${assembly} ./BLAST_assembly/${IDa}_plus_contigs.txt > ${IDa}_SPAdes_plus_contigs.fa

cat ${IDa}_SPAdes_plus_contigs.fa ${IDa}_metaSPAdes_plus.fa > ./BLAST_assembly/${ID}-metaSPAdes-assembly-plus.fasta
chmod -R 775 ./BLAST_assembly/${ID}-metaSPAdes-assembly-plus.fasta

		else
	continue
fi
done
done

rm *_metaSPAdes_plus.fa
rm *_SPAdes_plus_contigs.fa
rm ./BLAST_assembly/*_plus_contigs* 
rm ./BLAST_assembly/*_minus_contigs.txt
rm ${dirout}/*_BLASTn_results.tsv #Eliminar los primeros resulstados de BLASTn
#find /home/user/Analisis_corridas/SPAdes/virus/BLAST_assembly -type f -size 0 -exec rm -f {} \;

# ----------------------------------------------
# Correr BLASTn para confirmar las orientaciones
# ----------------------------------------------

echo -e "db = ${db}"

for ens in ${dirblas}/*fasta; do
    ID=$(basename ${ens} | cut -d '-' -f '1')

blastn -query ${ens} -db $Bn_DB_PATH/${db} -outfmt "6 qseqid salltitles sstrand pident qcovs bitscore evalue" \
       -max_target_seqs 1 -max_hsps 1 -culling_limit 1 -perc_identity 80 -evalue 1e-10 \
       -out ${dirblas}/${ID}_results.tsv

# --------------------------------
# Modificar los archivos de salida
# --------------------------------

cat ./BLAST_assembly/${ID}_results.tsv | tr " " "_" > ./BLAST_assembly/${ID}_results_2.tsv
cat ./BLAST_assembly/${ID}_results_2.tsv | awk '{print $1}' > ./BLAST_assembly/${ID}_nodos.txt
cat ./BLAST_assembly/${ID}_results_2.tsv | awk '{print $2}' > ./BLAST_assembly/${ID}_gen.txt
cat ./BLAST_assembly/${ID}_results_2.tsv | awk '{print $3}' > ./BLAST_assembly/${ID}_sentido.txt
cat ./BLAST_assembly/${ID}_results_2.tsv | awk '{print $4}' > ./BLAST_assembly/${ID}_ident.txt
cat ./BLAST_assembly/${ID}_results_2.tsv | awk '{print $5}' > ./BLAST_assembly/${ID}_cov.txt
cat ./BLAST_assembly/${ID}_results_2.tsv | awk '{print $6}' > ./BLAST_assembly/${ID}_bitsc.txt
cat ./BLAST_assembly/${ID}_results_2.tsv | awk '{print $7}' > ./BLAST_assembly/${ID}_eval.txt
cat ./BLAST_assembly/${ID}_gen.txt | cut -d ',' -f '1' > ./BLAST_assembly/${ID}_gen1.txt
cat ./BLAST_assembly/${ID}_gen1.txt | tr "_" " " > ./BLAST_assembly/${ID}_gen2.txt
cat ./BLAST_assembly/${ID}_gen2.txt | tr "( )" " | " > ./BLAST_assembly/${ID}_gen3.txt
paste ./BLAST_assembly/${ID}_nodos.txt ./BLAST_assembly/${ID}_sentido.txt ./BLAST_assembly/${ID}_ident.txt ./BLAST_assembly/${ID}_cov.txt ./BLAST_assembly/${ID}_bitsc.txt ./BLAST_assembly/${ID}_eval.txt ./BLAST_assembly/${ID}_gen3.txt > ./BLAST_assembly/${ID}_BLASTn_results_tmp.tsv
sed -i '1i Contig\tSentido\t%Identidad\t%Cobertura\tBitscore\te_value\tSec_ref' ./BLAST_assembly/${ID}_BLASTn_results_tmp.tsv
cat ./BLAST_assembly/${ID}_BLASTn_results_tmp.tsv | uniq > ./BLAST_assembly/${ID}_BLASTn_results.tsv

# -------------------------
# Remover archivos sin uso
# -------------------------

rm ./BLAST_assembly/${ID}_results.tsv
rm ./BLAST_assembly/${ID}_results_2.tsv
rm ./BLAST_assembly/*.txt*
rm ./BLAST_assembly/*tmp*

done

#------------------------------------------------------
# Concatenar los resultados de BLAST revisado en un solo archivo
#------------------------------------------------------

for file in ./BLAST_assembly/*BLASTn_results*; do
    ename=$(basename ${file} | cut -d '_' -f '1')
echo -e "\n########## ${ename} ########## \n$(cat ${file})"

   done >> ./BLAST_assembly/BLASTn_all_rev.tsv
rm ./BLAST_assembly/*_BLASTn_results.tsv

# -----------------------------------------------------------------------------
# Eliminar contigs pequeños para todos los generos virales menos para influenza
# -----------------------------------------------------------------------------
cd ${dirblas}

for file in ${dirkmer}/*spa; do
    genero=$(cat ${file} | sed -n '2p' | cut -d ' ' -f '2,3,4' | cut -d ',' -f '1'| tr ' ' '_')
    ID=$(basename ${file} | cut -d '_' -f '1')
echo -e "${genero}"

for assembly in *.fasta; do
    assembly_ID=$(basename ${assembly} | cut -d '-' -f '1')

if [[ ${ID} == ${assembly_ID} ]]; then
        echo -e "If control: ${ID} ${assembly_ID}"
if [[ ${genero} != "Influenza_A_virus" ]]; then
	echo -e "If control: ${genero}"

seqtk seq -L 50 ${assembly} > ${ID}-metaSPAdes-assembly-plus.fa

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
chmod -R 775 *-metaSPAdes-assembly-plus.fa
# ------------------------------------------------------------------------------------
# Mover los ensambles a una carpeta nombrada con el genero del organismo identificado
# ------------------------------------------------------------------------------------

cd ${dirkmer}

for file in *spa; do
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

rm ${dirkmer}/*spa

echo -e "############################################" "\n"
echo -e   "\t" ===== Fin: $(date) =====  "\n"
echo -e "############################################"  "\n"
