#!/bin/bash

echo -e "###########################################################################################" "\n"

echo -e ===== Identificación del gen HA de IA con BLASTx en ensambles obtenidos con SPAdes ===== "\n"

echo -e          "\t"                 ===== Inicio: $(date) ===== "\n"

echo -e "###########################################################################################" "\n"

#Crear la base de datos de BLASTn: makeblastdb -in archivo.fa -dbtype prot -out ./virus_prot
#La base de datos de virus se puede descargar en: https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=Influenza%20A%20virus,%20taxid:11320&HostLineage_ss=NOT%20Homo%20sapiens%20(human),%20taxid:HostId_i:*%20NOT%20HostId_i:9606

#------------------------------------------------------------------------
# Definir rutas como variables
dirfa="/home/bioinfocenasa/Analisis_corridas/SPAdes/virus"
dirout="/home/bioinfocenasa/Analisis_corridas/BLAST/virus/BLAST_assembly"
db="IA_H_prot"
#--------------------------------------------------------------------------

cd ${dirfa}

echo -e "db = ${db}"

for ensamble in *.fa; do
    ID=$(basename ${ensamble} | cut -d '-' -f '1')

# ---------------------------------------------------------------------------
# Ejecutar BLASTX sobre los ensambles para identificar los contigs de interés
# ---------------------------------------------------------------------------

echo -e "\n\033[1;36m========== Searching for gene HA ==========\033[0m\n"

blastx -query ${ensamble} \
       -db $BxIA_DB_PATH/${db} \
       -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend evalue bitscore qseq" \
       -num_threads 5 \
	   -max_target_seqs 1 -max_hsps 1 \
       -out ${dirout}/${ID}_results.tsv

awk '{print ">"$1"_HA""\n"$11}' ${dirout}/${ID}_results.tsv > ${dirout}/${ID}_HA_prot.fna

# -----------------------------
# Información del alineamiento
# -----------------------------

awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$9"\t"$10"\t""HA"}' ${dirout}/${ID}_results.tsv \
> ${dirout}/${ID}_HA_info.tsv

sed -i '1i Contig\tReferencia\tIdentidad\tAling_long\te-value\tbitscore\tGen' ${dirout}/${ID}_HA_info.tsv

	done

rm ${dirout}/*_results.tsv


echo -e "############################################" "\n"
echo -e   "\t" ===== Fin: $(date) =====  "\n"
echo -e "############################################"  "\n"
