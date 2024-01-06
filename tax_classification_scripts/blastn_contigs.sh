#$ -l tmem=5G
#$ -l h_vmem=5G
#$ -l h_rt=480:0:0
#$ -wd /SAN/ugi/HAP_VAP/pneumonia/metagenomic_scripts/logs
#$ -S /bin/bash
#$ -pe smp 8
#$ -j y
#$ -R y
#$ -N blastn

start=`date +%s`

source /SAN/ballouxlab/uk_bats_meta/miniconda3/etc/profile.d/conda.sh
conda activate metagenomics

wkdir=/SAN/ugi/HAP_VAP/pneumonia
query_dir=$wkdir/results/assembly_out/medaka_out/all_polished
out_dir=$wkdir/results/metagenomic_out/blast_out
db=$wkdir/databases/refseq_genomes_151223.blastn_db/refseq_genomes_151223.blastn_db

# Run blastn
query=$1
out=$(echo $query|sed "s|$query_dir|$out_dir|g"|sed "s|.polished.fasta|.tsv|g")

echo ${query}
echo ${out}
echo ${db}

blastn -db ${db} \
  -task dc-megablast \
  -query ${query} \
  -out ${out} \
  -qcov_hsp_perc 99 \
  -outfmt 6 \
  -num_threads 8

