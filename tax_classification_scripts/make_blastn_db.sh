#$ -l tmem=5G
#$ -l h_vmem=5G
#$ -l h_rt=480:0:0
#$ -wd /SAN/ugi/HAP_VAP/pneumonia/metagenomic_scripts/logs
#$ -S /bin/bash
#$ -pe smp 8
#$ -j y
#$ -R y
#$ -N make_blastdb

start=`date +%s`

source /SAN/ballouxlab/uk_bats_meta/miniconda3/etc/profile.d/conda.sh
conda activate metagenomics

wkdir=/SAN/ugi/HAP_VAP/pneumonia
query_dir=$wkdir/results/assembly_out/medaka_out/all_polished
out_dir=$wkdir/results/metagenomic_out/blast_out

## Make db ##
ref=$wkdir/data/genomes/refseq_genomes_151223/refseq_genomes_151223.fna
db=$wkdir/databases/refseq_genomes_151223.blastn_db/refseq_genomes_151223.blastn_db

makeblastdb -input_type fasta \
	-in $ref \
  -out $db \
	-parse_seqids \
	-title "refseq_genomes_151223.blastn_db" \
	-dbtype nucl
