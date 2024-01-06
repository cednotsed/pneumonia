#$ -l tmem=3G
#$ -l h_vmem=3G
#$ -l h_rt=480:0:0
#$ -wd /SAN/ugi/HAP_VAP/pneumonia/metagenomic_scripts/logs
#$ -S /bin/bash
#$ -pe smp 8
#$ -j y
#$ -R y
#$ -N diamond

start=`date +%s`

source /SAN/ballouxlab/uk_bats_meta/miniconda3/etc/profile.d/conda.sh
conda activate nanopore

WKDIR=/SAN/ugi/HAP_VAP/pneumonia
OUTDIR=$WKDIR/results/tax_classification_out/diamond_blast_out
ASS_DIR=$WKDIR/results/assembly_out/medaka_out/all_polished
DB=$WKDIR/databases/nr_131223/nr.dmnd
N_THREADS=8

#Classification
fna=$1

prefix=$(echo $fna|sed "s|$ASS_DIR/||g"| sed "s|.fasta||g"| sed "s|.polished||g")
out_path=$OUTDIR/$prefix.tsv

echo $prefix
echo $out_path

diamond blastx \
  -q $fna \
  -d $DB \
  -o $out_path \
  -F 15 \
  -f 6 \
  --range-culling \
  --top 10 \
  -p $N_THREADS
