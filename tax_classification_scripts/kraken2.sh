#$ -l tmem=10G
#$ -l h_vmem=10G
#$ -l h_rt=480:0:0
#$ -wd /SAN/ugi/HAP_VAP/pneumonia/metagenomic_scripts/logs
#$ -S /bin/bash
#$ -pe smp 8
#$ -j y
#$ -R y
#$ -N kraken2

start=`date +%s`

source /SAN/ballouxlab/uk_bats_meta/miniconda3/etc/profile.d/conda.sh
conda activate metagenomics

WKDIR=/SAN/ugi/HAP_VAP/pneumonia
RUN_LIST=$WKDIR/data/raw_fast5s/k2_to_do.txt
OUTDIR=$WKDIR/results/metagenomic_out/kraken2_out
FASTQ_DIR=$WKDIR/data/basecalled_fastqs
DB=$WKDIR/databases/k2_pluspf_20231009
TMP_DIR=$WKDIR/results/metagenomic_out/k2_temp
N_THREADS=8

mkdir -p $TMP_DIR
mkdir -p $OUTDIR

#Classification
for fq1 in $FASTQ_DIR/INHALE_FRESH_15*.fastq.gz
do
  prefix=$(echo $fq1|sed "s|$FASTQ_DIR/||g"| sed "s|.fastq.gz||g")
  out_path=$OUTDIR/$prefix.tsv

  echo $out_path

  kraken2 \
      --threads $N_THREADS\
      --db $DB \
      --report $out_path \
      $fq1 > /dev/null
done
