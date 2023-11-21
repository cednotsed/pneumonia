#$ -l tmem=8G
#$ -l h_vmem=8G
#$ -l h_rt=480:0:0
#$ -wd /SAN/ugi/HAP_VAP/pneumonia/metagenomic_scripts/logs
#$ -S /bin/bash
#$ -pe smp 12
#$ -j y
#$ -R y
#$ -N kraken2

start=`date +%s`

source /SAN/ballouxlab/uk_bats_meta/miniconda3/etc/profile.d/conda.sh
conda activate metagenomics

WKDIR=/SAN/ugi/HAP_VAP/pneumonia
OUTDIR=$WKDIR/results/metagenomic_out/kraken2_out
FASTQ_DIR=$WKDIR/data/basecalled_fastqs
DB=$WKDIR/databases/k2_pluspf_20231009
TMP_DIR=$WKDIR/results/metagenomic_out/k2_temp
N_THREADS=12

mkdir -p $TMP_DIR
mkdir -p $OUTDIR

#Classification
for fq1 in $FASTQ_DIR/*.fastq.gz
do
    prefix=$(echo $fq1|sed "s|$FASTQ_DIR/||g"| sed "s|.fastq.gz||g")
    out_path=$OUTDIR/$prefix.tsv

    echo $out_path

    kraken2 \
        --threads $N_THREADS\
        --db $DB \
        --report $out_path \
        --output $TMP_DIR \
        $fq1 $fq2
done
