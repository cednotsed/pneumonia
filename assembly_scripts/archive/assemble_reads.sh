#$ -l tmem=10G
#$ -l h_vmem=10G
#$ -l h_rt=480:0:0
#$ -wd /SAN/ugi/HAP_VAP/pneumonia/assembly_scripts/logs
#$ -S /bin/bash
#$ -pe smp 8
#$ -j y
#$ -R y
#$ -N assembly

start=`date +%s`

source /SAN/ballouxlab/uk_bats_meta/miniconda3/etc/profile.d/conda.sh
conda activate nanopore

wkdir=/SAN/ugi/HAP_VAP/pneumonia
in_dir=$wkdir/data/basecalled_fastqs
out_dir=$wkdir/results/assembly_out/unicycler_out

input=$1
output=$(echo $input| sed "s|$in_dir|$out_dir|g"| sed "s|.fastq.gz||g")

echo $input
echo $output

unicycler \
  -l $input \
  -o $output \
  -t 8 \
  --keep 3 \
  --mode normal \
  --verbosity 1
