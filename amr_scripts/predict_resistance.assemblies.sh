#$ -l tmem=3G
#$ -l h_vmem=3G
#$ -l h_rt=480:0:0
#$ -wd /SAN/ugi/HAP_VAP/pneumonia/metagenomic_scripts/logs
#$ -S /bin/bash
#$ -pe smp 8
#$ -j y
#$ -R y
#$ -N resfinder

source /SAN/ballouxlab/uk_bats_meta/miniconda3/etc/profile.d/conda.sh
conda activate resfinder

wkdir=/SAN/ugi/HAP_VAP/pneumonia
in_dir=$wkdir/results/assembly_out/medaka_out/all_polished
in_fna=$in_dir/INHALE_FRESH_31-barcode05.polished.fasta
db_res=$wkdir/databases/resfinder_151223/resfinder_db
out_basedir=$wkdir/results/amr_out/resfinder_out

echo $in_dir
echo $out_dir

for in_fna in $in_dir/*.fasta
do
  out_dir=$out_basedir/$(echo $in_fna|sed "s|$in_dir/||g"| sed "s|.polished.fasta||g")

  echo $in_fna
  echo $out_dir

  run_resfinder.py \
    -ifa $in_fna \
    -o $out_dir \
    -s "Other" \
    -db_res $db_res \
    --nanopore \
    -acq \
    --min_cov 0.5 \
    --threshold 0.9
done


