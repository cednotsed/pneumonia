#$ -l tmem=10G
#$ -l h_vmem=10G
#$ -l h_rt=480:0:0
#$ -wd /SAN/ugi/HAP_VAP/pneumonia/assembly_scripts
#$ -S /bin/bash
#$ -pe smp 8
#$ -j y
#$ -R y
#$ -N checkm

source /SAN/ballouxlab/uk_bats_meta/miniconda3/etc/profile.d/conda.sh
conda activate nanopore

wkdir=/SAN/ugi/HAP_VAP/pneumonia
#read_path=$1
#read_dir=$wkdir/data/basecalled_fastqs
in_dir=$wkdir/results/assembly_out/medaka_out/all_polished
#in_dir=$(echo $read_path| sed "s|$read_dir|$polish_dir|g"| sed "s|.fastq.gz||g")
#out_dir=$(echo $in_dir| sed "s|medaka_out|checkm_out|g")
out_dir=$wkdir/results/assembly_out/checkm_out

echo $in_dir
echo $out_dir

checkm lineage_wf -t 8 -x fasta $in_dir $out_dir
