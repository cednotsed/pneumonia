#$ -l tmem=10G
#$ -l h_vmem=10G
#$ -l h_rt=480:0:0
#$ -wd /SAN/ugi/HAP_VAP/pneumonia/metagenomic_scripts/logs
#$ -S /bin/bash
#$ -pe smp 8
#$ -j y
#$ -R y
#$ -N assembly

# MODEL CHOSEN BASED ON https://github.com/epi2me-labs/wf-bacterial-genomes/blob/master/data/medaka_models.tsv
# guppy printed workflow FLO-MIN106 SQK-RPB004 dna_r9.4.1_450bps_hac
# Dorado model taken from guppy config is dna_r9.4.1_e8_hac@v3.3
# Corresponding MEDAKA model is r941_min_hac_g507
# BASED ON https://github.com/epi2me-labs/wf-bacterial-genomes/blob/master/data/medaka_models.tsv

#source /SAN/ballouxlab/uk_bats_meta/miniconda3/etc/profile.d/conda.sh
#conda activate nanopore

wkdir=/SAN/ugi/HAP_VAP/pneumonia
read_path=$1
read_dir=$wkdir/data/basecalled_fastqs
ass_dir=$wkdir/results/assembly_out/unicycler_out
in_dir=$(echo $read_path| sed "s|$read_dir|$ass_dir|g"| sed "s|.fastq.gz||g")
out_dir=$(echo $in_dir| sed "s|unicycler_out|medaka_out|g")

echo $in_dir
echo $out_dir

medaka_consensus \
  -i $read_path \
  -d $in_dir/assembly.fasta \
  -o $out_dir \
  -t 8 \
  -m r941_min_hac_g507

