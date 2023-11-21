#$ -l tmem=10G
#$ -l h_rt=480:0:0
#$ -wd /SAN/ugi/HAP_VAP/pneumonia
#$ -S /bin/bash
#$ -j y
#$ -pe gpu 1
#$ -R y
#$ -N basecalling

base_dir=/SAN/ugi/HAP_VAP/pneumonia/data/raw_fast5s
threads=12

fast5_dir=$1
echo $fast5_dir

in_dir=$fast5_dir/multi_fast5s
out_dir=$fast5_dir/basecalled_fastqs_new

mkdir $out_dir

# Basecalling
echo BASECALLING...
/SAN/ugi/HAP_VAP/ont-guppy/bin/guppy_basecaller \
  --device "cuda:0" \
  --input_path $in_dir \
  --save_path $out_dir \
  --config dna_r9.4.1_450bps_hac.cfg \
  --barcode_kits "SQK-RPB004" \
  --detect_adapter \
  --detect_barcodes \
  --detect_primer \
  --trim_adapters \
  --trim_primers \
  --enable_trim_barcodes \
  --num_alignment_threads $threads \
  --barcode_nested_output_folder \
  --read_batch_size 64 \
  --compress_fastq \
  --num_callers 14 \
  --num_barcoding_threads $threads \
  --num_read_splitting_threads $threads \
  --gpu_runners_per_device 8 \
  --chunks_per_runner 768 \
  --chunk_size 1024 \
  --do_read_splitting \
  --detect_mid_strand_adapter \
  --detect_mid_strand_barcodes \
  --min_score_adapter 50 \
  --min_score_adapter_mid 50 \
  --min_score_barcode_front 50 \
  --min_score_barcode_mid 50 \
  --min_score_primer 50 \
  --min_score_read_splitting 50 \
  --max_read_split_depth 5

echo ALL FILES COMPLETED!!!
