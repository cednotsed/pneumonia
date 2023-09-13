#$ -l tmem=10G
#$ -l h_rt=480:0:0
#$ -wd /SAN/ugi/HAP_VAP/pneumonia
#$ -S /bin/bash
#$ -j y
#$ -pe gpu 1
#$ -R y
#$ -N basecalling

base_dir=/SAN/ugi/HAP_VAP/pneumonia/data/raw_fast5s
threads=6

for i in $base_dir/*.tar
do
  echo $i
  #i=/SAN/ugi/HAP_VAP/pneumonia/data/raw_fast5s/INHALE_FRESH_14.tar
  fast5_dir=$(echo $i|sed "s|.tar||g")
  echo $fast5_dir

  # Untar
  echo UNTARRING...
  mkdir $fast5_dir
  tar -xf $i -C $fast5_dir --strip-components=1

  # Shift fast5s
  echo PROCESSING FAST5s
  temp_dir=$fast5_dir/all_fast5s
  mkdir $temp_dir
  #find $fast5_dir -type f|grep \.fast5|xargs -P6 -I'{}' mv '{}' $temp_dir
  mv $fast5_dir/* $temp_dir

  # Make multi fast5s
  echo concatenate fast5s
  in_dir=$fast5_dir/multi_fast5s
  mkdir $in_dir

  single_to_multi_fast5 \
    --input_path $temp_dir \
    --save_path $in_dir \
    --threads $threads \
    --recursive \
    --batch_size 4096

  # Clear single fast5s
  find $temp_dir -type f| grep \.fast5| xargs -P $threads -I '{}' rm '{}'
  rm -r $temp_dir

  # Basecalling
  echo BASECALLING...
  out_dir=$fast5_dir/basecalled_fastqs
  /SAN/ugi/HAP_VAP/ont-guppy/bin/guppy_basecaller \
    --device "cuda:0" \
    --input_path $in_dir \
    --save_path $out_dir \
    --config dna_r9.4.1_450bps_hac.cfg \
    --barcode_kits "SQK-RBK004" \
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

  #rm -r $in_dir

  echo DONE $i
done

echo ALL FILES COMPLETED!!!
