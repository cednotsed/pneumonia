wkdir=/SAN/ugi/HAP_VAP/pneumonia
in_dir=$wkdir/data/raw_fast5s
out_dir=$wkdir/data/basecalled_fastqs

for i in $in_dir/*
do
  #in_files=$(find $i -type f|grep fastq.gz|grep fastq_pass)
  barcode_dirs=$(echo $i/basecalled_fastqs/barcode*)
  run_name=$(echo $i| tr '/' '\n' | tail -n 1)
  echo $run_name

  for dir in $barcode_dirs
  do
    barcode_name=$(echo $dir|tr '/' '\n' | tail -n 1)
    out_path=$out_dir/$run_name-$barcode_name.fastq.gz

    # Merge fastq files
    find $dir -type f|grep fastq.gz|grep fastq_pass|xargs -I'{}' cat '{}' >> $out_path
  done
done
