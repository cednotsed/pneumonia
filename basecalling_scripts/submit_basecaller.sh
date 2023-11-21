base_dir=/SAN/ugi/HAP_VAP/pneumonia/data/raw_fast5s
for i in $base_dir/*
do
  echo $i
  qsub basecall_single_run.sh "$i"
done
