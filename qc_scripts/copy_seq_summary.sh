base_dir=../data
read_dir=$base_dir/raw_fast5s
out_dir=$base_dir/sequencing_summaries
mkdir $out_dir

for in_path in $(ls $read_dir/INHALE_FRESH_15/basecalled_fastqs_new/sequencing_summary.txt)
do
  echo $in_path
  out_path=$(echo $in_path| sed "s|$read_dir|$out_dir|g"| sed "s|/basecalled_fastqs_new/|_|g")
  echo $out_path
  cp $in_path $out_path
done


