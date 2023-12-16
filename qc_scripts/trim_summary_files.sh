sum_dir=../data/sequencing_summaries

for sum_path in $sum_dir/INHALE_FRESH_15*.txt
do
  echo $sum_path
  out_path=$(echo $sum_path|sed "s|.txt|.parsed.csv|g")
  cat $sum_path |awk '{print $2,$11,$15,$16,$22}' > $out_path
done

