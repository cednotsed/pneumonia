wkdir=/SAN/ugi/HAP_VAP/pneumonia
results_dir=$wkdir/results/amr_out/resfinder_out
out_dir=$wkdir/results/amr_out/resfinder_out/all_phenotypes

for i in $results_dir/INHALE*
do
  in_path=$i/pheno_table.txt
  out_path=$(echo $in_path| sed "s|/pheno_table.txt|.pheno_table.txt|g"| sed "s|$results_dir|$out_dir|g")
  echo $in_path
  echo $out_path
  cp $in_path $out_path
done
