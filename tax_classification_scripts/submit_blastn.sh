wkdir=/SAN/ugi/HAP_VAP/pneumonia
query_dir=$wkdir/results/assembly_out/medaka_out/all_polished
out_dir=$wkdir/results/metagenomic_out/blast_out

for i in $query_dir/*
do
  echo $i
  qsub blastn_contigs.sh "$i"
done
