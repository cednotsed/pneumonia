WKDIR=/SAN/ugi/HAP_VAP/pneumonia
ASS_DIR=$WKDIR/results/assembly_out/medaka_out/all_polished

#Classification
for fna in $ASS_DIR/*.fasta
do
  qsub worker_diamond_blast.sh "$fna"
done


