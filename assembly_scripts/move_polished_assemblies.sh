wkdir=/SAN/ugi/HAP_VAP/pneumonia
polish_dir=$wkdir/results/assembly_out/medaka_out

for i in $polish_dir/INHALE*
do
  out_path=$(echo $i/consensus.fasta| sed "s|/consensus|.polished|g"| sed "s|medaka_out|medaka_out/all_polished|g")
  cp $i/consensus.fasta $out_path
done
