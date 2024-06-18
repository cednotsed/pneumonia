wkdir=/flask/scratch/matthewsp/pneumonia
base_dir=$wkdir/results/assembly_out/medaka_out
out_dir=$wkdir/results/assembly_out/contig_out

for i in $base_dir/*
do
    echo $i
    cp $i/consensus.fasta $(echo $i/consensus.fasta| sed "s|$base_dir|$out_dir|g"| sed "s|.no_human/consensus.fasta|.fna|g")
done
