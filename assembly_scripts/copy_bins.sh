wkdir=/flask/scratch/matthewsp/pneumonia
base_dir=$wkdir/results/assembly_out/medaka_out
out_dir=$wkdir/results/assembly_out/bin_out

find $base_dir -type f|grep bin_dir > temp.out


while read line
do
    cp $line $(echo $line| sed "s|.no_human/bin_dir/|-|g"| sed "s|.fa|.fna|g"| sed "s|$base_dir|$out_dir|g")

done < temp.out
