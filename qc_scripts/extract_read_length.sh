in_dir=../data/basecalled_fastqs/no_humans
out_path=../results/qc_out/read_lengths.txt

for i in $in_dir/*.fastq.gz
do
    echo $i
    zcat $i | awk -v file="$i" '{if(NR%4==2) {print length($1),file}}'|sed "s|$in_dir/||g" >> $out_path
done
