wkdir=/mnt/c/git_repos/pneumonia
ref_dir=$wkdir/data/genomes/irep_filt
map_dir=$wkdir/results/irep_out/mapping_out

# Iterate through fastqs
for ref in $ref_dir/*
do
    taxid=taxid$(echo $ref| sed "s|$ref_dir/||g")
    ref_name=$(ls $ref/*.fna.gz|sed "s|$ref/||g"| sed "s|.fna.gz||g")
    run_name=$(echo $fq|sed "s|$fastq_dir/||g"|sed "s|.fastq.gz||g")
    out_path=$out_dir/${run_name}.${ref_name}.${taxid}.sam
    echo $run_name $ref_name $taxid
done
