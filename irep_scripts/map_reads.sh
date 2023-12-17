wkdir=/mnt/c/git_repos/pneumonia
fastq_dir=/mnt/d/basecalled_fastqs
ref_dir=$wkdir/data/genomes/irep_filt
out_dir=$wkdir/results/irep_out/mapping_out

# Iterate through fastqs
for fq in $fastq_dir/*.fastq.gz
do
    for ref in $ref_dir/*
    do
        taxid=taxid$(echo $ref| sed "s|$ref_dir/||g")
        ref_name=$(ls $ref/*.fna.gz|sed "s|$ref/||g"| sed "s|.fna.gz||g")
        run_name=$(echo $fq|sed "s|$fastq_dir/||g"|sed "s|.fastq.gz||g")
        out_path=$out_dir/${run_name}.${ref_name}.${taxid}.sam
 #       echo $run_name $ref_name $taxid
        echo $out_path       

        minimap2 \
            -ax map-ont \
            -t 12 \
            $ref $fq > $out_path
    done
done
