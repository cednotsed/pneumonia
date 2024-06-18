wkdir=/flask/scratch/matthewsp/pneumonia
taxid_list=$wkdir/data/metadata/irep/references_to_download.txt
genome_basedir=$wkdir/data/genomes/irep

while read line
do
    taxid=$(echo $line|cut -d',' -f1)
    out_dir=$genome_basedir/$taxid
    mkdir $out_dir

    ncbi-genome-download \
        --taxids $taxid \
        --format fasta \
        --assembly-levels complete,chromosome \
        --flat-output \
        -o $out_dir \
        --parallel 8 \
        bacteria

done < $taxid_list
