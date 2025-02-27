wkdir=/mnt/c/git_repos/pneumonia
taxid_list=$wkdir/data/metadata/irep/taxids_to_download.txt
genome_basedir=$wkdir/data/genomes/irep_references

while read line
do
    taxid=$(echo $line|cut -d',' -f1)
    out_dir=$genome_basedir/$taxid
    mkdir $out_dir

    ncbi-genome-download \
        --taxids $taxid \
        --format fasta \
        --refseq-categories reference \
        --flat-output \
        -o $out_dir \
        --parallel 8 \
        bacteria

done < $taxid_list
