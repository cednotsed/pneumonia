wkdir=/mnt/c/git_repos/pneumonia
genome_basedir=$wkdir/data/genomes/irep
out_basedir=$wkdir/data/genomes/irep_filt

for dir in $genome_basedir/*
do
    taxid=$(echo $dir|sed "s|$genome_basedir/||g")
    out_dir=$out_basedir/$taxid
    echo $out_dir
    mkdir $out_dir   

    # Copy first reference only
    first_file=$(ls $dir/*.fna.gz|head -n 1)
    cp $first_file $out_dir
done
