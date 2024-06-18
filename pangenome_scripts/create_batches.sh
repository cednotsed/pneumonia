bug_prefix=m_catarrhalis
wkdir=/flask/scratch/matthewsp/pneumonia
fna_dir=$wkdir/data/genomes/hunt_et_al_v0.2/$bug_prefix.filt
meta_dir=$wkdir/data/metadata/bug_metadata/prokka_batches

mkdir $meta_dir/$bug_prefix

ls $fna_dir > $meta_dir/$bug_prefix.combined.txt

split -l 1000 $meta_dir/$bug_prefix.combined.txt $meta_dir/$bug_prefix/

