file_dir=../data/genomes/irep_filt

find $file_dir -type f|sed "s|../data/genomes/irep_filt/||g" > reference_paths.txt
