wkdir=/mnt/c/git_repos/pneumonia
ref_list=$wkdir/data/metadata/bugs.nanosim.meta.tsv
abun_list=$wkdir/data/metadata/bugs.nanosim.abundance.tsv
out_prefix=$wkdir/data/simulated_fastqs/test
threads=8

model_prefix=$wkdir/databases/nanosim_db/metagenome_ERR3152364_Even/training

~/NanoSim/src/simulator.py metagenome \
    --genome_list $ref_list \
    --abun $abun_list \
    --model_prefix $model_prefix \
    --median_len 1500 \
    --sd_len 100 \
    --seed 66 \
    --basecaller guppy \
    --fastq \
    --num_threads $threads

