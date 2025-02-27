bug_prefix=e_coli
in_dir=/mnt/c/git_repos/pneumonia/results/pangenome_out/$bug_prefix.filt/aligned_gene_sequences
n_files=$(ls $in_dir|wc -l)
echo $n_files

seq 1 $n_files|xargs -P2 -I{} Rscript mask_gene_alignments.R "$bug_prefix" {}
