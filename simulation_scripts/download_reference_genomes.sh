acc_list=/mnt/c/git_repos/pneumonia/data/metadata/bugs_to_download.curated.accessions_only.txt
out_dir=/mnt/c/git_repos/pneumonia/data/genomes/bug_references

ncbi-genome-download \
    --formats fasta \
    --flat-output \
    --progress-bar \
    --assembly-accessions $acc_list \
    --output-folder $out_dir \
    all
