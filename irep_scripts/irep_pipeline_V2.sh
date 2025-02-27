wkdir=/mnt/c/git_repos/pneumonia
idx=$1
config_file=$wkdir/data/metadata/irep/irep.config.tsv
out_dir=$wkdir/results/irep_out/mapping_out
irep_dir=$wkdir/results/irep_out/irep_results

mkdir $out_dir
mkdir $irep_dir

row=$(head -n $idx $config_file|tail -n 1)

prefix=$(echo "$row"|awk -F'\t' '{print $1}')
taxid=$(echo "$row"|awk -F'\t' '{print $3}')
ref_path=$(echo "$row"|awk -F'\t' '{print $4}')
fastq_path=$(echo "$row"|awk -F'\t' '{print $5}')

echo $prefix
echo $taxid
echo $ref_path

minimap2 \
    -ax map-ont \
    -t 4 \
    $ref_path $fastq_path | \
    samtools view -h -F4 -F2048 - |samtools sort - -o $out_dir/$prefix.$taxid.sam

bPTR \
    -f $ref_path \
    -s $out_dir/$prefix.$taxid.sam \
    -o $irep_dir/$prefix.$taxid.bPTR.tsv \
    -plot $irep_dir/$prefix.$taxid.pdf \
    -m coverage \
    -ff
