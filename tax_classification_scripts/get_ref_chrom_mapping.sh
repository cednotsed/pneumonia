wkdir=/SAN/ugi/HAP_VAP/pneumonia
ref_dir=$wkdir/data/genomes/refseq_genomes_151223
out_path=$wkdir/data/metadata/refseq_genomes_151223/ref_to_chrom_mapping.csv
#/SAN/ugi/HAP_VAP/pneumonia/data/genomes/refseq_genomes_151223/Bacteria/ncbi_dataset/data

# Create file
echo $'file,chrom\n' > $out_path

find $ref_dir -type f| \
  grep -v refseq_genomes_151223.fna| \
  grep .fna| \
  xargs -I'{}' sh get_ref_chrom_mapping.worker.sh '{}' "$out_path"

