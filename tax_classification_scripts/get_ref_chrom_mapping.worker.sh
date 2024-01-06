wkdir=/SAN/ugi/HAP_VAP/pneumonia
ref_dir=$wkdir/data/genomes/refseq_genomes_151223

#/SAN/ugi/HAP_VAP/pneumonia/data/genomes/refseq_genomes_151223/Bacteria/ncbi_dataset/data

file_path=$1
out_path=$2

accession=$(echo $file_path|cut -d'/' -f12)
cat $file_path| \
  grep \>| \
  sed "s|>||g"| \
  cut -d' ' -f1| \
  xargs -I'{}' echo $accession,'{}' >> $out_path

