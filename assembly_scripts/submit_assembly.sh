#for i in /SAN/ugi/HAP_VAP/pneumonia/data/basecalled_fastqs/*.fastq.gz
#do
#  echo $i
#  sh assemble_reads.sh $i
#done

ls /SAN/ugi/HAP_VAP/pneumonia/data/basecalled_fastqs/*.fastq.gz|xargs -P4 -I'{}' sh assemble_reads.sh '{}'
