#for i in /SAN/ugi/HAP_VAP/pneumonia/data/basecalled_fastqs/*.fastq.gz
#do
#  #echo $i
#  sh polish_assembly.sh $i
#done

ls /SAN/ugi/HAP_VAP/pneumonia/data/basecalled_fastqs/*.fastq.gz|xargs -P4 -I'{}' sh polish_assembly.sh '{}'
