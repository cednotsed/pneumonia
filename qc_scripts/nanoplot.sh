fastq_path=../data/basecalled_fastqs/SaminaPool1-barcode12.fastq.gz

NanoPlot \
  -t 2 \
  --fastq $fastq_path \
  --maxlength 100000 \
  --plots {kde,hex,dot} \
  --outdir ./test

