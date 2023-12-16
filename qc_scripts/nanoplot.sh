sum_path=../data/sequencing_summaries/INHALE_FRESH_5_sequencing_summary.txt

NanoPlot \
  -t 6 \
  --summary $sum_path \
  --loglength \
  --outdir ./test
#  --maxlength 100000 \
#  --plots {kde,hex,dot} \
#  --outdir ./test

