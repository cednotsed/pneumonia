db_res=$wkdir/databases/resfinder_151223/resfinder_db
  run_resfinder.py \
    -ifa complete.fna \
    -o ./ \
    -s "Other" \
    -db_res $db_res \
    --nanopore \
    -acq \
    --min_cov 0.5 \
    --threshold 0.9


