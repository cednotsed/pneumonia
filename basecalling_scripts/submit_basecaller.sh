base_dir=/mnt/c/git_repos/pneumonia/data/raw_fast5s
for i in $base_dir/*
do
  echo $i
  sh basecall_single_run.sh "$i"
done
