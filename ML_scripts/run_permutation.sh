config_file=runs.csv
SLURM_ARRAY_TASK_ID=$1
experiment=$(head -n $SLURM_ARRAY_TASK_ID $config_file|tail -n 1)

g1=$(echo $experiment|cut -d',' -f1)
g2=$(echo $experiment|cut -d',' -f2)
col_name=$(echo $experiment|cut -d',' -f3)
threshold=$(echo $experiment|cut -d',' -f4)

python permutation.py $g1 $g2 $col_name $threshold
