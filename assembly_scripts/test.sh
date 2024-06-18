for i in medaka*;
do
    if grep -q "oom" $i
    then
        cat $i|head -n 1|sed "s|/flask/scratch/matthewsp/pneumonia/data/basecalled_fastqs/no_humans/||g"
    fi
done
