#metadata_list=../data/metadata/bug_metadata/s_aureus.uk.txt
#in_dir=../data/genomes/hunt_et_al_v0.2/s_aureus
#out_dir=../data/genomes/hunt_et_al_v0.2/s_aureus.uk

bug_prefix=m_catarrhalis
metadata_list=../data/metadata/bug_metadata/$bug_prefix.uk.txt
in_dir=../data/genomes/hunt_et_al_v0.2/$bug_prefix
out_dir=../data/genomes/hunt_et_al_v0.2/$bug_prefix.uk


while read line;
do
    cp $in_dir/$line.fa $out_dir/$line.fna
done < $metadata_list



