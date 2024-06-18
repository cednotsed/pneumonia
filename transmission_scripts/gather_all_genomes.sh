#assembly_list=../data/metadata/bug_metadata/k_pneumoniae.assemblies.txt
#assembly_dir=../results/assembly_out/bin_out
#
#bg_list=../data/metadata/bug_metadata/k_pneumoniae_BG.mashtree.tree_trimmed_list_RTL_0.9
#bg_dir=../data/genomes/hunt_et_al_v0.2/k_pneumoniae.uk
#
#out_dir=../data/genomes/hunt_et_al_v0.2/k_pneumoniae.filt

bug_prefix=m_catarrhalis

assembly_list=../data/metadata/bug_metadata/$bug_prefix.assemblies.txt
assembly_dir=../results/assembly_out/bin_out

bg_list=../data/metadata/bug_metadata/${bug_prefix}_BG.mashtree.tree_trimmed_list_RTL_0.9
bg_dir=../data/genomes/hunt_et_al_v0.2/${bug_prefix}.uk

out_dir=../data/genomes/hunt_et_al_v0.2/${bug_prefix}.filt

mkdir $out_dir

# Get all background genome
while read line;
do
    cp $bg_dir/$line.fna $out_dir
done < $bg_list

# Get all assemblies
while read line;
do
    cp $assembly_dir/$line $out_dir
done < $assembly_list



