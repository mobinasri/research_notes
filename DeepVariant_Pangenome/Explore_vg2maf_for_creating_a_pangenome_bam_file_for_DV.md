### Comment 1: 01/19/2024
Personalized pangenomes has shown to be useful for acheiving higher accuracy in mapping reads to pangenome graphs, which can then helps variant calling for 
hard cases especially the ones suffering from read mismappings. I want to add the haplotype information from personalized pangenome to the DeepVariant examples, and
explore how this information can help DV for fixing the remaining errors. Since during creating personaliazed pangenome the haplotypes are selected based on how
similar they are to the sample's genome (kmer-based) this information has the potential to nagivate DeepVariant toward the correct call.
One feasible approach is converting the graph (in gbz format) to MAF and then to BAM which can then be passed to DV for example generation. How 
to represent these haplotypes in DV examples is a different subproject which will be explored somewhere else. 

Here I do initial tests for converting gbz to bam by using vg and vg2maf. vg2maf is a tool being developed by Glenn and it's available in this repo https://github.com/glennhickey/vg2maf/tree/development.
He said that I should use the devlopment branch for now. We need this conversion to be fast since this conversion should be repeated per sample (we have a 
separate personalized pangenome per sample) and its runtime will add up to the total runtime for the variant calling process. We can say that it's 
like a preprocessing step.

I wrote down a bash script that takes a gbz file and uses it to create vg, gam, dist files which are needed for creating a MAF file. This is the link to the bash script.
https://github.com/mobinasri/research_notes/blob/main/DeepVariant_Pangenome/files/gbz_to_maf.bash

```
# go to working directory
cd /private/groups/patenlab/masri/internship/convert_hprc_to_maf

# download vg2maf
wget http://public.gi.ucsc.edu/~hickey/debug/vg2maf-0369bf2b4b87ec4daed9848fd4e5bb9072d02cc7
mv vg2maf-0369bf2b4b87ec4daed9848fd4e5bb9072d02cc7 vg2maf_0369bf2b4b87ec4daed9848fd4e5bb9072d02cc7

# download vg
wget https://github.com/vgteam/vg/releases/download/v1.53.0/vg

# download gbz for HPRC v1.1 (https://github.com/human-pangenomics/hpp_pangenome_resources)
# this is not a personalized pangenome but the steps should be the same for conversion
mkdir hprc_graph_files && cd hprc_graph_files
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.gbz

cd ../

# run the script on Slurm
sbatch gbz_to_maf.bash $PWD/vg_v1.53.0 $PWD/vg2maf_0369bf2b4b87ec4daed9848fd4e5bb9072d02cc7 $PWD/hprc_graph_files/hprc-v1.1-mc-grch38.gbz

```


