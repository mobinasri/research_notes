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

### Comment 2 : 01/19/2024

first I ran the script linked in [comment 1](https://github.com/mobinasri/research_notes/blob/main/DeepVariant_Pangenome/Explore_vg2maf_for_creating_a_pangenome_bam_file_for_DV.md#comment-1-01192024) with 128 Gb memory
```
--mem=128gb
```
It could create all vg, gam and dist files but crashed in the middle of creating MAF. Here is error log from slurm.
```
masri@phoenix:/private/groups/patenlab/masri/internship/convert_hprc_to_maf$ cat gbz_to_maf.1683797.log
/private/groups/patenlab/masri/internship/convert_hprc_to_maf
phoenix-07.prism
Wed Jan 17 07:55:28 PM PST 2024
[Wed Jan 17 07:55:28 PM PST 2024] Convert hprc-v1.1-mc-grch38.gbz hprc-v1.1-mc-grch38.vg
[Wed Jan 17 09:20:07 PM PST 2024] Convert hprc-v1.1-mc-grch38.gbz to hprc-v1.1-mc-grch38.gam
[Thu Jan 18 12:27:26 AM PST 2024] Sort and index hprc-v1.1-mc-grch38.gam
 break into sorted chunks       [=======================]100.0%
 merge 596 files                [=======================]100.0%
[Thu Jan 18 09:48:55 AM PST 2024] Create hprc-v1.1-mc-grch38.dist
[Thu Jan 18 01:31:15 PM PST 2024] Create hprc-v1.1-mc-grch38.maf
[vg2maf]: Loaded graph
/var/lib/slurm/slurmd/job1683797/slurm_script: line 60: 1728046 Killed                  ${VG2MAF_BIN} ${DIRNAME}/${PREFIX}.vg -d ${DIRNAME}/${PREFIX}.dist -r ${REF_SAMPLE} -g ${DIRNAME}/${PREFIX}.sort.gam -p --threads 8 > ${DIRNAME}/${PREFIX}.maf
slurmstepd-phoenix-07: error: Detected 2 oom-kill event(s) in StepId=1683797.batch. Some of your processes may have been killed by the cgroup out-of-memory handler.
```

Therefore I increased the memory size to 256 and reran the script after commenting all lines expect the last command. 

#### Runtimes
vg2maf is working fine for about 22 hours so far but not finished yet. Here is the most recently updated log from slurm.
```
cat gbz_to_maf.1685027.log
/private/groups/patenlab/masri/internship/convert_hprc_to_maf
phoenix-20.prism
Thu Jan 18 05:54:04 PM PST 2024
[Thu Jan 18 05:54:04 PM PST 2024] Create hprc-v1.1-mc-grch38.maf
[vg2maf]: Loaded graph
[vg2maf]: Applied position overlay
[vg2maf]: Loaded distance index
[vg2maf]: Loaded GAM index
[vg2maf]: Converting chain 0
```
This long runtime for vg2maf is a bit concerning to me. Other than vg2maf all previous steps for converting gbz to vg, gbz to gam and vg to dist are also time consuming. I think maybe I have to try a personalized pangenome with 4 or 8 haplotypes instead of the whole HPRC graph. This way I can have a better sense of how long it would take for a typical personlized pangenome.

| From | Convert To / Task    | Time (hrs:mins)    |
| :---:   | :---: | :---: |
| GBZ | VG   | 1:25   |
| GBZ | GAM   | 3:00   |
| GAM | sort and index   | 9:20  |
| VG | DIST   | 4:45   |
| Total |    | 18:30  |
