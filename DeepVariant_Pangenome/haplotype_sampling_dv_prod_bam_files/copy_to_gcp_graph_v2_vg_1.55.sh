#!/bin/bash
#SBATCH --job-name=gcp_copy
#SBATCH --partition=long
#SBATCH --mail-user=masri@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=8gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --output=%x.%j.log
#SBATCH --time=72:00:00

pwd; hostname; date

set -e


DIR="/private/groups/patenlab/masri/haplotype_sampling/dv_prod_bam_files/non_diploid/runs_toil_slurm_graph_v2_vg_1.55"
for i in $(find ${DIR} | grep "\.gbz$");do echo $i; gsutil cp $i gs://pepper-deepvariant/mobinasri/haplotype_sampling/dv_prod_bam_files/non_diploid_sampling/graph_v2_vg_1.55/ ;done
