### Comment 1 : 04/11/2025 

### Remove chm13 from graph-v2.0 
### Creating haplotype-sampled HPRC-v2.0 graphs for 48 training production bam files 

```
cat input_mapping_graph_v2_vg_1.55.csv

input,type,value
HaplotypeSampling.IN_OUTPUT_NAME_PREFIX,scalar,$input.sample_id
HaplotypeSampling.DIPLOID,scalar,"false"
HaplotypeSampling.INPUT_READ_FILE_ARRAY,array,$input.alignment_bam_array
HaplotypeSampling.HAPLOTYPE_NUMBER_ARRAY,array,"[8, 16, 32, 64]"
HaplotypeSampling.REFERENCE_FASTA,scalar,"/private/groups/patenlab/masri/common/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"
HaplotypeSampling.HAPL_FILE,scalar,"/private/groups/patenlab/masri/haplotype_sampling/graph_v2.0/remove_chm13/runs_toil_slurm/hprc-v2.0-mc-grch38-eval/analysis/remove_sample_from_gbz_outputs/hprc-v2.0-mc-grch38-eval.CHM13_removed.hapl"
HaplotypeSampling.IN_DIST_FILE,scalar,"/private/groups/patenlab/masri/haplotype_sampling/graph_v2.0/remove_chm13/runs_toil_slurm/hprc-v2.0-mc-grch38-eval/analysis/remove_sample_from_gbz_outputs/hprc-v2.0-mc-grch38-eval.CHM13_removed.dist"
HaplotypeSampling.IN_GBZ_FILE,scalar,"/private/groups/patenlab/masri/haplotype_sampling/graph_v2.0/remove_chm13/runs_toil_slurm/hprc-v2.0-mc-grch38-eval/analysis/remove_sample_from_gbz_outputs/hprc-v2.0-mc-grch38-eval.CHM13_removed.gbz"
HaplotypeSampling.CORES,scalar,8
HaplotypeSampling.DOCKER_IMAGE,scalar,"quay.io/vgteam/vg:v1.55.0"
```


#### Create json files
```
# set working directory
WORKING_DIR="/private/groups/patenlab/masri/haplotype_sampling/dv_prod_bam_files/non_diploid"
cd ${WORKING_DIR}

## Save WDL path and name in environment variables
WDL_PATH="/private/groups/patenlab/masri/apps/vg_wdl/workflows/haplotype_sampling_customized.wdl"
WDL_FILENAME=$(basename ${WDL_PATH})
WDL_NAME=${WDL_FILENAME%%.wdl}


## Make a folder for saving files related to run e.g. input and output jsons
cd ${WORKING_DIR}
mkdir -p runs_toil_slurm_graph_v2_vg_1.55
cd runs_toil_slurm_graph_v2_vg_1.55

## Make a directory for saving input json files
mkdir -p ${WDL_NAME}_input_jsons
cd ${WDL_NAME}_input_jsons

python3 /private/groups/patenlab/masri/apps/hprc_intermediate_assembly/hpc/launch_from_table.py \
     --data_table ${WORKING_DIR}/data_table.csv \
     --field_mapping ${WORKING_DIR}/input_mapping_graph_v2_vg_1.55.csv \
     --workflow_name ${WDL_NAME}

```

#### Run workflow
```
## Make sure you are in the working directory
cd ${WORKING_DIR}

## Set environment variables for sbatch
USERNAME="masri"
EMAIL="masri@ucsc.edu"
TIME_LIMIT="48:00:00"

## Partition should be modifed based on the available partitions on the server
PARTITION="long"


## Go to the execution directory
mkdir -p runs_toil_slurm_graph_v2_vg_1.55/${WDL_NAME}_logs
cd runs_toil_slurm_graph_v2_vg_1.55

## Run jobs arrays
sbatch      --job-name=${WDL_NAME}_${USERNAME} \
            --cpus-per-task=48 \
            --mem=256G \
            --mail-user=${EMAIL} \
            --mail-type=FAIL,END \
            --output=${WDL_NAME}_logs/${WDL_NAME}_%A_%a.log \
            --array=[1-48]%10  \
            --time=${TIME_LIMIT} \
            --partition=${PARTITION} \
            /private/groups/hprc/qc_hmm_flagger/hprc_intermediate_assembly/hpc/toil_sbatch_single_machine.sh \
            --wdl ${WDL_PATH} \
            --sample_csv  ${WORKING_DIR}/data_table.csv \
            --input_json_path ${WORKING_DIR}/runs_toil_slurm_graph_v2_vg_1.55/${WDL_NAME}_input_jsons/\${SAMPLE_ID}_${WDL_NAME}.json

```

#### Copy to gcp bucket
```
sbatch copy_to_gcp_graph_v2_vg_1.55.sh
```
