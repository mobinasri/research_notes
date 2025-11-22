# Running Giraffe and pangenome-aware DV with `graph-v2-sep8` on HG001, HG002 and HG003 samples (Element and Illumina)

## Remove CHM13 from `graph-v2-sep8`

I forked the vg_wdl repo and added a WDL for removing a sample from the graph `workflows/remove_sample_from_gbz.wdl`.
- vg_wdl repo: `https://github.com/mobinasri/vg_wdl`  
- commit: `6af010df5da1a0207b75309baa947e4cc84a44c8`

I ran these commands to generate an input json file using the csv files `data_table.csv` and `input_mapping.csv`:

(Instructions for using `launch_from_table.py` is available here https://github.com/human-pangenomics/hprc_intermediate_assembly/tree/main/hpc)
```
# set working directory
WORKING_DIR="/private/groups/patenlab/masri/haplotype_sampling/graph_v2.0/remove_chm13_vg_1.68"
cd ${WORKING_DIR}

## Save WDL path and name in environment variables
WDL_PATH="/private/groups/patenlab/masri/apps/vg_wdl/workflows/remove_sample_from_gbz.wdl"
WDL_FILENAME=$(basename ${WDL_PATH})
WDL_NAME=${WDL_FILENAME%%.wdl}


## Make a folder for saving files related to run e.g. input and output jsons
cd ${WORKING_DIR}
mkdir -p runs_toil_slurm
cd runs_toil_slurm

## Make a directory for saving input json files
mkdir -p ${WDL_NAME}_input_jsons
cd ${WDL_NAME}_input_jsons

python3 /private/groups/patenlab/masri/apps/hprc_intermediate_assembly/hpc/launch_from_table.py \
     --data_table ${WORKING_DIR}/data_table.csv \
     --field_mapping ${WORKING_DIR}/input_mapping.csv \
     --workflow_name ${WDL_NAME}
```
The content of `data_table.csv`
```
sample_id,gbz
hprc-sep8-mc-grch38-eval,/private/home/ghickey/dev/work/chrom-ordering/hprc-sep8-mc-grch38-eval/hprc-sep8-mc-grch38-eval.gbz
```

The content of `input_mapping.csv`
```
input,type,value
RemoveSampleFromGraph.IN_GBZ_FILE,scalar,$input.gbz
RemoveSampleFromGraph.CREATE_INDEX_OPTIONS,scalar,"--snarl-limit 1"
RemoveSampleFromGraph.SAMPLE_NAME_TO_REMOVE,scalar,"CHM13"
RemoveSampleFromGraph.CORES,scalar,16
RemoveSampleFromGraph.DOCKER_IMAGE,scalar,"quay.io/vgteam/vg:v1.68.0"
```

I ran these commands to sumbit the job to SLURM with Toil for removing CHM13 from `graph-v2-sep8` and creating a new graph.
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
mkdir -p runs_toil_slurm/${WDL_NAME}_logs
cd runs_toil_slurm

## Run jobs arrays
sbatch      --job-name=${WDL_NAME}_${USERNAME} \
            --cpus-per-task=32 \
            --mem=128G \
            --mail-user=${EMAIL} \
            --mail-type=FAIL,END \
            --output=${WDL_NAME}_logs/${WDL_NAME}_%A_%a.log \
            --array=[1]%1  \
            --time=${TIME_LIMIT} \
            --partition=${PARTITION} \
            /private/groups/hprc/qc_hmm_flagger/hprc_intermediate_assembly/hpc/toil_sbatch_single_machine.sh \
            --wdl ${WDL_PATH} \
            --sample_csv  ${WORKING_DIR}/data_table.csv \
            --input_json_path ${WORKING_DIR}/runs_toil_slurm/${WDL_NAME}_input_jsons/\${SAMPLE_ID}_${WDL_NAME}.json

```

The output gbz file is located here:
```
/private/groups/patenlab/masri/haplotype_sampling/graph_v2.0/remove_chm13_vg_1.68/runs_toil_slurm/hprc-sep8-mc-grch38-eval/analysis/remove_sample_from_gbz_outputs/hprc-sep8-mc-grch38-eval.CHM13_removed.gbz
```

## Haplotype sampling to 32 haplotypes

Next I haplotype-sampled the graphs

## Mapping reads with Giraffe

## Running pangenome-aware DV

## Evaluating with Happy
