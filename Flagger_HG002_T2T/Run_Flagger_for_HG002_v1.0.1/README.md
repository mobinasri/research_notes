### Comment 1 : 01/19/2024

#### Input data
Nancy from the T2T team asked me to run Flagger on the HG002 T2T assembly v1.0.1. She sent me these alignments to use for this analysis.

```
# HiFi Revio ~50x per haploid
https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/polishing/HG002/v1.0/mapping/hifi_revio_3cell/hg002v1.0_hifi_revio_3cells.pri.bam
https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/polishing/HG002/v1.0/mapping/hifi_revio_3cell/hg002v1.0_hifi_revio_3cells.pri.bam.bai

# ONT R10 UL ~60x per haploid
https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/polishing/HG002/v1.0/mapping/ont_r10_ul_dorado/hg002v1.0_ont_r10_ul_dorado.pri.bam
https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/polishing/HG002/v1.0/mapping/ont_r10_ul_dorado/hg002v1.0_ont_r10_ul_dorado.pri.bam.bai

# ONT R10 duplex ~40x per haploid
https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/polishing/HG002/v1.0/mapping/ont_r10_duplex/hg002v1.0_ont_r10_duplex.pri.bam
https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/polishing/HG002/v1.0/mapping/ont_r10_duplex/hg002v1.0_ont_r10_duplex.pri.bam.bai
```

### Comment 2 : 01/21/2024

#### Map to chm13v2.0
First I need to map each haplotype to CHM13 to extract coverage biased blocks.

```
# set env variables
USER_NAME="masri"
EMAIL="masri@ucsc.edu"
WDL_PATH="/private/groups/patenlab/masri/apps/flagger/wdls/tasks/alignment/asm2asm_aligner.wdl"
WDL_NAME=$(basename ${WDL_PATH%%.wdl})
INPUT_MAPPING_CSV="/private/groups/patenlab/masri/flagger/T2T/v1.0.1/map_to_chm13/input_mapping.csv"
INPUT_DATA_TABLE_CSV="/private/groups/patenlab/masri/flagger/T2T/v1.0.1/map_to_chm13/asm_data.csv"
WORKING_DIR="/private/groups/patenlab/masri/flagger/T2T/v1.0.1/map_to_chm13"

mkdir -p ${WORKING_DIR}
cd ${WORKING_DIR}

mkdir ${WDL_NAME}_output_jsons
mkdir -p ${WDL_NAME}_output_jsons
mkdir -p ${WDL_NAME}_input_jsons

LAUNCH_FROM_TABLE_PY="/private/groups/patenlab/masri/hprc/hprc_intermediate_assembly/hpc/launch_from_table.py"
LAUNCH_WORKFLOW_JOB_ARRAY_BASH="/private/groups/patenlab/masri/hprc/hprc_intermediate_assembly/hpc/launch_workflow_job_array.sh"

# create input json files
cd ${WORKING_DIR}/${WDL_NAME}_input_jsons
python3  ${LAUNCH_FROM_TABLE_PY} \
    --data_table ${INPUT_DATA_TABLE_CSV} \
    --field_mapping ${INPUT_MAPPING_CSV} \
    --workflow_name ${WDL_NAME}
```

Take a look at a input json file:
```
cat HG002_v1.0.1_mat_asm2asm_aligner.json
{
  "asm2asmAlignment.queryAssemblyFastaGz": "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/hg002v1.0.1.mat.fasta.gz",
  "asm2asmAlignment.alignmentBam.threadCount": 32,
  "asm2asmAlignment.refAssemblyFastaGz": "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz",
  "asm2asmAlignment.alignmentBam.memSize": 64,
  "asm2asmAlignment.alignmentBam.diskSize": 128,
  "asm2asmAlignment.preset": "asm5",
  "asm2asmAlignment.suffix": "to_chm13v2.0",
  "asm2asmAlignment.alignmentBam.options": "--eqx --cs -L",
  "asm2asmAlignment.aligner": "minimap2"
}
```

Run the alignments
```
# run the job array
cd ${WORKING_DIR}
cd ${WDL_NAME}_output_jsons
TOIL_SLURM_ARGS="--time=3-0:00 --partition=long"
INPUT_JSON_DIR=${WORKING_DIR}/${WDL_NAME}_input_jsons
mkdir ${WDL_NAME}_logs

sbatch      --job-name=${WDL_NAME}_${USERNAME} \
            --mail-user=${EMAIL} \
            --output=${WDL_NAME}_logs/${WDL_NAME}_%x_%j_%A_%a.log \
            --array=1-2%2  \
            --partition=long  \
            ${LAUNCH_WORKFLOW_JOB_ARRAY_BASH} \
            ${INPUT_DATA_TABLE_CSV} \
            ${WDL_PATH} \
            "${TOIL_SLURM_ARGS}"  \
            ${INPUT_JSON_DIR}

```

I'll update once the alignments are finished.

