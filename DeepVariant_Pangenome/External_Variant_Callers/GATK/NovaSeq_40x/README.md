# Running GATK on NovaSeq 40x

### run GATK in regular mode (no DRAGEN and `use_gatk3_haplotype_caller = true`)
```
cd /private/groups/patenlab/masri/internship/external_callers/results/HG002/GATK_40x/regular_v3

WDL_PATH="/private/groups/patenlab/masri/internship/external_callers/apps/warp/pipelines/broad/dna_seq/germline/single_sample/wgs/WholeGenomeGermlineSingleSample.wdl"
INPUT_JSON_PATH="/private/groups/patenlab/masri/internship/external_callers/results/HG002/GATK_40x/regular_v3/gatk_v3_regular.json"
SAMPLE_NAME="HG002_40x_GATK_v3_regular_mode"
WDL_NAME=$(basename ${WDL_PATH%%.wdl})
mkdir -p ${WDL_NAME}_logs
EMAIL="masri@ucsc.edu"
USERNAME="masri"
RUN_WDL_BASH="/private/groups/patenlab/masri/internship/external_callers/apps/run_wdl_single_json.sh"

sbatch      --job-name=${WDL_NAME}_${USERNAME} \
            --cpus-per-task=64 \
            --mem=250G \
            --mail-user=${EMAIL} \
            --output=${WDL_NAME}_logs/${WDL_NAME}_%A_%a.log \
            --partition=long  \
            --time=72:00:00 \
            ${RUN_WDL_BASH} \
            --wdl ${WDL_PATH} \
            --sample_name ${SAMPLE_NAME} \
            --input_json ${INPUT_JSON_PATH}

```

### run GATK in regular mode (no DRAGEN and `use_gatk3_haplotype_caller = false`)
```
cd /private/groups/patenlab/masri/internship/external_callers/results/HG002/GATK_40x/regular_v4.6

WDL_PATH="/private/groups/patenlab/masri/internship/external_callers/apps/warp/pipelines/broad/dna_seq/germline/single_sample/wgs/WholeGenomeGermlineSingleSample.wdl"
INPUT_JSON_PATH="/private/groups/patenlab/masri/internship/external_callers/results/HG002/GATK_40x/regular_v4.6/gatk_v4.6_regular.json"
SAMPLE_NAME="HG002_40x_GATK_v4.6_regular_mode"
WDL_NAME=$(basename ${WDL_PATH%%.wdl})
mkdir -p ${WDL_NAME}_logs
EMAIL="masri@ucsc.edu"
USERNAME="masri"
RUN_WDL_BASH="/private/groups/patenlab/masri/internship/external_callers/apps/run_wdl_single_json.sh"

sbatch      --job-name=${WDL_NAME}_${USERNAME} \
            --cpus-per-task=64 \
            --mem=250G \
            --mail-user=${EMAIL} \
            --output=${WDL_NAME}_logs/${WDL_NAME}_%A_%a.log \
            --partition=long  \
            --time=72:00:00 \
            ${RUN_WDL_BASH} \
            --wdl ${WDL_PATH} \
            --sample_name ${SAMPLE_NAME} \
            --input_json ${INPUT_JSON_PATH}
```

### run GATK in `dragen_functional_equivalence` mode

```
cd /private/groups/patenlab/masri/internship/external_callers/results/HG002/GATK_40x/dragen_functional_equivalence

WDL_PATH="/private/groups/patenlab/masri/internship/external_callers/apps/warp/pipelines/broad/dna_seq/germline/single_sample/wgs/WholeGenomeGermlineSingleSample.wdl"
INPUT_JSON_PATH="/private/groups/patenlab/masri/internship/external_callers/results/HG002/GATK_40x/dragen_functional_equivalence/dragen_mode_functional_equivalence.json"
SAMPLE_NAME="HG002_40x_GATK_v4.6_dragen_functional_equivalence"
WDL_NAME=$(basename ${WDL_PATH%%.wdl})
mkdir -p ${WDL_NAME}_logs
EMAIL="masri@ucsc.edu"
USERNAME="masri"
RUN_WDL_BASH="/private/groups/patenlab/masri/internship/external_callers/apps/run_wdl_single_json.sh"

sbatch      --job-name=${WDL_NAME}_${USERNAME} \
            --cpus-per-task=64 \
            --mem=250G \
            --mail-user=${EMAIL} \
            --output=${WDL_NAME}_logs/${WDL_NAME}_%A_%a.log \
            --partition=long  \
            --time=72:00:00 \
            ${RUN_WDL_BASH} \
            --wdl ${WDL_PATH} \
            --sample_name ${SAMPLE_NAME} \
            --input_json ${INPUT_JSON_PATH}
```

### run GATK in `dragen_maximum_quality` mode
```
cd /private/groups/patenlab/masri/internship/external_callers/results/HG002/GATK_40x/dragen_maximum_quality

WDL_PATH="/private/groups/patenlab/masri/internship/external_callers/apps/warp/pipelines/broad/dna_seq/germline/single_sample/wgs/WholeGenomeGermlineSingleSample.wdl"
INPUT_JSON_PATH="/private/groups/patenlab/masri/internship/external_callers/results/HG002/GATK_40x/dragen_maximum_quality/dragen_mode_maximum_quality.json"
SAMPLE_NAME="HG002_40x_GATK_v4.6_dragen_maximum_quality"
WDL_NAME=$(basename ${WDL_PATH%%.wdl})
mkdir -p ${WDL_NAME}_logs
EMAIL="masri@ucsc.edu"
USERNAME="masri"
RUN_WDL_BASH="/private/groups/patenlab/masri/internship/external_callers/apps/run_wdl_single_json.sh"

sbatch      --job-name=${WDL_NAME}_${USERNAME} \
            --cpus-per-task=64 \
            --mem=250G \
            --mail-user=${EMAIL} \
            --output=${WDL_NAME}_logs/${WDL_NAME}_%A_%a.log \
            --partition=long  \
            --time=72:00:00 \
            ${RUN_WDL_BASH} \
            --wdl ${WDL_PATH} \
            --sample_name ${SAMPLE_NAME} \
            --input_json ${INPUT_JSON_PATH}
```
