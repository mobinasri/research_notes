# Run GATK on HG003 Element 1000bp_ins in four modes (3 modes of GATK-v4.6 + 1 mode of GATK-v3)
### (Repeating https://github.com/mobinasri/research_notes/edit/main/DeepVariant_Pangenome/External_Variant_Callers/GATK/README.md for Element data)

This doc contains wdls for running GATK. The wdls were updated last time with GATKv4.6.0
https://broadinstitute.github.io/warp/docs/Pipelines/Whole_Genome_Germline_Single_Sample_Pipeline/README/

Here is the main wdl we have to run to obtain VCF from uBam.
https://github.com/broadinstitute/warp/blob/develop/pipelines/broad/dna_seq/germline/single_sample/wgs/WholeGenomeGermlineSingleSample.wdl

Based on documentaion we can run this workflow in three modes:
1. **Regaular GATK (without DRAGEN modules)**: Use default values
2. **GATK-DRAGEN (functional equivalence mode)**: Only by setting `dragen_functional_equivalence_mode = true` In this mode the output of pipeline will be functionally equivalent to the one produced with the DRAGEN hardware.
3. **GATK-DRAGEN (maximum quality mode)**: Only by setting `dragen_maximum_quality_mode = true` In this mode the pipeline use DRAGEN mapper and variant calling steps however the final results are **NOT** functionally equivalent to the one produced with the DRAGEN hardware.

### Note: I will set `use_gatk3_haplotype_caller = false` in all modes since I want to use GATK-4.6.

## Convert fastq to uBam (since GATK workflow takes only uBAM as input)
### Get WDL
```
APPS_DIR="/private/groups/patenlab/masri/internship/external_callers/apps"
mkdir -p  ${APPS_DIR}
cd ${APPS_DIR}
git clone https://github.com/gatk-workflows/seq-format-conversion
```

### Prepare json
```
WORK_DIR="/private/groups/patenlab/masri/internship/external_callers/data/element/uBAM"
cd ${WORK_DIR}

cat paired-fastq-to-unmapped-bam.inputs.json 
{
  "ConvertPairedFastQsToUnmappedBamWf.readgroup_name": "HG003_Element",
  "ConvertPairedFastQsToUnmappedBamWf.sample_name": "HG003",
  "ConvertPairedFastQsToUnmappedBamWf.fastq_1": "https://storage.googleapis.com/brain-genomics-public/research/element/cloudbreak_wgs/HG003.element.cloudbreak.1000bp_ins.R1.fastq.gz",
  "ConvertPairedFastQsToUnmappedBamWf.fastq_2": "https://storage.googleapis.com/brain-genomics-public/research/element/cloudbreak_wgs/HG003.element.cloudbreak.1000bp_ins.R2.fastq.gz", 
  "ConvertPairedFastQsToUnmappedBamWf.library_name": "HG003_Element_1000bp_ins",
  "ConvertPairedFastQsToUnmappedBamWf.platform_unit": "unit1",
  "ConvertPairedFastQsToUnmappedBamWf.run_date": "2024",
  "ConvertPairedFastQsToUnmappedBamWf.platform_name": "element",
  "ConvertPairedFastQsToUnmappedBamWf.sequencing_center": "XX",
  "ConvertPairedFastQsToUnmappedBamWf.make_fofn": true  
}

```

### Run Fastq to uBAM with Toil
```
cd /private/groups/patenlab/masri/internship/external_callers/data/element/uBAM

WDL_PATH="/private/groups/patenlab/masri/internship/external_callers/apps/seq-format-conversion/paired-fastq-to-unmapped-bam.wdl"
INPUT_JSON_PATH="/private/groups/patenlab/masri/internship/external_callers/data/element/uBAM/paired-fastq-to-unmapped-bam.inputs.json"
SAMPLE_NAME="HG003-Element"
WDL_NAME=$(basename ${WDL_PATH%%.wdl})
mkdir -p ${WDL_NAME}_logs
EMAIL="masri@ucsc.edu"
USERNAME="masri"
RUN_WDL_BASH="/private/groups/patenlab/masri/internship/external_callers/apps/run_wdl_single_json.sh"

sbatch      --job-name=${WDL_NAME}_${USERNAME} \
            --cpus-per-task=32 \
            --mem=32G \
            --mail-user=${EMAIL} \
            --output=${WDL_NAME}_logs/${WDL_NAME}_%A_%a.log \
            --partition=long  \
            --time=6:00:00 \
            ${RUN_WDL_BASH} \
            --wdl ${WDL_PATH} \
            --sample_name ${SAMPLE_NAME} \
            --input_json ${INPUT_JSON_PATH}
```

### Output uBAM
```
ls -lh /private/groups/patenlab/masri/internship/external_callers/data/element/uBAM/HG003-Element/analysis/paired-fastq-to-unmapped-bam_outputs/
total 86G
-rwxr--r-- 1 masri patenlab 86G Dec 23 16:20 HG003_Element.unmapped.bam
-rwxr--r-- 1 masri patenlab  78 Dec 23 16:16 HG003.ubam.list
```


### run GATK in regular mode (no DRAGEN and `use_gatk3_haplotype_caller = true`)

```
cd /private/groups/patenlab/masri/internship/external_callers/results/GATK_Element/regular_v3

WDL_PATH="/private/groups/patenlab/masri/internship/external_callers/apps/warp/pipelines/broad/dna_seq/germline/single_sample/wgs/WholeGenomeGermlineSingleSample.wdl"
INPUT_JSON_PATH="/private/groups/patenlab/masri/internship/external_callers/results/GATK_Element/regular_v3/gatk_v3_regular.json"
SAMPLE_NAME="HG003_Element_GATK_v3_regular_mode"
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
cd /private/groups/patenlab/masri/internship/external_callers/results/GATK_Element/regular_v4.6

WDL_PATH="/private/groups/patenlab/masri/internship/external_callers/apps/warp/pipelines/broad/dna_seq/germline/single_sample/wgs/WholeGenomeGermlineSingleSample.wdl"
INPUT_JSON_PATH="/private/groups/patenlab/masri/internship/external_callers/results/GATK_Element/regular_v4.6/gatk_v4.6_regular.json"
SAMPLE_NAME="HG003_Element_GATK_v4.6_regular_mode"
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
cd /private/groups/patenlab/masri/internship/external_callers/results/GATK_Element/dragen_functional_equivalence

WDL_PATH="/private/groups/patenlab/masri/internship/external_callers/apps/warp/pipelines/broad/dna_seq/germline/single_sample/wgs/WholeGenomeGermlineSingleSample.wdl"
INPUT_JSON_PATH="/private/groups/patenlab/masri/internship/external_callers/results/GATK_Element/dragen_functional_equivalence/dragen_mode_functional_equivalence.json"
SAMPLE_NAME="HG003_Element_GATK_v4.6_dragen_functional_equivalence"
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
cd /private/groups/patenlab/masri/internship/external_callers/results/GATK/dragen_maximum_quality

WDL_PATH="/private/groups/patenlab/masri/internship/external_callers/apps/warp/pipelines/broad/dna_seq/germline/single_sample/wgs/WholeGenomeGermlineSingleSample.wdl"
INPUT_JSON_PATH="/private/groups/patenlab/masri/internship/external_callers/results/GATK_Element/dragen_maximum_quality/dragen_mode_maximum_quality.json"
SAMPLE_NAME="HG003_Element_GATK_v4.6_dragen_maximum_quality"
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

## Copying vcf files to gs bucket
```

```
