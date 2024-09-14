# Run GATK on HG003 Novaseq 35x in three modes

This doc contains wdls for running GATK. The wdls were updated last time with GATKv4.6.0
https://broadinstitute.github.io/warp/docs/Pipelines/Whole_Genome_Germline_Single_Sample_Pipeline/README/

Here is the main wdl we have to run to obtain VCF from uBam.
https://github.com/broadinstitute/warp/blob/develop/pipelines/broad/dna_seq/germline/single_sample/wgs/WholeGenomeGermlineSingleSample.wdl

Based on documentaion we can run this workflow in three modes:
1. **Regaular GATK (without DRAGEN modules)**: Use default values
2. **GATK-DRAGEN (functional equivalence mode)**: Only by setting `dragen_functional_equivalence_mode = true` In this mode the output of pipeline will be functionally equivalent to the one produced with the DRAGEN hardware.
3. **GATK-DRAGEN (maximum quality mode)**: Only by setting `dragen_maximum_quality_mode = true` In this mode the pipeline use DRAGEN mapper and variant calling steps however the final results are **NOT** functionally equivalent to the one produced with the DRAGEN hardware.

### Note: I will set `use_gatk3_haplotype_caller = false` in all modes since I want to use GATK-4.6.

## Get Data
```
mkdir -p /private/groups/patenlab/masri/internship/external_callers/data
cd /private/groups/patenlab/masri/internship/external_callers/data

wget https://storage.googleapis.com/brain-genomics-public/research/sequencing/fastq/novaseq/wgs_pcr_free/35x/HG003.novaseq.pcr-free.35x.R1.fastq.gz
wget https://storage.googleapis.com/brain-genomics-public/research/sequencing/fastq/novaseq/wgs_pcr_free/35x/HG003.novaseq.pcr-free.35x.R2.fastq.gz
```

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
WORK_DIR="/private/groups/patenlab/masri/internship/external_callers/data/uBAM"
cd ${WORK_DIR}

cat paired-fastq-to-unmapped-bam.inputs.json 
{
  "ConvertPairedFastQsToUnmappedBamWf.readgroup_name": "HG003_Seq",
  "ConvertPairedFastQsToUnmappedBamWf.sample_name": "HG003",
  "ConvertPairedFastQsToUnmappedBamWf.fastq_1": "https://storage.googleapis.com/brain-genomics-public/research/sequencing/fastq/novaseq/wgs_pcr_free/35x/HG003.novaseq.pcr-free.35x.R1.fastq.gz",
  "ConvertPairedFastQsToUnmappedBamWf.fastq_2": "https://storage.googleapis.com/brain-genomics-public/research/sequencing/fastq/novaseq/wgs_pcr_free/35x/HG003.novaseq.pcr-free.35x.R2.fastq.gz", 
  "ConvertPairedFastQsToUnmappedBamWf.library_name": "Novaseq-HG003-35x",
  "ConvertPairedFastQsToUnmappedBamWf.platform_unit": "unit1",
  "ConvertPairedFastQsToUnmappedBamWf.run_date": "2024",
  "ConvertPairedFastQsToUnmappedBamWf.platform_name": "illumina",
  "ConvertPairedFastQsToUnmappedBamWf.sequencing_center": "XX",
  "ConvertPairedFastQsToUnmappedBamWf.make_fofn": true  
}

```

### Run Fastq to uBAM with Toil
```
cd /private/groups/patenlab/masri/internship/external_callers/data/uBAM

WDL_PATH="/private/groups/patenlab/masri/internship/external_callers/apps/seq-format-conversion/paired-fastq-to-unmapped-bam.wdl"
INPUT_JSON_PATH="/private/groups/patenlab/masri/internship/external_callers/data/uBAM/paired-fastq-to-unmapped-bam.inputs.json"
SAMPLE_NAME="HG003"
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
ls -lh /private/groups/patenlab/masri/internship/external_callers/data/uBAM/HG003/analysis/paired-fastq-to-unmapped-bam_outputs/HG003_Seq.unmapped.bam
-rwxr--r-- 1 masri patenlab 59G Sep 12 00:31 /private/groups/patenlab/masri/internship/external_callers/data/uBAM/HG003/analysis/paired-fastq-to-unmapped-bam_outputs/HG003_Seq.unmapped.bam
```

## Run GTAK variant caller in 3 modes

### Get and prepare WDL

```
mkdir -p /private/groups/patenlab/masri/internship/external_callers/apps
cd /private/groups/patenlab/masri/internship/external_callers/apps

# commit bc8c97d4db33c698d3256649ff817187389e78a4
git clone https://github.com/broadinstitute/warp
```

### Replace all `mv` with `ln -s`

Toil has some issues with `mv` if the file that is going to be moved is an input to the WDL task.
```
cd /private/groups/patenlab/masri/internship/external_callers/apps/warp

for wdl_file in $(find ./tasks/broad/ | grep ".wdl$");do sed -i 's|mv|ln \-s|g' ${wdl_file} ;done
for wdl_file in $(find ./pipelines/broad/ | grep ".wdl$");do sed -i 's|mv|ln \-s|g' ${wdl_file};done
```

### Remove `"` around cpu numbers because of an error in Toil
```
cd /private/groups/patenlab/masri/internship/external_callers/apps/warp

for wdl_file in $(find ./tasks/broad/ | grep ".wdl$");do sed -i '/cpu:/s/\"//g' ${wdl_file} ;done
for wdl_file in $(find ./pipelines/broad/ | grep ".wdl$");do sed -i '/cpu:/s/\"//g' ${wdl_file};done
```

### Fix indentation for the python codes written between `python3 <<CODE` and `CODE`
```
for wdl_file in $(find ./tasks/broad/ | grep ".wdl$");
do
  awk '
  /<<CODE/ {
      match($0, /^ */);            # Count leading spaces on the python line
      spaces = RLENGTH;            # Store the number of leading spaces
      in_block = 1;                # Set flag to indicate we are in the block between XXX-Y and ZZZ-Y
  }
  in_block && / CODE/ {
      sub("^ {" spaces "}", "");   # Remove the leading spaces from the CODE line
      print;                       # Print the CODE line
      in_block = 0;                # Reset the flag as we are out of the block
      next;
  }
  in_block {
      sub("^ {" spaces "}", "");   # Remove the leading spaces from lines within the block
  }
  {
      print;                       # Print all lines (modified or not)
  }' wdl_file
done

for wdl_file in $(find ./pipelines/broad/ | grep ".wdl$");
do
  awk '
  /<<CODE/ {
      match($0, /^ */);            # Count leading spaces on the python line
      spaces = RLENGTH;            # Store the number of leading spaces
      in_block = 1;                # Set flag to indicate we are in the block between XXX-Y and ZZZ-Y
  }
  in_block && / CODE/ {
      sub("^ {" spaces "}", "");   # Remove the leading spaces from the CODE line
      print;                       # Print the CODE line
      in_block = 0;                # Reset the flag as we are out of the block
      next;
  }
  in_block {
      sub("^ {" spaces "}", "");   # Remove the leading spaces from lines within the block
  }
  {
      print;                       # Print all lines (modified or not)
  }' wdl_file
done

```

### Make 3 input json files for 3 modes

```
cp /private/groups/patenlab/masri/internship/external_callers/apps/warp/pipelines/broad/dna_seq/germline/single_sample/wgs/test_inputs/Plumbing/dragen_mode_functional_equivalence.json .
sed 's|gs://|https://storage.googleapis.com/|g' dragen_mode_functional_equivalence.json 
```
I removed these two lines since they are specific for sample NA12878
```
  "WholeGenomeGermlineSingleSample.fingerprint_genotypes_file": "gs://broad-gotc-test-storage/single_sample/plumbing/bams/G96830.NA12878/G96830.NA12878.hg38.reference.fingerprint.vcf.gz",
  "WholeGenomeGermlineSingleSample.fingerprint_genotypes_index": "gs://broad-gotc-test-storage/single_sample/plumbing/bams/G96830.NA12878/G96830.NA12878.hg38.reference.fingerprint.vcf.gz.tbi",
```
The final jsons are uploaded here:

### run GATK in regular mode (no DRAGEN)

```
cd /private/groups/patenlab/masri/internship/external_callers/results/GATK/regular

WDL_PATH="/private/groups/patenlab/masri/internship/external_callers/apps/warp/pipelines/broad/dna_seq/germline/single_sample/wgs/WholeGenomeGermlineSingleSample.wdl"
INPUT_JSON_PATH="/private/groups/patenlab/masri/internship/external_callers/results/GATK/regular/gatk_regular.json"
SAMPLE_NAME="HG003_GATK_v4.6_regular_mode"
WDL_NAME=$(basename ${WDL_PATH%%.wdl})
mkdir -p ${WDL_NAME}_logs
EMAIL="masri@ucsc.edu"
USERNAME="masri"
RUN_WDL_BASH="/private/groups/patenlab/masri/internship/external_callers/apps/run_wdl_single_json.sh"

sbatch      --job-name=${WDL_NAME}_${USERNAME} \
            --cpus-per-task=64 \
            --mem=256G \
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
cd /private/groups/patenlab/masri/internship/external_callers/results/GATK/dragen_functional_equivalence

WDL_PATH="/private/groups/patenlab/masri/internship/external_callers/apps/warp/pipelines/broad/dna_seq/germline/single_sample/wgs/WholeGenomeGermlineSingleSample.wdl"
INPUT_JSON_PATH="/private/groups/patenlab/masri/internship/external_callers/results/GATK/dragen_functional_equivalence/dragen_mode_functional_equivalence.json"
SAMPLE_NAME="HG003_GATK_v4.6_dragen_functional_equivalence"
WDL_NAME=$(basename ${WDL_PATH%%.wdl})
mkdir -p ${WDL_NAME}_logs
EMAIL="masri@ucsc.edu"
USERNAME="masri"
RUN_WDL_BASH="/private/groups/patenlab/masri/internship/external_callers/apps/run_wdl_single_json.sh"

sbatch      --job-name=${WDL_NAME}_${USERNAME} \
            --cpus-per-task=64 \
            --mem=256G \
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
INPUT_JSON_PATH="/private/groups/patenlab/masri/internship/external_callers/results/GATK/dragen_maximum_quality/dragen_mode_maximum_quality.json"
SAMPLE_NAME="HG003_GATK_v4.6_dragen_maximum_quality"
WDL_NAME=$(basename ${WDL_PATH%%.wdl})
mkdir -p ${WDL_NAME}_logs
EMAIL="masri@ucsc.edu"
USERNAME="masri"
RUN_WDL_BASH="/private/groups/patenlab/masri/internship/external_callers/apps/run_wdl_single_json.sh"

sbatch      --job-name=${WDL_NAME}_${USERNAME} \
            --cpus-per-task=64 \
            --mem=256G \
            --mail-user=${EMAIL} \
            --output=${WDL_NAME}_logs/${WDL_NAME}_%A_%a.log \
            --partition=long  \
            --time=72:00:00 \
            ${RUN_WDL_BASH} \
            --wdl ${WDL_PATH} \
            --sample_name ${SAMPLE_NAME} \
            --input_json ${INPUT_JSON_PATH}
```
