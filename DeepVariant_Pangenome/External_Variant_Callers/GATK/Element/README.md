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

### Platform should be uppercase in ubam header
This is due to this bug I received in `regular_v3` and `regular_v4.6` modes

```
ERROR::INVALID_PLATFORM_VALUE:Read name HG003_Element, The platform (PL) attribute (element) + was not one of the valid values for read group
```
It seems that I should change `element` to `ELEMENT` to have a valid PL based on sam specification.
```
# get header and modify PL
cd /private/groups/patenlab/masri/internship/external_callers/data/element/uBAM/HG003-Element/analysis/paired-fastq-to-unmapped-bam_outputs
samtools view -H HG003_Element.unmapped.bam > HG003_Element.unmapped.header.sam
sed 's|PL:element|PL:ELEMENT|g' HG003_Element.unmapped.header.sam > HG003_Element.unmapped.header.modified.sam

# get docker and run interactively
docker run -it --rm -v/private/groups/patenlab/masri:/private/groups/patenlab/masri -u $(id -u):$(id -g) us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10

cd /private/groups/patenlab/masri/internship/external_callers/data/element/uBAM/HG003-Element/analysis/paired-fastq-to-unmapped-bam_outputs

# make a new bam with new header
java -jar /usr/picard/picard.jar ReplaceSamHeader \
    I=HG003_Element.unmapped.bam \
    HEADER=HG003_Element.unmapped.header.modified.sam \
    O=HG003_Element.unmapped.PL_ELEMENT.bam
```

### run GATK in regular mode (no DRAGEN and `use_gatk3_haplotype_caller = true`)

```
cd /private/groups/patenlab/masri/internship/external_callers/results/GATK_Element/regular_v3/rerun

WDL_PATH="/private/groups/patenlab/masri/internship/external_callers/apps/warp/pipelines/broad/dna_seq/germline/single_sample/wgs/WholeGenomeGermlineSingleSample.wdl"
INPUT_JSON_PATH="/private/groups/patenlab/masri/internship/external_callers/results/GATK_Element/regular_v3/rerun/gatk_v3_regular.json"
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
```
cd /private/groups/patenlab/masri/internship/external_callers/results/GATK_Element/regular_v3/rerun/HG003_Element_GATK_v3_regular_mode/analysis/WholeGenomeGermlineSingleSample_outputs

bcftools view -Oz \
    -f 'PASS,.' \
    HG003_Element_1000bp_ins_GATK_Regular_v3.rb.g.vcf.gz \
    -o HG003_Element_1000bp_ins_GATK_Regular_v3.rb.g.PASS.vcf.gz

tabix  HG003_Element_1000bp_ins_GATK_Regular_v3.rb.g.PASS.vcf.gz
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

```
cd /private/groups/patenlab/masri/internship/external_callers/results/GATK_Element/regular_v4.6/HG003_Element_GATK_v4.6_regular_mode/analysis/WholeGenomeGermlineSingleSample_outputs

bcftools view -Oz \
    -f 'PASS,.' \
    HG003_Element_1000bp_ins_GATK_Regular_v4.6.rb.g.vcf.gz \
    -o HG003_Element_1000bp_ins_GATK_Regular_v4.6.rb.g.PASS.vcf.gz

tabix  HG003_Element_1000bp_ins_GATK_Regular_v4.6.rb.g.PASS.vcf.gz
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
cd /private/groups/patenlab/masri/internship/external_callers/results/GATK_Element/dragen_maximum_quality

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

## Bug while running `GATK_v4.6_dragen_functional_equivalence`:
```
	[139869693069056]	ERROR: This thread caught an exception first
	When maskLen < 15, the function ssw_align doesn't return 2nd best alignment information.
	When maskLen < 15, the function ssw_align doesn't return 2nd best alignment information.
	When maskLen < 15, the function ssw_align doesn't return 2nd best alignment information.
	When maskLen < 15, the function ssw_align doesn't return 2nd best alignment information.
	When maskLen < 15, the function ssw_align doesn't return 2nd best alignment information.
	When maskLen < 15, the function ssw_align doesn't return 2nd best alignment information.
	When maskLen < 15, the function ssw_align doesn't return 2nd best alignment information.
	When maskLen < 15, the function ssw_align doesn't return 2nd best alignment information.
	When maskLen < 15, the function ssw_align doesn't return 2nd best alignment information.
	When maskLen < 15, the function ssw_align doesn't return 2nd best alignment information.
	When maskLen < 15, the function ssw_align doesn't return 2nd best alignment information.
	When maskLen < 15, the function ssw_align doesn't return 2nd best alignment information.
	When maskLen < 15, the function ssw_align doesn't return 2nd best alignment information.
	When maskLen < 15, the function ssw_align doesn't return 2nd best alignment information.
	When maskLen < 15, the function ssw_align doesn't return 2nd best alignment information.
	When maskLen < 15, the function ssw_align doesn't return 2nd best alignment information.
	When maskLen < 15, the function ssw_align doesn't return 2nd best alignment information.
	When maskLen < 15, the function ssw_align doesn't return 2nd best alignment information.
	Error: : Invalid argument: /usr/gitc/temp/DRAGMAP-1.2.1/src/lib/sequences/Seed.cpp(64): Throw in function dragenos::sequences::Seed::Data dragenos::sequences::Seed::getPrimaryData(bool) const
	Dynamic exception type: boost::wrapexcept<dragenos::common::PreConditionException>
	std::exception::what: Requesting primary data for an invalid seed
	: Requesting primary data for an invalid seed

```
It seems that this problem is mentioned and fixed in [DRAGMAP github](https://github.com/Illumina/DRAGMAP/issues/23) repo but it is not yet updated in [warp repo](https://github.com/broadinstitute/warp/blob/2fffa8750e4c42832aaafe65df63d06a44c31f96/tasks/broad/DragmapAlignment.wdl#L34) (the docker should be updated with v1.3 of DRAGMAP).

## Copying vcf files to gs bucket
```
gsutil ls -R gs://pepper-deepvariant/mobinasri/pangenome_paper/external_callers_element | grep ".vcf.gz"

gs://pepper-deepvariant/mobinasri/pangenome_paper/external_callers_element/GATK/regular_v3/HG003_GATK_v3_regular_mode/analysis/WholeGenomeGermlineSingleSample_outputs/HG003_Element_1000bp_ins_GATK_Regular_v3.rb.g.PASS.vcf.gz
gs://pepper-deepvariant/mobinasri/pangenome_paper/external_callers_element/GATK/regular_v3/HG003_GATK_v3_regular_mode/analysis/WholeGenomeGermlineSingleSample_outputs/HG003_Element_1000bp_ins_GATK_Regular_v3.rb.g.PASS.vcf.gz.tbi
gs://pepper-deepvariant/mobinasri/pangenome_paper/external_callers_element/GATK/regular_v3/HG003_GATK_v3_regular_mode/analysis/WholeGenomeGermlineSingleSample_outputs/HG003_Element_1000bp_ins_GATK_Regular_v3.rb.g.vcf.gz
gs://pepper-deepvariant/mobinasri/pangenome_paper/external_callers_element/GATK/regular_v3/HG003_GATK_v3_regular_mode/analysis/WholeGenomeGermlineSingleSample_outputs/HG003_Element_1000bp_ins_GATK_Regular_v3.rb.g.vcf.gz.tbi
gs://pepper-deepvariant/mobinasri/pangenome_paper/external_callers_element/GATK/regular_v4.6/HG003_GATK_v4.6_regular_mode/analysis/WholeGenomeGermlineSingleSample_outputs/HG003_Element_1000bp_ins_GATK_Regular_v4.6.rb.g.PASS.vcf.gz
gs://pepper-deepvariant/mobinasri/pangenome_paper/external_callers_element/GATK/regular_v4.6/HG003_GATK_v4.6_regular_mode/analysis/WholeGenomeGermlineSingleSample_outputs/HG003_Element_1000bp_ins_GATK_Regular_v4.6.rb.g.PASS.vcf.gz.tbi
gs://pepper-deepvariant/mobinasri/pangenome_paper/external_callers_element/GATK/regular_v4.6/HG003_GATK_v4.6_regular_mode/analysis/WholeGenomeGermlineSingleSample_outputs/HG003_Element_1000bp_ins_GATK_Regular_v4.6.rb.g.vcf.gz
gs://pepper-deepvariant/mobinasri/pangenome_paper/external_callers_element/GATK/regular_v4.6/HG003_GATK_v4.6_regular_mode/analysis/WholeGenomeGermlineSingleSample_outputs/HG003_Element_1000bp_ins_GATK_Regular_v4.6.rb.g.vcf.gz.tbi
gs://pepper-deepvariant/mobinasri/pangenome_paper/external_callers_element/octopus_v0.7.4/HG003.element.cloudbreak.1000bp_ins.bwa_mem.hg38.octopus.no_forest.vcf.gz
gs://pepper-deepvariant/mobinasri/pangenome_paper/external_callers_element/octopus_v0.7.4/HG003.element.cloudbreak.1000bp_ins.bwa_mem.hg38.octopus.no_forest.vcf.gz.tbi
gs://pepper-deepvariant/mobinasri/pangenome_paper/external_callers_element/octopus_v0.7.4/HG003.element.cloudbreak.1000bp_ins.vg.grch38.octopus.no_forest.vcf.gz
gs://pepper-deepvariant/mobinasri/pangenome_paper/external_callers_element/octopus_v0.7.4/HG003.element.cloudbreak.1000bp_ins.vg.grch38.octopus.no_forest.vcf.gz.tbi
```
