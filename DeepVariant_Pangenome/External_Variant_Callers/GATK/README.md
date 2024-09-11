# Run GATK on HG003 Novaseq 35x in three modes

This doc contains wdls for running GATK. The wdls were updated last time with GATKv4.5
https://broadinstitute.github.io/warp/docs/Pipelines/Whole_Genome_Germline_Single_Sample_Pipeline/README/

Here is the main wdl we have to run to obtain VCF from uBam.
https://github.com/broadinstitute/warp/blob/develop/pipelines/broad/dna_seq/germline/single_sample/wgs/WholeGenomeGermlineSingleSample.wdl

Based on documentaion we can run this workflow in three modes:
1. **Regaular GATK (without DRAGEN modules)**: Use default values
2. **GATK-DRAGEN (functional equivalence mode)**: Only by setting `dragen_functional_equivalence_mode = true` In this mode the output of pipeline will be functionally equivalent to the one produced with the DRAGEN hardware.
3. **GATK-DRAGEN (maximum quality mode)**: Only by setting `dragen_maximum_quality_mode = true` In this mode the pipeline use DRAGEN mapper and variant calling steps however the final results are **NOT** functionally equivalent to the one produced with the DRAGEN hardware.

### Note: I will set `use_gatk3_haplotype_caller = false` in all modes since I want to use GATK-4.5.

## Get Data
```
mkdir -p /private/groups/patenlab/masri/internship/external_callers/data
cd /private/groups/patenlab/masri/internship/external_callers/data

wget https://storage.googleapis.com/brain-genomics-public/research/sequencing/fastq/novaseq/wgs_pcr_free/35x/HG003.novaseq.pcr-free.35x.R1.fastq.gz
wget https://storage.googleapis.com/brain-genomics-public/research/sequencing/fastq/novaseq/wgs_pcr_free/35x/HG003.novaseq.pcr-free.35x.R2.fastq.gz
```

## Get WDL

```

```


