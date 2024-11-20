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
for wdl_file in $(find ./tasks/broad/ ./pipelines/broad/ | grep ".wdl$");
do
  awk '
  /CODE/ && /python/ {
      match($0, /^ */);            # Count leading spaces on the python line
      spaces = RLENGTH;            # Store the number of leading spaces
      in_block = 1;                # Set flag to indicate we are in the block between python \<\< CODE and CODE
  }
  in_block && /CODE/ && !/python/ {
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
  }' ${wdl_file} > temp && mv temp ${wdl_file}
done
```

### Increase `memory_multiplier` for HaplotypeCaller_GATK4_VCF
```
task HaplotypeCaller_GATK4_VCF {
  input {
    Int memory_multiplier = 2
  }
```

### It seems like Toil does not evaluate the outputs of a task sequentially (?)

I got this bug in the output log of Toil while running GATK in the regular mode.
```
	[2024-09-16T04:03:08-0700] [MainThread] [E] [toil.wdl.wdltoil] Expression evaluation failed for is_outlier_data: duplication_rate > max_duplication_in_reasonable_sample || chimerism_rate > max_chimerism_in_reasonable_sample
	Traceback (most recent call last):
	  File "/private/home/masri/.local/lib/python3.10/site-packages/WDL/Expr.py", line 123, in eval
	    ans = self._eval(env, stdlib)
	  File "/private/home/masri/.local/lib/python3.10/site-packages/WDL/Expr.py", line 863, in _eval
	    return env[self.name]
	  File "/private/home/masri/.local/lib/python3.10/site-packages/WDL/Env.py", line 127, in __getitem__
	    return self.resolve(name)
	  File "/private/home/masri/.local/lib/python3.10/site-packages/WDL/Env.py", line 114, in resolve
	    return self.resolve_binding(name).value
	  File "/private/home/masri/.local/lib/python3.10/site-packages/WDL/Env.py", line 106, in resolve_binding
	    raise KeyError()
	KeyError
	
	The above exception was the direct cause of the following exception:
```

#### I changed this code in `tasks/broad/Qc.wdl`
```
  output {
    Float duplication_rate = read_float("duplication_value.txt")
    Float chimerism_rate = read_float("chimerism_value.txt")
    Boolean is_outlier_data = duplication_rate > max_duplication_in_reasonable_sample || chimerism_rate > max_chimerism_in_reasonable_sample
  }
}
```
#### To
```
  output {
    Float duplication_rate = read_float("duplication_value.txt")
    Float chimerism_rate = read_float("chimerism_value.txt")
    Boolean is_outlier_data = (read_float("duplication_value.txt") > max_duplication_in_reasonable_sample) || (read_float("chimerism_value.txt") > max_chimerism_in_reasonable_sample)
  }
```

### Fix allocating high memory to Java. (Thanks to Adam for catching this)

I tried to run the workflows multiple times and they keep failing because of running out of memory even after increasing `memory_multiplier`. The output of `toil stats` was not reporting a significant amount of memory usage for any task (espeically for `HaplotypeCallerGATK4`). For example here it reports a total memory of only 1.5GB. 
```
 WholeGenomeGermlineSingleSample.BamToGvcf.HaplotypeCallerGATK4
    Total Cores: 6.0
    Count |                                        Time* |                                             Clock |                                               Wait |                                   Memory 
        n |      min     med*     ave      max     total |      min       med       ave       max      total |         min     med        ave     max       total |      min     med     ave     max   total 
        6 |     0.06     0.07    0.07     0.07      0.39 |     0.06      0.07      0.07      0.07       0.39 |       -0.00    0.00       0.00    0.00        0.00 |  241336K 241616K 241693K 242024K1450160K
```

I asked Adam about it and he pointed me to this line:
https://github.com/broadinstitute/warp/blob/fe21c077ffc9c48ae7368e44bb7f6e88ce303b24/tasks/broad/GermlineVariantDiscovery.wdl#L124

It seems like that this task is using `free` to get the total available memory but `free` does not consider the memory limitation enforced by Toil. Using `free` on a Pheonix node took the full 2TB of available memory. I will use `memory_size_mb` instead of `free` so that hopefully that out-of-memory issue is fixed.

#### I replaced this line
```
    available_memory_mb=$(free -m | awk '/^Mem/ {print $2}')
```
#### with this
```
    #$(free -m | awk '/^Mem/ {print $2}')
    available_memory_mb=~{memory_size_mb} 
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

### run GATK in regular mode (no DRAGEN and `use_gatk3_haplotype_caller = true`)

```
cd /private/groups/patenlab/masri/internship/external_callers/results/GATK/regular_v3

WDL_PATH="/private/groups/patenlab/masri/internship/external_callers/apps/warp/pipelines/broad/dna_seq/germline/single_sample/wgs/WholeGenomeGermlineSingleSample.wdl"
INPUT_JSON_PATH="/private/groups/patenlab/masri/internship/external_callers/results/GATK/regular_v3/gatk_v3_regular.json"
SAMPLE_NAME="HG003_GATK_v3_regular_mode"
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
cd /private/groups/patenlab/masri/internship/external_callers/results/GATK/regular_v4.6

WDL_PATH="/private/groups/patenlab/masri/internship/external_callers/apps/warp/pipelines/broad/dna_seq/germline/single_sample/wgs/WholeGenomeGermlineSingleSample.wdl"
INPUT_JSON_PATH="/private/groups/patenlab/masri/internship/external_callers/results/GATK/regular_v4.6/gatk_v4.6_regular.json"
SAMPLE_NAME="HG003_GATK_v4.6_regular_mode"
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
INPUT_JSON_PATH="/private/groups/patenlab/masri/internship/external_callers/results/GATK/dragen_maximum_quality/dragen_mode_maximum_quality.json"
SAMPLE_NAME="HG003_GATK_v4.6_dragen_maximum_quality"
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

## Results

### Initial stats for `regular` mode (v4.6)
```
cd /private/groups/patenlab/masri/internship/external_callers/results/GATK/regular_v4.6/HG003_GATK_v4.6_regular_mode/analysis/WholeGenomeGermlineSingleSample_outputs

bcftools view -Oz \
    -f 'PASS,.' \
    HG003_Novaseq_35x_GATK_Regular.rb.g.vcf.gz \
    -o HG003_Novaseq_35x_GATK_Regular.rb.g.PASS.vcf.gz

tabix HG003_Novaseq_35x_GATK_Regular.rb.g.PASS.vcf.gz
```
```
## bcftools stats

# SN	[2]id	[3]key	[4]value
SN	0	number of samples:	1
SN	0	number of records:	45055671
SN	0	number of no-ALTs:	39974778
SN	0	number of SNPs:	4121887
SN	0	number of MNPs:	0
SN	0	number of indels:	964861
SN	0	number of others:	0
SN	0	number of multiallelic sites:	5080893
SN	0	number of multiallelic SNP sites:	4106810
```

### Run hap.py for `regular` mode (v4.6)

```
QUERY_VCF="/private/groups/patenlab/masri/internship/external_callers/results/GATK/regular_v4.6/HG003_GATK_v4.6_regular_mode/analysis/WholeGenomeGermlineSingleSample_outputs/HG003_Novaseq_35x_GATK_Regular.rb.g.PASS.vcf.gz"
TRUTH_VCF="/private/groups/patenlab/masri/internship/external_callers/data/giab/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
CONF_BED="/private/groups/patenlab/masri/internship/external_callers/data/giab/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
REF_FASTA="/private/groups/patenlab/masri/internship/external_callers/data/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"
OUT_DIR="/private/groups/patenlab/masri/internship/external_callers/results/GATK/regular_v4.6/happy"
OUT_PREFIX="${OUT_DIR}/regular_v4.6"
BASE_DIR="/private/groups/patenlab/masri/internship/external_callers"

RUN_HAPPY_BASH="/private/groups/patenlab/masri/internship/external_callers/apps/run_happy.sh"

EMAIL="masri@ucsc.edu"

mkdir -p ${OUT_DIR}
sbatch      --job-name=run_happy_regular_v4.6 \
            --cpus-per-task=16 \
            --mem=16G \
            --mail-user=${EMAIL} \
            --output=${OUT_DIR}/run_happy_%A_%a.log \
            --partition=long  \
            --time=6:00:00 \
            ${RUN_HAPPY_BASH} \
            --truth ${TRUTH_VCF} \
            --query ${QUERY_VCF} \
            --bed ${CONF_BED} \
            --ref ${REF_FASTA} \
            --output_prefix ${OUT_PREFIX} \
            --base_dir ${BASE_DIR}
```

### Initial stats for `regular` mode (v3)
```
cd /private/groups/patenlab/masri/internship/external_callers/results/GATK/regular_v3/HG003_GATK_v3_regular_mode/analysis/WholeGenomeGermlineSingleSample_outputs

bcftools view -Oz \
    -f 'PASS,.' \
    HG003_Novaseq_35x_GATK_Regular.rb.g.vcf.gz \
    -o HG003_Novaseq_35x_GATK_Regular.rb.g.PASS.vcf.gz

tabix HG003_Novaseq_35x_GATK_Regular.rb.g.PASS.vcf.gz
```

```
## bcftools stats

# SN	[2]id	[3]key	[4]value
SN	0	number of samples:	1
SN	0	number of records:	48731937
SN	0	number of no-ALTs:	43664349
SN	0	number of SNPs:	4111882
SN	0	number of MNPs:	0
SN	0	number of indels:	961392
SN	0	number of others:	0
SN	0	number of multiallelic sites:	5067588
SN	0	number of multiallelic SNP sites:	4106196
```
### Run hap.py for `regular` mode (v3)
```
QUERY_VCF="/private/groups/patenlab/masri/internship/external_callers/results/GATK/regular_v3/HG003_GATK_v3_regular_mode/analysis/WholeGenomeGermlineSingleSample_outputs/HG003_Novaseq_35x_GATK_Regular.rb.g.PASS.vcf.gz"
TRUTH_VCF="/private/groups/patenlab/masri/internship/external_callers/data/giab/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
CONF_BED="/private/groups/patenlab/masri/internship/external_callers/data/giab/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
REF_FASTA="/private/groups/patenlab/masri/internship/external_callers/data/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"
OUT_DIR="/private/groups/patenlab/masri/internship/external_callers/results/GATK/regular_v3/happy"
OUT_PREFIX="${OUT_DIR}/regular_v3"
BASE_DIR="/private/groups/patenlab/masri/internship/external_callers"

RUN_HAPPY_BASH="/private/groups/patenlab/masri/internship/external_callers/apps/run_happy.sh"

EMAIL="masri@ucsc.edu"

mkdir -p ${OUT_DIR}
sbatch      --job-name=run_happy_regular_v4.6 \
            --cpus-per-task=16 \
            --mem=16G \
            --mail-user=${EMAIL} \
            --output=${OUT_DIR}/run_happy_%A_%a.log \
            --partition=long  \
            --time=6:00:00 \
            ${RUN_HAPPY_BASH} \
            --truth ${TRUTH_VCF} \
            --query ${QUERY_VCF} \
            --bed ${CONF_BED} \
            --ref ${REF_FASTA} \
            --output_prefix ${OUT_PREFIX} \
            --base_dir ${BASE_DIR}
```

### Initial stats for `dragen_functional_equivalence` mode

```
cd /private/groups/patenlab/masri/internship/external_callers/results/GATK/dragen_functional_equivalence/HG003_GATK_v4.6_dragen_functional_equivalence/analysis/WholeGenomeGermlineSingleSample_outputs

bcftools view -Oz \
    -f 'PASS,.' \
    HG003_Novaseq_35x_GATK_Dragen_functional_equivalence.hard-filtered.rb.g.vcf.gz \
    -o HG003_Novaseq_35x_GATK_Dragen_functional_equivalence.hard-filtered.rb.g.PASS.vcf.gz

tabix HG003_Novaseq_35x_GATK_Dragen_functional_equivalence.hard-filtered.rb.g.PASS.vcf.gz
```
```
## bcftools stats

# SN	[2]id	[3]key	[4]value
SN	0	number of samples:	1
SN	0	number of records:	5410903
SN	0	number of no-ALTs:	488062
SN	0	number of SNPs:	3970930
SN	0	number of MNPs:	0
SN	0	number of indels:	957573
SN	0	number of others:	0
SN	0	number of multiallelic sites:	4922841
SN	0	number of multiallelic SNP sites:	3965268
```
### Run hap.py for `dragen_functional_equivalence` mode

```
QUERY_VCF="/private/groups/patenlab/masri/internship/external_callers/results/GATK/dragen_functional_equivalence/HG003_GATK_v4.6_dragen_functional_equivalence/analysis/WholeGenomeGermlineSingleSample_outputs/HG003_Novaseq_35x_GATK_Dragen_functional_equivalence.hard-filtered.rb.g.PASS.vcf.gz"
TRUTH_VCF="/private/groups/patenlab/masri/internship/external_callers/data/giab/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
CONF_BED="/private/groups/patenlab/masri/internship/external_callers/data/giab/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
REF_FASTA="/private/groups/patenlab/masri/internship/external_callers/data/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"
OUT_DIR="/private/groups/patenlab/masri/internship/external_callers/results/GATK/dragen_functional_equivalence/happy"
OUT_PREFIX="${OUT_DIR}/dragen_functional_equivalence"
BASE_DIR="/private/groups/patenlab/masri/internship/external_callers"

RUN_HAPPY_BASH="/private/groups/patenlab/masri/internship/external_callers/apps/run_happy.sh"

EMAIL="masri@ucsc.edu"

mkdir -p ${OUT_DIR}
sbatch      --job-name=run_happy_dragen_functional_equivalence \
            --cpus-per-task=16 \
            --mem=16G \
            --mail-user=${EMAIL} \
            --output=${OUT_DIR}/run_happy_%A_%a.log \
            --partition=long  \
            --time=6:00:00 \
            ${RUN_HAPPY_BASH} \
            --truth ${TRUTH_VCF} \
            --query ${QUERY_VCF} \
            --bed ${CONF_BED} \
            --ref ${REF_FASTA} \
            --output_prefix ${OUT_PREFIX} \
            --base_dir ${BASE_DIR}
```

### Initial stats for `dragen_maximum_quality` mode

```
cd /private/groups/patenlab/masri/internship/external_callers/results/GATK/dragen_maximum_quality/HG003_GATK_v4.6_dragen_maximum_quality/analysis/WholeGenomeGermlineSingleSample_outputs

bcftools view -Oz \
    -f 'PASS,.' \
    HG003_Novaseq_35x_GATK_Dragen_max_qual.hard-filtered.rb.g.vcf.gz \
    -o HG003_Novaseq_35x_GATK_Dragen_max_qual.hard-filtered.rb.g.PASS.vcf.gz

tabix HG003_Novaseq_35x_GATK_Dragen_max_qual.hard-filtered.rb.g.PASS.vcf.gz
```
```
## bcftools stats

# SN	[2]id	[3]key	[4]value
SN	0	number of samples:	1
SN	0	number of records:	5412432
SN	0	number of no-ALTs:	487153
SN	0	number of SNPs:	3972167
SN	0	number of MNPs:	0
SN	0	number of indels:	958344
SN	0	number of others:	0
SN	0	number of multiallelic sites:	4925279
SN	0	number of multiallelic SNP sites:	3957532

```

### Run hap.py for `dragen_maximum_quality` mode

```
QUERY_VCF="/private/groups/patenlab/masri/internship/external_callers/results/GATK/dragen_maximum_quality/HG003_GATK_v4.6_dragen_maximum_quality/analysis/WholeGenomeGermlineSingleSample_outputs/HG003_Novaseq_35x_GATK_Dragen_max_qual.hard-filtered.rb.g.PASS.vcf.gz"
TRUTH_VCF="/private/groups/patenlab/masri/internship/external_callers/data/giab/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
CONF_BED="/private/groups/patenlab/masri/internship/external_callers/data/giab/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
REF_FASTA="/private/groups/patenlab/masri/internship/external_callers/data/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"
OUT_DIR="/private/groups/patenlab/masri/internship/external_callers/results/GATK/dragen_maximum_quality/happy"
OUT_PREFIX="${OUT_DIR}/dragen_maximum_quality"
BASE_DIR="/private/groups/patenlab/masri/internship/external_callers"

RUN_HAPPY_BASH="/private/groups/patenlab/masri/internship/external_callers/apps/run_happy.sh"

EMAIL="masri@ucsc.edu"

mkdir -p ${OUT_DIR}
sbatch      --job-name=run_happy_dragen_maximum_quality \
            --cpus-per-task=16 \
            --mem=16G \
            --mail-user=${EMAIL} \
            --output=${OUT_DIR}/run_happy_%A_%a.log \
            --partition=long  \
            --time=6:00:00 \
            ${RUN_HAPPY_BASH} \
            --truth ${TRUTH_VCF} \
            --query ${QUERY_VCF} \
            --bed ${CONF_BED} \
            --ref ${REF_FASTA} \
            --output_prefix ${OUT_PREFIX} \
            --base_dir ${BASE_DIR}
```


