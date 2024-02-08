## Investigating errors in polishing HPRC y2 assemblies
### Comment 1 : 02/07/2024

After Mira noticed QV drops after polishing some HPRC-Y2 samples we decided to investigate what is causing that. Here is Mira's slide deck (look at page 2 and 3)

https://docs.google.com/presentation/d/1ntvL7ozycyiJV5BLPv9Cd-EHkAet0nUh6R5TW9pwhj0/edit#slide=id.g2b781253d9e_0_99

Here are the steps we decided to go through:
- Map polished assemblies to raw assemblies
- Use that mapping to project the FP kmer blocks from polished to raw assemblies
- Now that we have FP kmers before and after polishing in the same coordinates I'll look for the places where polishing is generating FP kmers
- Select 10 random cases with induced FP kmers ( maybe prioritizing the places with higher density of induced FP kmers)
- Map hap1 and hap2 of the raw assemblies against each other then find the equivalent blocks of the selected cases on the other haplotype. This way we can examine both haplotypes at the same time.
- Look at PHARAOH HiFi alignments and DeepPolisher vcf at those locations and see if there is something wrong with the alignments or DeepPolisher
- Mira also suggested mapping short reads to both haplotypes (all-to-one) and viewing those along with PHARAOH HiFi alignments

We selected the sample `HG04115` for closer investigation since it was one of the samples with highest QV drop both whole genome and in confident regions.

#### Made a wdl for aligning short reads with bwa and added to Flagger repo
WDL is added here:

https://github.com/mobinasri/flagger/blob/dev-0.3.0/wdls/tasks/alignment/bwa.wdl

#### Where the files are located 

Here are the places where the files are located on pheonix
```
# polished assemblies, pharaoh alignments are here
ls /private/groups/hprc/polishing/batch3/HG04115/hprc_DeepPolisher_outputs/
HG04115_Hap1.polished.fasta	 HG04115_Hap2.polished.fasta	  HG04115.hifi.to.diploid.asm.PHARAOH.bam      polisher_output.vcf.gz
HG04115_Hap1.polished.fasta.fai  HG04115_Hap2.polished.fasta.fai  HG04115.hifi.to.diploid.asm.PHARAOH.bam.bai

# raw assemblies are here
ls /private/groups/hprc/assembly/batch2/HG04115/analysis/assembly/
HG03894.yak  HG03896.yak  HG04115.binFiles.tar.gz  HG04115.mat.contig_gfa.tar.gz  HG04115.mat.fa.gz  HG04115.pat.contig_gfa.tar.gz  HG04115.pat.fa.gz  HG04115.raw_unitig_gfa.tar.gz
```

I put the input files here
```
/private/groups/patenlab/masri/hprc/polishing/investigating_Y2_results/HG04115$ tree
.
├── assemblies
│   ├── HG04115_Hap1.polished.fasta
│   ├── HG04115_Hap1.polished.fasta.fai
│   ├── HG04115_Hap2.polished.fasta
│   ├── HG04115_Hap2.polished.fasta.fai
│   ├── HG04115.mat.fa
│   ├── HG04115.mat.fa.fai
│   ├── HG04115.mat.fa.gz
│   ├── HG04115.pat.fa
│   ├── HG04115.pat.fa.fai
│   └── HG04115.pat.fa.gz
└── short_reads
    └── alignment
        ├── data_table.csv
        └── input_mapping.csv
```

inside `data_table.csv`:
```
sample_id,suffix,read_cram,assembly_fasta_gz
HG04115_pat_raw,"short_reads","s3://human-pangenomics/submissions/325b4b1c-9f20-49be-b03a-596da89c466e--1000G_CHILDREN/HG04115/1000G_data/HG04115.final.cram","/private/groups/patenlab/masri/hprc/polishing/investigating_Y2_results/HG04115/assemblies/HG04115.pat.fa.gz"
HG04115_mat_raw,"short_reads","s3://human-pangenomics/submissions/325b4b1c-9f20-49be-b03a-596da89c466e--1000G_CHILDREN/HG04115/1000G_data/HG04115.final.cram","/private/groups/patenlab/masri/hprc/polishing/investigating_Y2_results/HG04115/assemblies/HG04115.mat.fa.gz"
```

inside `input_mapping.csv`
```
input,type,value
"bwaAlignment.referenceFasta",scalar,"/private/groups/patenlab/masri/hprc/polishing/HG002/assemblies/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
"bwaAlignment.buildBwaIndex.memSize",scalar,16
"bwaAlignment.BwaAlignment.threadCount",scalar,32
"bwaAlignment.BwaAlignment.bwaParams",scalar,"-p"
"bwaAlignment.sampleName",scalar,$input.sample_id
"bwaAlignment.buildBwaIndex.threadCount",scalar,4
"bwaAlignment.assembly",scalar,$input.assembly_fasta_gz
"bwaAlignment.cramFile",scalar,$input.read_cram
"bwaAlignment.suffix",scalar,$input.suffix
```

#### Run short read alignment (all to one)

Prepare json
```
USER_NAME="masri"
EMAIL="masri@ucsc.edu"

WDL_PATH="/private/groups/patenlab/masri/apps/flagger/wdls/tasks/alignment/bwa.wdl"
WDL_NAME=$(basename ${WDL_PATH%%.wdl})

INPUT_MAPPING_CSV="/private/groups/patenlab/masri/hprc/polishing/investigating_Y2_results/HG04115/short_reads/alignment/input_mapping.csv"
INPUT_DATA_TABLE_CSV="/private/groups/patenlab/masri/hprc/polishing/investigating_Y2_results/HG04115/short_reads/alignment/data_table.csv"

WORKING_DIR="/private/groups/patenlab/masri/hprc/polishing/investigating_Y2_results/HG04115/short_reads/alignment"

mkdir -p ${WORKING_DIR}
cd ${WORKING_DIR}

mkdir ${WDL_NAME}_output_jsons
mkdir ${WDL_NAME}_input_jsons


LAUNCH_FROM_TABLE_PY="/private/groups/patenlab/masri/apps/hprc_intermediate_assembly/hpc/launch_from_table.py"
LAUNCH_WORKFLOW_JOB_ARRAY_BASH="/private/groups/patenlab/masri/apps/hprc_intermediate_assembly/hpc/launch_workflow_job_array_single_machine.sh"

cd ${WORKING_DIR}/${WDL_NAME}_input_jsons
python3  ${LAUNCH_FROM_TABLE_PY} \
    --data_table ${INPUT_DATA_TABLE_CSV} \
    --field_mapping ${INPUT_MAPPING_CSV} \
    --workflow_name ${WDL_NAME}
```

Look at input json files
```
cat HG04115_mat_raw_bwa.json
{
  "bwaAlignment.referenceFasta": "/private/groups/patenlab/masri/hprc/polishing/HG002/assemblies/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
  "bwaAlignment.buildBwaIndex.memSize": 16,
  "bwaAlignment.BwaAlignment.threadCount": 32,
  "bwaAlignment.BwaAlignment.bwaParams": "-p",
  "bwaAlignment.sampleName": "HG04115_mat_raw",
  "bwaAlignment.buildBwaIndex.threadCount": 4,
  "bwaAlignment.assembly": "/private/groups/patenlab/masri/hprc/polishing/investigating_Y2_results/HG04115/assemblies/HG04115.mat.fa.gz",
  "bwaAlignment.cramFile": "s3://human-pangenomics/submissions/325b4b1c-9f20-49be-b03a-596da89c466e--1000G_CHILDREN/HG04115/1000G_data/HG04115.final.cram",
  "bwaAlignment.suffix": "short_reads"
}


cat HG04115_pat_raw_bwa.json
{
  "bwaAlignment.referenceFasta": "/private/groups/patenlab/masri/hprc/polishing/HG002/assemblies/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
  "bwaAlignment.buildBwaIndex.memSize": 16,
  "bwaAlignment.BwaAlignment.threadCount": 32,
  "bwaAlignment.BwaAlignment.bwaParams": "-p",
  "bwaAlignment.sampleName": "HG04115_pat_raw",
  "bwaAlignment.buildBwaIndex.threadCount": 4,
  "bwaAlignment.assembly": "/private/groups/patenlab/masri/hprc/polishing/investigating_Y2_results/HG04115/assemblies/HG04115.pat.fa.gz",
  "bwaAlignment.cramFile": "s3://human-pangenomics/submissions/325b4b1c-9f20-49be-b03a-596da89c466e--1000G_CHILDREN/HG04115/1000G_data/HG04115.final.cram",
  "bwaAlignment.suffix": "short_reads"

```

Run the job
```
cd ${WORKING_DIR}
cd ${WDL_NAME}_output_jsons
INPUT_JSON_DIR=${WORKING_DIR}/${WDL_NAME}_input_jsons
mkdir -p ${WDL_NAME}_logs

sbatch      --job-name=${WDL_NAME}_${USERNAME} \
            --cpus-per-task=32 \
            --mem=64G \
            --mail-user=${EMAIL} \
            --output=${WDL_NAME}_logs/${WDL_NAME}_%A_%a.log \
            --array=1-2%2  \
            --partition=medium  \
            --time=12:00:00 \
            ${LAUNCH_WORKFLOW_JOB_ARRAY_BASH} \
            --wdl ${WDL_PATH} \
            --sample_csv  ${INPUT_DATA_TABLE_CSV} \
            --input_json_dir ${INPUT_JSON_DIR}
```

### Comment 2 : 02/07/2024

#### Make asm2asm alignments

We need three assembly to assembly alignments for this analysis:
- polished pat to raw pat
- polished mat to raw mat
- raw mat to raw pat

`data_table.csv`
```
sample_id,ref_assembly_fasta_gz,query_assembly_fasta_gz
HG04115_pat_polished_to_raw,"/private/groups/patenlab/masri/hprc/polishing/investigating_Y2_results/HG04115/assemblies/HG04115.pat.fa.gz","/private/groups/patenlab/masri/hprc/polishing/investigating_Y2_results/HG04115/assemblies/HG04115_Hap1.polished.fa.gz"
HG04115_mat_polished_to_raw,"/private/groups/patenlab/masri/hprc/polishing/investigating_Y2_results/HG04115/assemblies/HG04115.mat.fa.gz","/private/groups/patenlab/masri/hprc/polishing/investigating_Y2_results/HG04115/assemblies/HG04115_Hap2.polished.fa.gz"
HG04115_polished_mat_to_pat,"/private/groups/patenlab/masri/hprc/polishing/investigating_Y2_results/HG04115/assemblies/HG04115.pat.fa.gz","/private/groups/patenlab/masri/hprc/polishing/investigating_Y2_results/HG04115/assemblies/HG04115.mat.fa.gz"
```

`input_mapping.csv`
```
input,type,value
asm2asmAlignment.queryAssemblyFastaGz,scalar,$input.query_assembly_fasta_gz
asm2asmAlignment.alignmentBam.threadCount,scalar,32
asm2asmAlignment.refAssemblyFastaGz,scalar,$input.ref_assembly_fasta_gz
asm2asmAlignment.alignmentBam.memSize,scalar,64
asm2asmAlignment.alignmentBam.diskSize,scalar,64
asm2asmAlignment.preset,scalar,"asm5"
asm2asmAlignment.suffix,scalar,$input.sample_id
asm2asmAlignment.alignmentBam.options,scalar,"--eqx --cs -L"
asm2asmAlignment.aligner,scalar,"minimap2"
```

Make json files:
```
USER_NAME="masri"
EMAIL="masri@ucsc.edu"

WDL_PATH="/private/groups/patenlab/masri/apps/flagger/wdls/tasks/alignment/asm2asm_aligner.wdl"
WDL_NAME=$(basename ${WDL_PATH%%.wdl})

INPUT_MAPPING_CSV="/private/groups/patenlab/masri/hprc/polishing/investigating_Y2_results/HG04115/asm_alignment/input_mapping.csv"
INPUT_DATA_TABLE_CSV="/private/groups/patenlab/masri/hprc/polishing/investigating_Y2_results/HG04115/asm_alignment/data_table.csv"

WORKING_DIR="/private/groups/patenlab/masri/hprc/polishing/investigating_Y2_results/HG04115/asm_alignment"

mkdir -p ${WORKING_DIR}
cd ${WORKING_DIR}

mkdir ${WDL_NAME}_output_jsons
mkdir ${WDL_NAME}_input_jsons


LAUNCH_FROM_TABLE_PY="/private/groups/patenlab/masri/apps/hprc_intermediate_assembly/hpc/launch_from_table.py"
LAUNCH_WORKFLOW_JOB_ARRAY_BASH="/private/groups/patenlab/masri/apps/hprc_intermediate_assembly/hpc/launch_workflow_job_array_single_machine.sh"

cd ${WORKING_DIR}/${WDL_NAME}_input_jsons
python3  ${LAUNCH_FROM_TABLE_PY} \
    --data_table ${INPUT_DATA_TABLE_CSV} \
    --field_mapping ${INPUT_MAPPING_CSV} \
    --workflow_name ${WDL_NAME}
```

Look at json files:
```
cat HG04115_mat_polished_to_raw_asm2asm_aligner.json
{
  "asm2asmAlignment.queryAssemblyFastaGz": "/private/groups/patenlab/masri/hprc/polishing/investigating_Y2_results/HG04115/assemblies/HG04115_Hap2.polished.fa.gz",
  "asm2asmAlignment.alignmentBam.threadCount": 32,
  "asm2asmAlignment.refAssemblyFastaGz": "/private/groups/patenlab/masri/hprc/polishing/investigating_Y2_results/HG04115/assemblies/HG04115.mat.fa.gz",
  "asm2asmAlignment.alignmentBam.memSize": 64,
  "asm2asmAlignment.alignmentBam.diskSize": 64,
  "asm2asmAlignment.preset": "asm5",
  "asm2asmAlignment.suffix": "HG04115_mat_polished_to_raw",
  "asm2asmAlignment.alignmentBam.options": "--eqx --cs -L",
  "asm2asmAlignment.aligner": "minimap2"
}


cat HG04115_pat_polished_to_raw_asm2asm_aligner.json
{
  "asm2asmAlignment.queryAssemblyFastaGz": "/private/groups/patenlab/masri/hprc/polishing/investigating_Y2_results/HG04115/assemblies/HG04115_Hap1.polished.fa.gz",
  "asm2asmAlignment.alignmentBam.threadCount": 32,
  "asm2asmAlignment.refAssemblyFastaGz": "/private/groups/patenlab/masri/hprc/polishing/investigating_Y2_results/HG04115/assemblies/HG04115.pat.fa.gz",
  "asm2asmAlignment.alignmentBam.memSize": 64,
  "asm2asmAlignment.alignmentBam.diskSize": 64,
  "asm2asmAlignment.preset": "asm5",
  "asm2asmAlignment.suffix": "HG04115_pat_polished_to_raw",
  "asm2asmAlignment.alignmentBam.options": "--eqx --cs -L",
  "asm2asmAlignment.aligner": "minimap2"
}


cat HG04115_polished_mat_to_pat_asm2asm_aligner.json
{
  "asm2asmAlignment.queryAssemblyFastaGz": "/private/groups/patenlab/masri/hprc/polishing/investigating_Y2_results/HG04115/assemblies/HG04115.mat.fa.gz",
  "asm2asmAlignment.alignmentBam.threadCount": 32,
  "asm2asmAlignment.refAssemblyFastaGz": "/private/groups/patenlab/masri/hprc/polishing/investigating_Y2_results/HG04115/assemblies/HG04115.pat.fa.gz",
  "asm2asmAlignment.alignmentBam.memSize": 64,
  "asm2asmAlignment.alignmentBam.diskSize": 64,
  "asm2asmAlignment.preset": "asm5",
  "asm2asmAlignment.suffix": "HG04115_polished_mat_to_pat",
  "asm2asmAlignment.alignmentBam.options": "--eqx --cs -L",
  "asm2asmAlignment.aligner": "minimap2"
}
```

Run the job
```
cd ${WORKING_DIR}
cd ${WDL_NAME}_output_jsons
INPUT_JSON_DIR=${WORKING_DIR}/${WDL_NAME}_input_jsons
mkdir -p ${WDL_NAME}_logs

sbatch      --job-name=${WDL_NAME}_${USERNAME} \
            --cpus-per-task=32 \
            --mem=64G \
            --mail-user=${EMAIL} \
            --output=${WDL_NAME}_logs/${WDL_NAME}_%A_%a.log \
            --array=1-3%3  \
            --partition=medium  \
            --time=12:00:00 \
            ${LAUNCH_WORKFLOW_JOB_ARRAY_BASH} \
            --wdl ${WDL_PATH} \
            --sample_csv  ${INPUT_DATA_TABLE_CSV} \
            --input_json_dir ${INPUT_JSON_DIR}
```
