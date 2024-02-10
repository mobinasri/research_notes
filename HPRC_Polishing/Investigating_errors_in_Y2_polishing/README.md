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

Get stats for the bam file
```
samtools flagstats HG04115_mat_raw_short_reads.sorted.bam
726438176 + 0 in total (QC-passed reads + QC-failed reads)
721374894 + 0 primary
0 + 0 secondary
5063282 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
724289710 + 0 mapped (99.70% : N/A)
719226428 + 0 primary mapped (99.70% : N/A)
2666994 + 0 paired in sequencing
1333497 + 0 read1
1333497 + 0 read2
739770 + 0 properly paired (27.74% : N/A)
1144398 + 0 with itself and mate mapped
714647 + 0 singletons (26.80% : N/A)
321688 + 0 with mate mapped to a different chr
164568 + 0 with mate mapped to a different chr (mapQ>=5)
```

#### Fix the issue with low percentage of properly paired reads
This problem arised because of not sorting reads by name while extracting reads from cram and saving them in fastq. Becasuse of this the paired-end reads may not end up being beside each other in the output fastq, which is neccessary for using `bwa mem` with `-p` parameter.
I sorted reads by name in cram and them save them in fastq. It fixed the issue and now the percentage of the properly paired reads is much higher than before.
```
samtools flagstats HG04115_mat_raw_short_reads.sorted.bam
726380027 + 0 in total (QC-passed reads + QC-failed reads)
721374894 + 0 primary
0 + 0 secondary
5005133 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
725024515 + 0 mapped (99.81% : N/A)
720019382 + 0 primary mapped (99.81% : N/A)
721374894 + 0 paired in sequencing
360687447 + 0 read1
360687447 + 0 read2
709529170 + 0 properly paired (98.36% : N/A)
719217508 + 0 with itself and mate mapped
801874 + 0 singletons (0.11% : N/A)
8399022 + 0 with mate mapped to a different chr
4398355 + 0 with mate mapped to a different chr (mapQ>=5)
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


### Comment 3 : 02/08/2024
#### Convert BAM to PAF

My projection script works with paf format so I have to convert all bam files to paf
```
# get some resources on pheonix
salloc  -c 8 --mem 32G bash -c 'ssh -Y $(scontrol show hostnames | head -n 1)'

# run the docker image with paftools
docker run -u$(id -u):$(id -g) --rm -it -v/private/groups/patenlab/masri/:/private/groups/patenlab/masri/ mobinasri/long_read_aligner:v0.3.3

cd /private/groups/patenlab/masri/hprc/polishing/investigating_Y2_results/HG04115/asm_alignment/asm2asm_aligner_output_jsons/

## convert to paf

# for polished to raw maternal
k8 /home/apps/minimap2-2.26/misc/paftools.js sam2paf -p <(samtools view -h HG04115_mat_polished_to_raw/asm2asm_aligner_outputs/HG04115_Hap2.polished.HG04115_mat_polished_to_raw.sorted.bam) > HG04115_mat_polished_to_raw/asm2asm_aligner_outputs/HG04115_Hap2.polished.HG04115_mat_polished_to_raw.sorted.pri.paf

# for polished to raw paternal
k8 /home/apps/minimap2-2.26/misc/paftools.js sam2paf -p <(samtools view -h HG04115_pat_polished_to_raw/asm2asm_aligner_outputs/HG04115_Hap1.polished.HG04115_pat_polished_to_raw.sorted.bam) > HG04115_pat_polished_to_raw/asm2asm_aligner_outputs/HG04115_Hap1.polished.HG04115_pat_polished_to_raw.sorted.pri.paf

# for maternal to paternal raw
k8 /home/apps/minimap2-2.26/misc/paftools.js sam2paf -p <(samtools view -h HG04115_raw_mat_to_pat/asm2asm_aligner_outputs/HG04115.mat.HG04115_raw_mat_to_pat.sorted.bam) > HG04115_raw_mat_to_pat/asm2asm_aligner_outputs/HG04115.mat.HG04115_raw_mat_to_pat.sorted.pri.paf
```

#### Get FP kmer Bed files and polisher output from Mira's directory

Copy merqury files
```
mkdir -p /private/groups/patenlab/masri/hprc/polishing/investigating_Y2_results/HG04115/fp_kmers/raw
cp /private/groups/hprc/polishing/batch3/hprc_polishing_QC/HG04115/hprc_polishing_QC_outputs/HG04115.merqury.tar.gz /private/groups/patenlab/masri/hprc/polishing/investigating_Y2_results/HG04115/fp_kmers/raw/HG04115.merqury.tar.gz
mkdir -p /private/groups/patenlab/masri/hprc/polishing/investigating_Y2_results/HG04115/fp_kmers/polished
cp /private/groups/hprc/polishing/batch3/hprc_polishing_QC/HG04115/hprc_polishing_QC_outputs/HG04115.polished.merqury.tar.gz /private/groups/patenlab/masri/hprc/polishing/investigating_Y2_results/HG04115/fp_kmers/polished/HG04115.polished.merqury.tar.gz
```

Extract files
```
cd /private/groups/patenlab/masri/hprc/polishing/investigating_Y2_results/HG04115/fp_kmers/polished/
tar -zvxf HG04115.polished.merqury.tar.gz 

cd /private/groups/patenlab/masri/hprc/polishing/investigating_Y2_results/HG04115/fp_kmers/raw
tar -xvzf HG04115.merqury.tar.gz 
```

Merged fp kmer tracks with counts
```
cd /private/groups/patenlab/masri/hprc/polishing/investigating_Y2_results/HG04115/fp_kmers/raw
cat HG04115.merqury.asm_only.bed | bedtools merge -o count -c 1 > ../HG04115.hap1.raw.fp_kmers.merged_count.bed
cat HG04115.merqury.altHap_only.bed | bedtools merge -o count -c 1 > ../HG04115.hap2.raw.fp_kmers.merged_count.bed

cd /private/groups/patenlab/masri/hprc/polishing/investigating_Y2_results/HG04115/fp_kmers/polished/
cat HG04115.polished.merqury.asm_only.bed | bedtools merge -o count -c 1 > ../HG04115.hap1.polished.fp_kmers.merged_count.bed
cat HG04115.polished.merqury.altHap_only.bed | bedtools merge -o count -c 1 > ../HG04115.hap2.polished.fp_kmers.merged_count.bed


# cat two haps in a single bed
# for polished
cat HG04115.hap1.polished.fp_kmers.merged_count.bed HG04115.hap2.polished.fp_kmers.merged_count.bed > HG04115.dip.polished.fp_kmers.merged_count.bed
# for raw
cat HG04115.hap1.raw.fp_kmers.merged_count.bed HG04115.hap2.raw.fp_kmers.merged_count.bed > HG04115.dip.raw.fp_kmers.merged_count.bed
```

Make a merged paf file for projecting from polished to raw assembly
```
cd /private/groups/patenlab/masri/hprc/polishing/investigating_Y2_results/HG04115/asm_alignment/asm2asm_aligner_output_jsons/HG04115_dip_polished_to_raw
cat ../HG04115_mat_polished_to_raw/asm2asm_aligner_outputs/HG04115_Hap2.polished.HG04115_mat_polished_to_raw.sorted.pri.paf ../HG04115_pat_polished_to_raw/asm2asm_aligner_outputs/HG04115_Hap1.polished.HG04115_pat_polished_to_raw.sorted.pri.paf  > HG04115_dip.polished_to_raw.paf
```

#### Project FP kmers from polished assembly to raw

```
# run docker with the projection script
docker run -u$(id -u):$(id -g) --rm -it -v/private/groups/patenlab/masri/:/private/groups/patenlab/masri/ mobinasri/flagger@sha256:5d738412b56bac5a64227569c1d6e57e7920e3d3e5724c17ab233f92279bcff6

mkdir -p /private/groups/patenlab/masri/hprc/polishing/investigating_Y2_results/HG04115/fp_kmers/projections
cd /private/groups/patenlab/masri/hprc/polishing/investigating_Y2_results/HG04115/fp_kmers/projections

# run projection
python3 /home/programs/src/project_blocks_multi_thread.py --mode asm2ref --paf ../../asm_alignment/asm2asm_aligner_output_jsons/HG04115_dip_polished_to_raw/HG04115_dip.polished_to_raw.paf --blocks ../HG04115.dip.polished.fp_kmers.merged_count.bed --outputProjectable HG04115.dip.polished.fp_kmers.merged_count.projectable.bed --outputProjection HG04115.dip.polished.fp_kmers.merged_count.projection_to_raw.bed --threads 8
```

#### Look at some numbers

```
# Number of fp kmers in the polished assembly
cat HG04115.dip.polished.fp_kmers.merged_count.bed | awk '{s+=$4}END{print s}'
121376

# Number of fp kmers in the raw assembly
cat HG04115.dip.raw.fp_kmers.merged_count.bed | awk '{s+=$4}END{print s}'
83423

# Number of fp kmers only in the polished assembly with no overlap
# with any contiguous FP block in the raw assembly by using parameter (-A)
bedtools subtract -a projections/HG04115.dip.polished.fp_kmers.merged_count.projection_to_raw.sorted.bed \
    -b HG04115.dip.raw.fp_kmers.merged_count.bed -A | \
    awk '{s+=$4}END{print s}'
61148

# Number of fp kmers only in the raw assembly with no overlap
# with any contiguous FP block in the polished assembly by using parameter (-A)
bedtools subtract -a HG04115.dip.raw.fp_kmers.merged_count.bed 
-b projections/HG04115.dip.polished.fp_kmers.merged_count.projection_to_raw.sorted.bed -A | \
    awk '{s+=$4}END{print s}'
22973

# how many FP kmers are located in FP kmer block with more than 10 kmers
bedtools subtract -a projections/HG04115.dip.polished.fp_kmers.merged_count.projection_to_raw.sorted.bed  \
    -b HG04115.dip.raw.fp_kmers.merged_count.bed -A | \
    awk '$4>10{s+=$4}END{print s}'
31978

```

So we are fixing about 22k FP kmers and adding 61k FP kmers by polishing. Nearly half of the induced FP kmers are located in the FP kmer blocks with more than 10 kmers.

### List of regions with induced FP kmers

```
# top 10 induced FP kmer blocks
bedtools subtract -a projections/HG04115.dip.polished.fp_kmers.merged_count.projection_to_raw.sorted.bed  -b HG04115.dip.raw.fp_kmers.merged_count.bed -A | sort -k4,4 -nr | head

h2tg000017l	44007725	44007744	75
h2tg000017l	44007596	44007615	75
h1tg000025l	11115931	11115967	46
h2tg000004l	54543133	54543171	39
h2tg000009l	13643	13676	33
h1tg000010l	14340	14367	32
h1tg000006l	103375874	103375936	32
h1tg000021l	81115646	81115725	29
h1tg000006l	128659088	128659124	29
h1tg000005l	91526502	91526513	24

# select 10 random FP kmer blocks with more than 10 FP kmers added after polishing
bedtools subtract -a projections/HG04115.dip.polished.fp_kmers.merged_count.projection_to_raw.sorted.bed     -b HG04115.dip.raw.fp_kmers.merged_count.bed -A | awk '$4>10' | shuf -n10

h1tg000010l	39501395	39501429	15
h1tg000023l	37458874	37458906	13
h1tg000010l	58275971	58276003	13
h1tg000007l	65505462	65505492	11
h1tg000001l	58898253	58898283	11
h2tg000004l	20621992	20622023	11
h2tg000029l	8755529	8755560	12
h1tg000013l	102495348	102495378	11
h2tg000002l	83450466	83450497	12
h1tg000012l	40544128	40544158	11
```

### List of regions with fixed FP kmers
```
# top 10 fixed FP kmer blocks
bedtools subtract -b projections/HG04115.dip.polished.fp_kmers.merged_count.projection_to_raw.sorted.bed  -a HG04115.dip.raw.fp_kmers.merged_count.bed -A | sort -k4,4 -nr | head
h2tg000029l	1322059	1322280	176
h1tg000028l	593951	594085	86
h1tg000005l	189129735	189129843	60
h1tg000016l	67184872	67184959	59
h1tg000010l	30912006	30912079	52
h2tg000015l	130537683	130537764	51
h1tg000001l	87103240	87103308	45
h1tg000016l	67186145	67186217	42
h2tg000015l	32495310	32495369	39
h1tg000016l	67186020	67186087	39

# select 10 random FP kmer blocks with more than 10 FP kmers removed after polishing
bedtools subtract -b projections/HG04115.dip.polished.fp_kmers.merged_count.projection_to_raw.sorted.bed  -a HG04115.dip.raw.fp_kmers.merged_count.bed -A | awk '$4>10' | shuf -n10
h1tg000006l	64301616	64301657	21
h2tg000006l	25484793	25484830	17
h2tg000001l	4485075	4485108	13
h2tg000014l	152627	152666	19
h2tg000015l	14798062	14798100	18
h1tg000018l	163142884	163142924	20
h1tg000007l	64781285	64781326	20
h1tg000020l	17544148	17544187	19
h1tg000037l	3132060	3132098	18
h1tg000008l	30551358	30551398	20

```

### Project blocks from hap1 to hap2 and reverse

#### Top 10 fixed FP kmer blocks
I made two bed files pointing to the top 10 fixed blocks; one for hap1 and the other for hap2. I mapped raw-hap2 to raw-hap1 above so for projection from hap1 to hap2 we need to use `--mode ref2asm` and for projection from hap2 to hap1 `--mode asm2ref`.
```
cd /private/groups/patenlab/masri/hprc/polishing/investigating_Y2_results/HG04115/analysis/top_10_fixed_FP_kmer_block

cat top_10_fixed_FP_kmer_blocks.hap1.bed
h1tg000001l	87103240	87103308	45
h1tg000005l	189129735	189129843	60
h1tg000010l	30912006	30912079	52
h1tg000016l	67184872	67184959	59
h1tg000016l	67186020	67186087	39
h1tg000016l	67186145	67186217	42
h1tg000028l	593951	594085	86

cat top_10_fixed_FP_kmer_blocks.hap2.bed
h2tg000015l	32495310	32495369	39
h2tg000015l	130537683	130537764	51
h2tg000029l	1322059	1322280	176
```

#### Project top 10 fixed FP kmer blocks
```
mkdir -p /private/groups/patenlab/masri/hprc/polishing/investigating_Y2_results/HG04115/analysis/top_10_fixed_FP_kmer_block/projections

cd /private/groups/patenlab/masri/hprc/polishing/investigating_Y2_results/HG04115/analysis/top_10_fixed_FP_kmer_block/projections


# hap2 with asm2ref
python3 /home/programs/src/project_blocks_multi_thread.py \
    --mode asm2ref \
    --paf ../../../asm_alignment/asm2asm_aligner_output_jsons/HG04115_raw_mat_to_pat/asm2asm_aligner_outputs/HG04115.mat.HG04115_raw_mat_to_pat.sorted.pri.paf \
    --blocks ../top_10_fixed_FP_kmer_blocks.hap2.bed \
    --outputProjectable top_10_fixed_FP_kmer_blocks.hap2.projectable.bed \
    --outputProjection top_10_fixed_FP_kmer_blocks.hap2.projection.bed \
    --threads 8

# hap1 with ref2asm
python3 /home/programs/src/project_blocks_multi_thread.py \
    --mode ref2asm \
    --paf ../../../asm_alignment/asm2asm_aligner_output_jsons/HG04115_raw_mat_to_pat/asm2asm_aligner_outputs/HG04115.mat.HG04115_raw_mat_to_pat.sorted.pri.paf \
    --blocks ../top_10_fixed_FP_kmer_blocks.hap1.bed \
    --outputProjectable top_10_fixed_FP_kmer_blocks.hap1.projectable.bed \
    --outputProjection top_10_fixed_FP_kmer_blocks.hap1.projection.bed \
    --threads 8
```

Take a look at projectables and projections:
```
# for hap1
cat top_10_fixed_FP_kmer_blocks.hap1.projectable.bed
h1tg000005l	189129735	189129755	60
h1tg000005l	189129823	189129843	60
h1tg000010l	30912006	30912079	52
h1tg000016l	67186145	67186217	42
h1tg000016l	67186020	67186087	39
h1tg000016l	67184872	67184959	59

cat top_10_fixed_FP_kmer_blocks.hap1.projection.bed
h2tg000003l	187651139	187651159	60
h2tg000003l	187651227	187651247	60
h2tg000025l	49252606	49252676	52
h2tg000018l	45192569	45192643	42
h2tg000018l	45192701	45192766	39
h2tg000018l	45193825	45193910	59

# for hap2
cat top_10_fixed_FP_kmer_blocks.hap2.projectable.bed
h2tg000015l	32495310	32495369	39
h2tg000029l	1322260	1322280	176

cat top_10_fixed_FP_kmer_blocks.hap2.projection.bed
h1tg000007l	114453511	114453569	39
h1tg000013l	135811879	135811899	176

```

#### 10 random fixed FP kmer blocks

```
cd /private/groups/patenlab/masri/hprc/polishing/investigating_Y2_results/HG04115/analysis/random_10_fixed_FP_kmer_block

cat random_10_fixed_FP_kmer_blocks.hap1.bed
h1tg000006l	64301616	64301657	21
h1tg000007l	64781285	64781326	20
h1tg000008l	30551358	30551398	20
h1tg000018l	163142884	163142924	20
h1tg000020l	17544148	17544187	19
h1tg000037l	3132060	3132098	18

cat random_10_fixed_FP_kmer_blocks.hap2.bed
h2tg000001l	4485075	4485108	13
h2tg000006l	25484793	25484830	17
h2tg000014l	152627	152666	19
h2tg000015l	14798062	14798100	18
```

#### Project 10 random fixed FP kmer blocks
```
mkdir -p /private/groups/patenlab/masri/hprc/polishing/investigating_Y2_results/HG04115/analysis/random_10_fixed_FP_kmer_block/projections

cd /private/groups/patenlab/masri/hprc/polishing/investigating_Y2_results/HG04115/analysis/random_10_fixed_FP_kmer_block/projections

# project hap2 with asm2ref
python3 /home/programs/src/project_blocks_multi_thread.py \
    --mode asm2ref \
    --paf ../../../asm_alignment/asm2asm_aligner_output_jsons/HG04115_raw_mat_to_pat/asm2asm_aligner_outputs/HG04115.mat.HG04115_raw_mat_to_pat.sorted.pri.paf \
    --blocks ../random_10_fixed_FP_kmer_blocks.hap2.bed \
    --outputProjectable random_10_fixed_FP_kmer_blocks.hap2.projectable.bed \
    --outputProjection random_10_fixed_FP_kmer_blocks.hap2.projection.bed \
    --threads 8

# project hap1 with ref2asm
python3 /home/programs/src/project_blocks_multi_thread.py \
    --mode ref2asm \
    --paf ../../../asm_alignment/asm2asm_aligner_output_jsons/HG04115_raw_mat_to_pat/asm2asm_aligner_outputs/HG04115.mat.HG04115_raw_mat_to_pat.sorted.pri.paf \
    --blocks ../random_10_fixed_FP_kmer_blocks.hap1.bed \
    --outputProjectable random_10_fixed_FP_kmer_blocks.hap1.projectable.bed \
    --outputProjection random_10_fixed_FP_kmer_blocks.hap1.projection.bed \
    --threads 8
```

Take a look at projectables and projections:
```
# for hap1
cat random_10_fixed_FP_kmer_blocks.hap1.projectable.bed
h1tg000006l	64301616	64301657	21
h1tg000007l	64781285	64781326	20
h1tg000008l	30551378	30551398	20
h1tg000018l	163142884	163142924	20
h1tg000020l	17544148	17544187	19
h1tg000037l	3132060	3132098	18
h1tg000037l	3132060	3132098	18
h1tg000037l	3132060	3132098	18

cat random_10_fixed_FP_kmer_blocks.hap1.projection.bed
h2tg000007l	176613905	176613946	21
h2tg000015l	82595412	82595453	20
h2tg000001l	30587378	30587398	20
h2tg000008l	161802339	161802371	20
h2tg000011l	18515617	18515667	19
h2tg000056l	28131	28168	18
h2tg000085l	6414	6451	18
h2tg000040l	70561	70598	18

# for hap2
cat random_10_fixed_FP_kmer_blocks.hap2.projectable.bed 
h2tg000014l	152627	152666	19
h2tg000015l	14798062	14798100	18
h2tg000001l	4485075	4485108	13

cat random_10_fixed_FP_kmer_blocks.hap2.projection.bed
h1tg000001l	92915720	92915769	19
h1tg000007l	132161530	132161568	18
h1tg000008l	4447992	4448025	13
```

#### Top 10 induced FP kmer blocks
```
mkdir -p /private/groups/patenlab/masri/hprc/polishing/investigating_Y2_results/HG04115/analysis/top_10_induced_FP_kmer_block

cd /private/groups/patenlab/masri/hprc/polishing/investigating_Y2_results/HG04115/analysis/top_10_induced_FP_kmer_block

cat top_10_induced_FP_kmer_blocks.hap1.bed
h1tg000005l	91526502	91526513	24
h1tg000006l	103375874	103375936	32
h1tg000006l	128659088	128659124	29
h1tg000010l	14340	14367	32
h1tg000021l	81115646	81115725	29
h1tg000025l	11115931	11115967	46

cat top_10_induced_FP_kmer_blocks.hap2.bed
h2tg000004l	54543133	54543171	39
h2tg000009l	13643	13676	33
h2tg000017l	44007596	44007615	75
h2tg000017l	44007725	44007744	75
```

#### Project top 10 induced FP kmer blocks
```
mkdir -p /private/groups/patenlab/masri/hprc/polishing/investigating_Y2_results/HG04115/analysis/top_10_induced_FP_kmer_block/projections

cd /private/groups/patenlab/masri/hprc/polishing/investigating_Y2_results/HG04115/analysis/top_10_induced_FP_kmer_block/projections

python3 /home/programs/src/project_blocks_multi_thread.py \
    --mode asm2ref \
    --paf ../../../asm_alignment/asm2asm_aligner_output_jsons/HG04115_raw_mat_to_pat/asm2asm_aligner_outputs/HG04115.mat.HG04115_raw_mat_to_pat.sorted.pri.paf \
    --blocks ../top_10_induced_FP_kmer_blocks.hap2.bed \
    --outputProjectable top_10_induced_FP_kmer_blocks.hap2.projectable.bed \
    --outputProjection top_10_induced_FP_kmer_blocks.hap2.projection.bed \
    --threads 8

python3 /home/programs/src/project_blocks_multi_thread.py \
    --mode ref2asm \
    --paf ../../../asm_alignment/asm2asm_aligner_output_jsons/HG04115_raw_mat_to_pat/asm2asm_aligner_outputs/HG04115.mat.HG04115_raw_mat_to_pat.sorted.pri.paf \
    --blocks ../top_10_induced_FP_kmer_blocks.hap1.bed \
    --outputProjectable top_10_induced_FP_kmer_blocks.hap1.projectable.bed \
    --outputProjection top_10_induced_FP_kmer_blocks.hap1.projection.bed \
    --threads 8

```

Take a look at projectables and projections:
```
# for hap1

cat top_10_induced_FP_kmer_blocks.hap1.projectable.bed
h1tg000005l	91526502	91526513	24
h1tg000006l	103375874	103375936	32
h1tg000006l	128659088	128659124	29
h1tg000010l	14340	14367	32
h1tg000021l	81115646	81115725	29
h1tg000025l	11115931	11115967	46

cat top_10_induced_FP_kmer_blocks.hap1.projection.bed
h2tg000003l	90149103	90149114	24
h2tg000007l	138488216	138488278	32
h2tg000007l	113193927	113193981	29
h2tg000025l	80009747	80009774	32
h2tg000005l	80494787	80494876	29
h2tg000022l	122731862	122731941	46

# for hap2
cat top_10_induced_FP_kmer_blocks.hap2.projectable.bed
h2tg000017l	44007596	44007615	75
h2tg000017l	44007725	44007744	75
h2tg000009l	13643	13676	33

cat top_10_induced_FP_kmer_blocks.hap2.projection.bed
h1tg000004l	44101416	44101435	75
h1tg000004l	44101545	44101564	75
h1tg000028l	14193	14226	33
```

#### 10 random induced FP kmer blocks

```
mkdir -p /private/groups/patenlab/masri/hprc/polishing/investigating_Y2_results/HG04115/analysis/random_10_induced_FP_kmer_block

cd /private/groups/patenlab/masri/hprc/polishing/investigating_Y2_results/HG04115/analysis/random_10_induced_FP_kmer_block

cat random_10_induced_FP_kmer_blocks.hap1.bed
h1tg000001l	58898253	58898283	11
h1tg000007l	65505462	65505492	11
h1tg000010l	39501395	39501429	15
h1tg000010l	58275971	58276003	13
h1tg000012l	40544128	40544158	11
h1tg000013l	102495348	102495378	11
h1tg000023l	37458874	37458906	13

cat random_10_induced_FP_kmer_blocks.hap2.bed
h2tg000002l	83450466	83450497	12
h2tg000004l	20621992	20622023	11
h2tg000029l	8755529	8755560	12
```

#### Project 10 random induced FP kmer blocks
```
mkdir -p /private/groups/patenlab/masri/hprc/polishing/investigating_Y2_results/HG04115/analysis/random_10_induced_FP_kmer_block/projections

cd /private/groups/patenlab/masri/hprc/polishing/investigating_Y2_results/HG04115/analysis/random_10_induced_FP_kmer_block/projections

python3 /home/programs/src/project_blocks_multi_thread.py \
    --mode asm2ref \
    --paf ../../../asm_alignment/asm2asm_aligner_output_jsons/HG04115_raw_mat_to_pat/asm2asm_aligner_outputs/HG04115.mat.HG04115_raw_mat_to_pat.sorted.pri.paf \
    --blocks ../random_10_induced_FP_kmer_blocks.hap2.bed \
    --outputProjectable random_10_induced_FP_kmer_blocks.hap2.projectable.bed \
    --outputProjection random_10_induced_FP_kmer_blocks.hap2.projection.bed \
    --threads 8

python3 /home/programs/src/project_blocks_multi_thread.py \
    --mode ref2asm \
    --paf ../../../asm_alignment/asm2asm_aligner_output_jsons/HG04115_raw_mat_to_pat/asm2asm_aligner_outputs/HG04115.mat.HG04115_raw_mat_to_pat.sorted.pri.paf \
    --blocks ../random_10_induced_FP_kmer_blocks.hap1.bed \
    --outputProjectable random_10_induced_FP_kmer_blocks.hap1.projectable.bed \
    --outputProjection random_10_induced_FP_kmer_blocks.hap1.projection.bed \
    --threads 8
```

Take a look at projectables and projections:
```
# for hap1
cat random_10_induced_FP_kmer_blocks.hap1.projectable.bed
h1tg000001l	58898253	58898273	11
h1tg000007l	65505462	65505492	11
h1tg000010l	39501395	39501429	15
h1tg000010l	58275971	58276003	13
h1tg000012l	40544128	40544158	11
h1tg000013l	102495348	102495378	11
h1tg000023l	37458874	37458906	13

cat random_10_induced_FP_kmer_blocks.hap1.projection.bed
h2tg000024l	58859726	58859746	11
h2tg000015l	81872185	81872217	11
h2tg000025l	40636928	40636962	15
h2tg000025l	21754013	21754057	13
h2tg000016l	40541562	40541592	11
h2tg000029l	34668000	34668036	11
h2tg000020l	36393706	36393738	13

# for hap2
cat random_10_induced_FP_kmer_blocks.hap2.projectable.bed
h2tg000004l	20621992	20622023	11
h2tg000029l	8755529	8755560	12
h2tg000002l	83450466	83450497	12

cat random_10_induced_FP_kmer_blocks.hap2.projection.bed
h1tg000002l	140513931	140513962	11
h1tg000013l	128410074	128410105	12
h1tg000017l	21162439	21162470	12
```

Comment 4 : 02/09/2024

#### Project FP kmer block to chm13v2.0 coordinates

Make the paf files
```
# get the docker image with paftools
docker run -u$(id -u):$(id -g) --rm -it -v/private/groups/patenlab/masri/:/private/groups/patenlab/masri/ mobinasri/long_read_aligner:v0.3.3

cd /private/groups/patenlab/masri/hprc/polishing/investigating_Y2_results/HG04115/asm_alignment/asm2asm_aligner_output_jsons

# convert bam to paf to maternal alignment
k8 /home/apps/minimap2-2.26/misc/paftools.js sam2paf -p <(samtools view -h HG04115_raw_mat_to_chm13v2/asm2asm_aligner_outputs/HG04115.mat.HG04115_raw_mat_to_chm13v2.sorted.bam) > HG04115_raw_mat_to_chm13v2/asm2asm_aligner_outputs/HG04115.mat.HG04115_raw_mat_to_chm13v2.sorted.pri.paf

# convert bam to paf to paternal alignment
k8 /home/apps/minimap2-2.26/misc/paftools.js sam2paf -p <(samtools view -h HG04115_raw_pat_to_chm13v2/asm2asm_aligner_outputs/HG04115.pat.HG04115_raw_pat_to_chm13v2.sorted.bam) > HG04115_raw_pat_to_chm13v2/asm2asm_aligner_outputs/HG04115.pat.HG04115_raw_pat_to_chm13v2.sorted.pri.paf
```

Make bed files with all induced and fixed FP kmer blocks in the coordinates of the raw assembly
```
# get the docker image with projection script
docker run -u$(id -u):$(id -g) --rm -it -v/private/groups/patenlab/masri/:/private/groups/patenlab/masri/ mobinasri/flagger@sha256:5d738412b56bac5a64227569c1d6e57e7920e3d3e5724c17ab233f92279bcff6

cd /private/groups/patenlab/masri/hprc/polishing/investigating_Y2_results/HG04115/fp_kmers/analysis

# induced FP kmer block for hap1
bedtools subtract -a ../projections/HG04115.dip.polished.fp_kmers.merged_count.projection_to_raw.sorted.bed  -b ../HG04115.dip.raw.fp_kmers.merged_count.bed -A | grep "h1tg" | bedtools sort -i - > induced_fp_kmer_blocks.hap1.bed

# induced FP kmer block for hap2
bedtools subtract -a ../projections/HG04115.dip.polished.fp_kmers.merged_count.projection_to_raw.sorted.bed  -b ../HG04115.dip.raw.fp_kmers.merged_count.bed -A | grep "h2tg" | bedtools sort -i - > induced_fp_kmer_blocks.hap2.bed

# fixed FP kmer block for hap1
bedtools subtract -b ../projections/HG04115.dip.polished.fp_kmers.merged_count.projection_to_raw.sorted.bed  -a ../HG04115.dip.raw.fp_kmers.merged_count.bed -A | grep "h1tg" | bedtools sort -i - > fixed_fp_kmer_blocks.hap1.bed

# fixed FP kmer block for hap2
bedtools subtract -b ../projections/HG04115.dip.polished.fp_kmers.merged_count.projection_to_raw.sorted.bed  -a ../HG04115.dip.raw.fp_kmers.merged_count.bed -A | grep "h2tg" | bedtools sort -i - > fixed_fp_kmer_blocks.hap2.bed
```

Project BED files to chm13 coordinates
```

# fixed hap1
python3 /home/programs/src/project_blocks_multi_thread.py \
    --mode asm2ref \
    --blocks fixed_fp_kmer_blocks.hap1.bed \
    --paf ../../asm_alignment/asm2asm_aligner_output_jsons/HG04115_raw_pat_to_chm13v2/asm2asm_aligner_outputs/HG04115.pat.HG04115_raw_pat_to_chm13v2.sorted.pri.paf \
    --outputProjectable fixed_fp_kmer_blocks.hap1.to_chm13v2.projectable.bed \
    --outputProjection fixed_fp_kmer_blocks.hap1.to_chm13v2.projection.bed \
    --threads 8

bedtools sort -i fixed_fp_kmer_blocks.hap1.to_chm13v2.projection.bed > fixed_fp_kmer_blocks.hap1.to_chm13v2.projection.sorted.bed


# fixed hap2

python3 /home/programs/src/project_blocks_multi_thread.py \
    --mode asm2ref \
    --blocks fixed_fp_kmer_blocks.hap2.bed \
    --paf ../../asm_alignment/asm2asm_aligner_output_jsons/HG04115_raw_mat_to_chm13v2/asm2asm_aligner_outputs/HG04115.mat.HG04115_raw_mat_to_chm13v2.sorted.pri.paf \
    --outputProjectable fixed_fp_kmer_blocks.hap2.to_chm13v2.projectable.bed \
    --outputProjection fixed_fp_kmer_blocks.hap2.to_chm13v2.projection.bed \
    --threads 8

bedtools sort -i fixed_fp_kmer_blocks.hap2.to_chm13v2.projection.bed > fixed_fp_kmer_blocks.hap2.to_chm13v2.projection.sorted.bed


# induced hap1
python3 /home/programs/src/project_blocks_multi_thread.py \
    --mode asm2ref \
    --blocks induced_fp_kmer_blocks.hap1.bed \
    --paf ../../asm_alignment/asm2asm_aligner_output_jsons/HG04115_raw_pat_to_chm13v2/asm2asm_aligner_outputs/HG04115.pat.HG04115_raw_pat_to_chm13v2.sorted.pri.paf \
    --outputProjectable induced_fp_kmer_blocks.hap1.to_chm13v2.projectable.bed \
    --outputProjection induced_fp_kmer_blocks.hap1.to_chm13v2.projection.bed \
    --threads 8

bedtools sort -i induced_fp_kmer_blocks.hap1.to_chm13v2.projection.bed > induced_fp_kmer_blocks.hap1.to_chm13v2.projection.sorted.bed

# induced hap2
python3 /home/programs/src/project_blocks_multi_thread.py \
    --mode asm2ref \
    --blocks induced_fp_kmer_blocks.hap2.bed \
    --paf ../../asm_alignment/asm2asm_aligner_output_jsons/HG04115_raw_mat_to_chm13v2/asm2asm_aligner_outputs/HG04115.mat.HG04115_raw_mat_to_chm13v2.sorted.pri.paf \
    --outputProjectable induced_fp_kmer_blocks.hap2.to_chm13v2.projectable.bed \
    --outputProjection induced_fp_kmer_blocks.hap2.to_chm13v2.projection.bed \
    --threads 8

bedtools sort -i induced_fp_kmer_blocks.hap2.to_chm13v2.projection.bed > induced_fp_kmer_blocks.hap2.to_chm13v2.projection.sorted.bed

```

Make cov files for each projection bed file:
```
TYPE="fixed"
HAP="hap1"

bedtools subtract -a chm13v2.0.bed -b ../${TYPE}_fp_kmer_blocks.${HAP}.to_chm13v2.projection.sorted.bed | awk '{print $1"\t"$2"\t"$3"\t"0}' > ${TYPE}_fp_kmer_blocks.${HAP}.to_chm13v2.projection.sorted.only_zero.bed

# make a bed file covering the whole genome
cat ../${TYPE}_fp_kmer_blocks.${HAP}.to_chm13v2.projection.sorted.bed ${TYPE}_fp_kmer_blocks.${HAP}.to_chm13v2.projection.sorted.only_zero.bed | bedtools sort -i - > ${TYPE}_fp_kmer_blocks.${HAP}.to_chm13v2.projection.sorted.complete.bed

# make a cov file
awk 'BEGIN{ctg=""}FNR==NR{ctg_len[$1]=$2;next} FNR!=NR{if(ctg != $1){ctg=$1; print ">"ctg" "ctg_len[ctg]}; print $2+1"\t"$3"\t"$4}' chm13v2.0.fa.fai ${TYPE}_fp_kmer_blocks.${HAP}.to_chm13v2.projection.sorted.complete.bed > ${TYPE}_fp_kmer_blocks.${HAP}.to_chm13v2.projection.sorted.complete.cov

# convert cov to wig
cov2wig -i ${TYPE}_fp_kmer_blocks.${HAP}.to_chm13v2.projection.sorted.complete.cov -f chm13v2.0.fa.fai -s 100 -o ${TYPE}_fp_kmer_blocks.${HAP}.to_chm13v2.projection.sorted.complete.wig -n ${TYPE}_${HAP}

```
