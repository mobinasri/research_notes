#!/bin/bash
#SBATCH --job-name=cram2fq
#SBATCH --partition=long
#SBATCH --mail-user=masri@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=16gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --output=%x.%j.log
#SBATCH --time=72:00:00

pwd; hostname; date

set -x
set -e


json=$1

SAMTOOLS_BIN="/private/groups/patenlab/masri/apps/samtools-1.19.2/samtools"

CRAM_URL=$(jq -r ".cram_url" $json)
REF_FASTA=$(jq -r ".ref_fasta" $json)
OUTPUT_PREFIX=$(jq -r ".output_prefix" $json)
OUTPUT_DIR=$(jq -r ".output_dir" $json)
OPTIONS=$(jq -r ".other_options" $json)


mkdir -p ${OUTPUT_DIR}
chmod 777 ${OUTPUT_DIR}
cd ${OUTPUT_DIR}

# download haplotype assemblies
wget ${CRAM_URL}

CRAM_FILENAME=$(basename ${CRAM_URL})

${SAMTOOLS_BIN} fastq -@8 --reference ${REF_FASTA} -1 ${OUTPUT_PREFIX}.R1.fastq.gz -2 ${OUTPUT_PREFIX}.R2.fastq.gz ${CRAM_FILENAME}

rm -rf ${CRAM_FILENAME}

wait
date
