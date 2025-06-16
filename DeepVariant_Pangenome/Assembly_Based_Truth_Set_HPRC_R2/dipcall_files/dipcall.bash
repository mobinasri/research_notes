#!/bin/bash
#SBATCH --job-name=dipcall
#SBATCH --partition=long
#SBATCH --mail-user=masri@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=64gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --output=%x.%j.log
#SBATCH --time=72:00:00

pwd; hostname; date

set -x 
set -e


json=$1

DIPCALL_KIT_DIR=$(jq -r ".dipcall_kit_dir" $json)
THREADS=$(jq -r ".threads" $json)
PAT_FA_URL=$(jq -r ".paternal_fasta_url" $json)
MAT_FA_URL=$(jq -r ".maternal_fasta_url" $json)
REF_FA=$(jq -r ".reference_fasta" $json)
PAR_X_BED=$(jq -r ".PAR_X_bed" $json)
OUTPUT_PREFIX=$(jq -r ".output_prefix" $json)
OUTPUT_DIR=$(jq -r ".output_dir" $json)
OPTIONS=$(jq -r ".other_options" $json)


mkdir -p ${OUTPUT_DIR}
chmod 777 ${OUTPUT_DIR}
cd ${OUTPUT_DIR}

# download haplotype assemblies
wget ${PAT_FA_URL}
wget ${MAT_FA_URL}

PAT_FA=$(basename ${PAT_FA_URL})
MAT_FA=$(basename ${MAT_FA_URL})

if [ -z ${PAR_X_BED} ]; then
	# for female
	${DIPCALL_KIT_DIR}/run-dipcall ${OPTIONS} -t ${THREADS} ${OUTPUT_PREFIX} ${REF_FA} ${PAT_FA} ${MAT_FA} > prefix.mak
else
	# for male
	${DIPCALL_KIT_DIR}/run-dipcall ${OPTIONS} -t ${THREADS} -x ${PAR_X_BED} ${OUTPUT_PREFIX} ${REF_FA} ${PAT_FA} ${MAT_FA} > prefix.mak
fi

make -j2 -f prefix.mak

rm -rf ${PAT_FA} ${MAT_FA}

wait
date
