#!/bin/bash
#SBATCH --job-name=png-aware-dv
#SBATCH --partition=long
#SBATCH --mail-user=masri@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=260gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --output=%x.%j.log
#SBATCH --time=72:00:00

pwd; hostname; date

set -x 
set -e

json=$1

BAM=$(jq -r ".bam" $json)
PANGENOME=$(jq -r ".pangenome" $json)
MODEL_CKPT=$(jq -r ".model_ckpt" $json)
THREADS=$(jq -r ".threads" $json)
BIN_VERSION=$(jq -r ".bin_version" $json)
REF_FA=$(jq -r ".reference_fasta" $json)
OUTPUT_PREFIX=$(jq -r ".output_prefix" $json)
PANGENOME_HEIGHT=$(jq -r ".pangenome_height" $json)
OUTPUT_DIR=$(jq -r ".output_dir" $json)
MOUNT_DIR=$(jq -r ".mount_dir" $json)

mkdir -p ${OUTPUT_DIR}/intermediate_results_dir
chmod 777 ${OUTPUT_DIR}
cd ${OUTPUT_DIR}


docker run \
--rm \
-v ${MOUNT_DIR}:${MOUNT_DIR} \
--shm-size 15gb \
google/deepvariant:"${BIN_VERSION}" \
    /opt/deepvariant/bin/run_pangenome_aware_deepvariant \
    --customized_model ${MODEL_CKPT} \
    --model_type WGS \
    --ref ${REF_FA} \
    --reads ${BAM} \
    --gbz_shared_memory_size_gb 15 \
    --make_examples_extra_args "pileup_image_height_pangenome=${PANGENOME_HEIGHT}" \
    --pangenome ${PANGENOME} \
    --output_vcf ${OUTPUT_DIR}/${OUTPUT_PREFIX}.vcf.gz \
    --output_gvcf ${OUTPUT_DIR}/${OUTPUT_PREFIX}.g.vcf.gz \
    --num_shards ${THREADS} \
    --intermediate_results_dir ${OUTPUT_DIR}/intermediate_results_dir
