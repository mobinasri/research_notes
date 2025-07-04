#!/bin/bash
#SBATCH --job-name=happy
#SBATCH --partition=long
#SBATCH --mail-user=masri@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=16gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --output=%x.%j.log
#SBATCH --time=72:00:00

pwd; hostname; date

set -x 
set -e

json=$1

TRUTH_VCF=$(jq -r ".truth_vcf" $json)
TRUTH_CONF_BED=$(jq -r ".truth_conf_bed" $json)
QUERY_VCF=$(jq -r ".query_vcf" $json)
REF_FA=$(jq -r ".reference_fasta" $json)
OUTPUT_PREFIX=$(jq -r ".output_prefix" $json)
OUTPUT_DIR=$(jq -r ".output_dir" $json)
MOUNT_DIR=$(jq -r ".mount_dir" $json)
# male or female
GENDER=$(jq -r ".gender" $json)
THREADS=$(jq -r ".threads" $json)

mkdir -p ${OUTPUT_DIR}
chmod 777 ${OUTPUT_DIR}
cd ${OUTPUT_DIR}

HAPPY_VERSION="v0.3.12"

# male for HG002/HG003 and female for HG001 GENDER="male"
docker run -i -v${MOUNT_DIR}:${MOUNT_DIR} jmcdani20/hap.py:${HAPPY_VERSION} /opt/hap.py/bin/hap.py \
	${TRUTH_VCF} \
	${QUERY_VCF} \
	-f ${TRUTH_CONF_BED} \
	-r ${REF_FA} \
	-o ${OUTPUT_DIR}/${OUTPUT_PREFIX} \
	--gender ${GENDER} \
	--engine=vcfeval \
	--pass-only \
	--threads=${THREADS}
