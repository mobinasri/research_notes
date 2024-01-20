#!/bin/bash
#SBATCH --job-name=gbz_to_maf
#SBATCH --partition=long
#SBATCH --mail-user=masri@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=256gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --output=%x.%j.log
#SBATCH --time=72:00:00

pwd; hostname; date

set -e

USER_NAME=$(jq -r ".user_name" $json)
USER_UID="$(id -u ${USER_NAME})"
USER_GID="$(id -g ${USER_NAME})"


VG_BIN=$1
VG2MAF_BIN=$2
INPUT_GBZ=$3

FILENAME=$(basename ${INPUT_GBZ})
DIRNAME=$(dirname ${INPUT_GBZ})
PREFIX=${FILENAME%%.gbz}
ROOT_DIR="/private/groups/patenlab/masri"
REF_SAMPLE="GRCh38"

echo "[$(date)] Convert ${PREFIX}.gbz ${PREFIX}.vg"
# create vg file
${VG_BIN} convert ${INPUT_GBZ} -p > ${DIRNAME}/${PREFIX}.vg

echo "[$(date)] Convert ${PREFIX}.gbz to ${PREFIX}.gam"
# create gam file
${VG_BIN} paths -x ${INPUT_GBZ} -H --extract-gam > ${DIRNAME}/${PREFIX}.gam

echo "[$(date)] Sort and index ${PREFIX}.gam"
# sort and index gam
${VG_BIN} gamsort \
	${DIRNAME}/${PREFIX}.gam \
  -i ${DIRNAME}/${PREFIX}.sort.gam.gai \
	--threads 8 \
	-p > ${DIRNAME}/${PREFIX}.sort.gam

echo "[$(date)] Create ${PREFIX}.dist"
# make the distance index (if you have one for the gbz, you can just use that)
${VG_BIN} index \
	${DIRNAME}/${PREFIX}.vg \
	-j ${DIRNAME}/${PREFIX}.dist \
	--threads 8

echo "[$(date)] Create ${PREFIX}.maf"
${VG2MAF_BIN} \
	${DIRNAME}/${PREFIX}.vg \
	-d ${DIRNAME}/${PREFIX}.dist \
	-r ${REF_SAMPLE} \
	-g ${DIRNAME}/${PREFIX}.sort.gam \
	-p --threads 8 > ${DIRNAME}/${PREFIX}.maf
