#!/bin/bash
#SBATCH --job-name=map_vg_short_read
#SBATCH --partition=high_priority
#SBATCH --mail-user=masri@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=256gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --output=%x.%j.log
#SBATCH --time=72:00:00

pwd; hostname; date

set -x 
set -e

json=$1

THREADS=$(jq -r ".threads" $json)
FASTQ_URL_R1=$(jq -r ".fastq_url_r1" $json)
FASTQ_URL_R2=$(jq -r ".fastq_url_r2" $json)
GBZ=$(jq -r ".gbz" $json)
SAMPLE_NAME=$(jq -r ".sample_name" $json)
OUTPUT_SUFFIX=$(jq -r ".output_suffix" $json)
OUTPUT_DIR=$(jq -r ".output_dir" $json)
VG_BIN=$(jq -r ".vg_bin" $json)
OPTIONS=$(jq -r ".options" $json)

mkdir -p ${OUTPUT_DIR}
chmod 777 ${OUTPUT_DIR}
cd ${OUTPUT_DIR}

# download files
wget -q ${FASTQ_URL_R1}
wget -q ${FASTQ_URL_R2}

# Extract the filename from the URL
FASTQ_R1=${PWD}/$(basename ${FASTQ_URL_R1})
FASTQ_R2=${PWD}/$(basename ${FASTQ_URL_R2})

cat > ${SAMPLE_NAME}.fq.paths <<- EOM
${FASTQ_R1}
${FASTQ_R2}
EOM


${VG_BIN} paths \
  -x ${GBZ} \
  -L -Q GRCh38 > GRCh38.path_list.txt

echo "mapping reads ..."

time ${VG_BIN} giraffe --progress \
  --read-group "ID:1 LB:lib1 SM:${SAMPLE_NAME} PL:illumina PU:unit1" \
  --sample ${SAMPLE_NAME} \
  -o BAM --ref-paths GRCh38.path_list.txt \
  -P -L 3000 ${OPTIONS} \
  -f ${FASTQ_R1} \
  -f ${FASTQ_R2} \
  -Z ${GBZ} \
  -t ${THREADS} > ${SAMPLE_NAME}.unsorted.bam

echo "sorting bam ..."
BAM=${SAMPLE_NAME}.${OUTPUT_SUFFIX}.bam
time samtools view -h ${SAMPLE_NAME}.unsorted.bam | sed -e "s/GRCh38#0#//g" | samtools sort --threads 8 -m 2G -O BAM > ${BAM}

echo "removing unsorted bam and fastq files..."
rm -rf ${SAMPLE_NAME}.unsorted.bam
rm -rf ${FASTQ_R1} ${FASTQ_R2}


echo "indexing bam ..."
# Index the BAM.
samtools index -@8 ${BAM}
