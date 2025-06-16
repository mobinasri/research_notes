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
FASTQ_R1=$(jq -r ".fastq_r1" $json)
FASTQ_R2=$(jq -r ".fastq_r2" $json)
GBZ=$(jq -r ".gbz" $json)
HAPL=$(jq -r ".hapl" $json)
SAMPLE_NAME=$(jq -r ".sample_name" $json)
OUTPUT_SUFFIX=$(jq -r ".output_suffix" $json)
OUTPUT_DIR=$(jq -r ".output_dir" $json)
KMC_BIN=$(jq -r ".kmc_bin" $json)
VG_BIN=$(jq -r ".vg_bin" $json)

mkdir -p ${OUTPUT_DIR}
chmod 777 ${OUTPUT_DIR}
cd ${OUTPUT_DIR}

cat > ${SAMPLE_NAME}.fq.paths <<- EOM
${FASTQ_R1}
${FASTQ_R2}
EOM

echo "getting kmers ..."

TMPDIR=$(mktemp -d)
time ${KMC_BIN} -k29 -m128 -okff -t${THREADS} @${SAMPLE_NAME}.fq.paths ${SAMPLE_NAME}.fq $TMPDIR


${VG_BIN} paths \
  -x ${GBZ} \
  -L -Q GRCh38 > GRCh38.path_list.txt

echo "mapping reads ..."
time ${VG_BIN} giraffe --progress \
  --read-group "ID:1 LB:lib1 SM:${SAMPLE_NAME} PL:illumina PU:unit1" \
  --sample ${SAMPLE_NAME} \
  -o BAM --ref-paths GRCh38.path_list.txt \
  -P -L 3000 \
  -f ${FASTQ_R1} \
  -f ${FASTQ_R2} \
  -Z ${GBZ} \
  --kff-name ${SAMPLE_NAME}.fq.kff \
  --haplotype-name ${HAPL} \
  -t ${THREADS} > ${SAMPLE_NAME}.unsorted.bam

echo "sorting bam ..."
BAM=${SAMPLE_NAME}.${OUTPUT_SUFFIX}.bam
time samtools view -h ${SAMPLE_NAME}.unsorted.bam | sed -e "s/GRCh38#0#//g" | samtools sort --threads 8 -m 2G -O BAM > ${BAM}

echo "removing unsorted bam..."
rm -rf ${SAMPLE_NAME}.unsorted.bam

echo "indexing bam ..."
# Index the BAM.
samtools index -@8 ${BAM}
