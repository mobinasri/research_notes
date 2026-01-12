#!/bin/bash
#SBATCH --job-name=run_filter_hom_ref_masri
#SBATCH --partition=long
#SBATCH --mail-user=masri@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mail-type=FAIL,END
#SBATCH --mem=8gb
#SBATCH --cpus-per-task=16
#SBATCH --output=%x.%j.log
#SBATCH --time=48:00:00

set -euo pipefail

source /private/groups/patenlab/masri/haplotype_sampling/pangenome_aware_dv_paper/1KG_analysis/pangenome_aware_dv/vcf_env/bin/activate

N_PARALLEL=4
THREADS=4
WORKING_DIR="/private/groups/patenlab/masri/haplotype_sampling/pangenome_aware_dv_paper/1KG_analysis/linear_ref_based_dv/5_HPRC_samples"

cd "${WORKING_DIR}"

SAMPLES="HG01255,HG02280,HG02984,HG03831,HG04184"
GS_PREFIX="gs://brain-genomics-public/research/cohort/1KGP/vg/dv_grch38"

create_vcf_for_chrom() {
  local chrom="$1"

  local in="chr${chrom}.vcf.gz"
  local out="chr${chrom}.1kg_cohort.vg.linear_ref_based_dv.grch38.5_HPRC_samples.with_non_ref_allele.vcf.gz"

  echo "[$(date)] START chr${chrom}"

  gsutil -m cp "${GS_PREFIX}/${in}" ${WORKING_DIR}/
  gsutil -m cp "${GS_PREFIX}/${in}.tbi" ${WORKING_DIR}/

  python3 ${WORKING_DIR}/filter_hom_ref.py \
    --input "${WORKING_DIR}/${in}" \
    --output "${WORKING_DIR}/${out}" \
    --samples "${SAMPLES}" \
    --threads "${THREADS}"

  rm -f ${WORKING_DIR}/${in}
  rm -f ${WORKING_DIR}/${in}.tbi

  echo "[$(date)] DONE chr${chrom}"
}

export WORKING_DIR
export -f create_vcf_for_chrom
export THREADS SAMPLES GS_PREFIX

# Run up to N_PARALLEL chromosomes concurrently
printf "%s\n" $(seq 1 22) X Y \
  | xargs -n 1 -P "${N_PARALLEL}" -I {} bash -lc 'create_vcf_for_chrom "$@"' _ {}

echo "All chromosomes completed."

