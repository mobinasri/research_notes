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
WORKING_DIR="/private/groups/patenlab/masri/haplotype_sampling/pangenome_aware_dv_paper/1KG_Trio_analysis/pangenome_aware_dv"

cd "${WORKING_DIR}"

SAMPLES_CHR_12_22_XY=$(cat ../selected_trio_samples_for_pangenome_aware_dv.txt)
SAMPLES_CHR_1_11=$(cat ../selected_trio_samples.txt)

GS_PREFIX="gs://brain-genomics-public/research/cohort/1KGP/vg/pangenome_aware_to_grch38/glnexus_output"

create_vcf_for_chrom() {
  local chrom="$1"

  # Select sample list based on chromosome
  if [[ "$chrom" =~ ^([1-9]|1[01])$ ]]; then
    local SAMPLES="${SAMPLES_CHR_1_11}"
  else
    local SAMPLES="${SAMPLES_CHR_12_22_XY}"
  fi

  local in="chr${chrom}.1kg_cohort.vg.deepvariant.grch38.vcf.gz"
  local out="chr${chrom}.1kg_cohort.vg.pangenome_aware_dv.grch38.20_trios.with_non_ref_allele.vcf.gz"

  echo "[$(date)] START chr${chrom}"

  gsutil -m cp "${GS_PREFIX}/${in}" ${WORKING_DIR}/
  gsutil -m cp "${GS_PREFIX}/${in}.tbi" ${WORKING_DIR}/

  python3 ${WORKING_DIR}/../filter_hom_ref.py \
    --input "${WORKING_DIR}/${in}" \
    --output "${WORKING_DIR}/${out}" \
    --samples "${SAMPLES}" \
    --threads "${THREADS}"

  rm -f ${WORKING_DIR}/${in}
  rm -f ${WORKING_DIR}/${in}.tbi

  tabix ${WORKING_DIR}/${out}

  echo "[$(date)] DONE chr${chrom}"
}

export WORKING_DIR
export -f create_vcf_for_chrom
export THREADS SAMPLES_CHR_1_11 SAMPLES_CHR_12_22_XY GS_PREFIX

# Run up to N_PARALLEL chromosomes concurrently
printf "%s\n" $(seq 1 22) X Y \
  | xargs -n 1 -P "${N_PARALLEL}" -I {} bash -lc 'create_vcf_for_chrom "$@"' _ {}

echo "All chromosomes completed."

