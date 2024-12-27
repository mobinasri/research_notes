#!/bin/bash
#SBATCH --job-name=tune_read_frac_hmm_flagger
#SBATCH --partition=long
#SBATCH --mail-user=masri@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mail-type=FAIL,END
#SBATCH --mem=16gb
#SBATCH --cpus-per-task=16
#SBATCH --output=%x.%j.log
#SBATCH --time=10:00:00

WORKING_DIR=$(cat $1 | jq -r ".working_dir")
PLATFORM=$(cat $1 | jq -r ".platform")
FAI=$(cat $1 | jq -r ".fai")
COV_FILE_PATH=$(cat $1 | jq -r ".cov_file_path")
WINDOW_LEN=$(cat $1 | jq -r ".window_len")
MISASSEMBLY_TRUTH_BED=$(cat $1 | jq -r ".misassembly_truth_bed")
MODEL_TYPE=$(cat $1 | jq -r ".model_type")
FLAGGER_DOCKER=$(cat $1 | jq -r ".flagger_docker")

CHUNK_LEN=20000000
THREAD_COUNT=16
PREFIX=$(basename ${COV_FILE_PATH%%.cov.gz})

# go to working directory
cd ${WORKING_DIR}

# Extract collapsed blocks
cat ${MISASSEMBLY_TRUTH_BED} | grep "Col" > Col_TRUTH.bed

# Extract duplicated blocks
cat ${MISASSEMBLY_TRUTH_BED} | grep "Dup" > Dup_TRUTH.bed

# Extract erroneous blocks
cat ${MISASSEMBLY_TRUTH_BED} | grep "Err" | grep -v "Col" > Err_TRUTH.bed

# Make a bed file that covers only the first/last 50kb of the validation contigs
cat ${FAI} | \
    grep -e chr17 -e chr18 -e chr19 -e chr20 -e chr21 -e chr22 | \
    awk '$2<=50000{print $1"\t0\t"$2;next} {print $1"\t0\t"50e3; print $1"\t"$2-50e3"\t"$2}' > validation.ends_50kb.bed

docker run -u $(id -u):$(id -g) --rm -v/private/groups/patenlab/masri:/private/groups/patenlab/masri ${FLAGGER_DOCKER} \
	coverage_format_converter \
	-c ${CHUNK_LEN} \
	-w ${WINDOW_LEN} \
	-i ${COV_FILE_PATH} \
	-f ${FAI} \
	-t ${THREAD_COUNT} \
	-o ${WORKING_DIR}/${PREFIX}.bin

# print the column names of the final table
printf "#PLATFORM\tREAD_FRAC\tErr_TP\tErr_FP\tErr_FN\tDup_TP\tDup_FP\tDup_FN\tCol_TP\tCol_FP\tCol_FN\tTotal_FN_FP\n" > tune_read_frac_table_${PLATFORM}.tsv

for READ_FRAC in $(seq 0.0 0.05 1.0); do
	mkdir -p min_read_frac_${READ_FRAC}
	docker run -u $(id -u):$(id -g) --rm -v/private/groups/patenlab/masri:/private/groups/patenlab/masri ${FLAGGER_DOCKER} \
	hmm_flagger \
            --input ${WORKING_DIR}/${PREFIX}.bin \
            --outputDir ${WORKING_DIR}/min_read_frac_${READ_FRAC}  \
            --chunkLen ${CHUNK_LEN} \
            --windowLen ${WINDOW_LEN} \
            --modelType ${MODEL_TYPE} \
            --trackName  min_read_frac_${READ_FRAC} \
            --threads ${THREAD_COUNT} \
	    --adjustContigEnds \
	    --minReadFractionAtEnds ${READ_FRAC} 2>&1 | tee ${WORKING_DIR}/min_read_frac_${READ_FRAC}/hmm_flagger.log

	printf "${PLATFORM}\t${READ_FRAC}" >> tune_read_frac_table_${PLATFORM}.tsv 
	tot=0
	for LABEL in Err Dup Col
	do
		cat min_read_frac_${READ_FRAC}/final_flagger_prediction.bed | \
			grep ${LABEL} | \
			bedtools sort -i - | \
			bedtools intersect -a - -b ${LABEL}_TRUTH.bed | \
			bedtools intersect -a - -b validation.ends_50kb.bed | \
			awk 'BEGIN{s=0}{s+=$3-$2}END{print s}' > min_read_frac_${READ_FRAC}/${LABEL}_TP.txt
		
		cat min_read_frac_${READ_FRAC}/final_flagger_prediction.bed | \
			grep ${LABEL} | \
			bedtools sort -i - | \
			bedtools subtract -a - -b ${LABEL}_TRUTH.bed | \
			bedtools intersect -a - -b validation.ends_50kb.bed | \
			awk 'BEGIN{s=0}{s+=$3-$2}END{print s}' > min_read_frac_${READ_FRAC}/${LABEL}_FP.txt

		cat min_read_frac_${READ_FRAC}/final_flagger_prediction.bed | \
			grep ${LABEL} | \
			bedtools sort -i - | \
			bedtools subtract -a ${LABEL}_TRUTH.bed -b - | \
			bedtools intersect -a - -b validation.ends_50kb.bed | \
			awk 'BEGIN{s=0}{s+=$3-$2}END{print s}' > min_read_frac_${READ_FRAC}/${LABEL}_FN.txt

		tp=$(cat min_read_frac_${READ_FRAC}/${LABEL}_TP.txt)
		fp=$(cat min_read_frac_${READ_FRAC}/${LABEL}_FP.txt)
		fn=$(cat min_read_frac_${READ_FRAC}/${LABEL}_FN.txt)
		tot=$((tot + fn + fp))

		printf "\t"${tp}"\t"${fp}"\t"${fn} >> tune_read_frac_table_${PLATFORM}.tsv
	done
	printf "\t${tot}\n" >> tune_read_frac_table_${PLATFORM}.tsv

done
