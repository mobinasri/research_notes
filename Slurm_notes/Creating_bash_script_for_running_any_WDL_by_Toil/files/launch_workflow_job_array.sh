#!/bin/bash

# Author      : Mobin Asri, masri@ucsc.edu
# Description : Launch toil job submission for any WDL-based workflow using Slurm arrays (Based on Julian's script)
# Usage       : sbatch launch_wdl_job_array.sh sample_file.csv
#               	sample_file.csv should have a header (otherwised first sample will be skipped)
#					and the sample names should be in the first column

#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mail-type=FAIL,END
#SBATCH --mem=16gb
#SBATCH --time=60:00

## Pull samples names from CSV passed to script
sample_file=$1

## Pull WDL file to be used for creating job array
wdl_path=$2
wdl_name=$(basename ${wdl_path%%.wdl})

## Pull extra arguments for Toil
## It can be something like "--time=3-0:00 --partition=high_priority"
export TOIL_SLURM_ARGS="$3"

## The directory that contains the input json files
input_json_dir=$4


# Read the CSV file and extract the sample ID for the current job array task
# Skip first row to avoid the header
sample_id=$(awk -F ',' -v task_id=${SLURM_ARRAY_TASK_ID} 'NR>1 && NR==task_id+1 {print $1}' "${sample_file}")

# Ensure a sample ID is obtained
if [ -z "${sample_id}" ]; then
    echo "Error: Failed to retrieve a valid sample ID. Exiting."
    exit 1
fi

echo "${sample_id}"

## Create then change into sample directory...
mkdir -p ${sample_id}
cd ${sample_id}


mkdir ${wdl_name}_logs 
mkdir analysis

export SINGULARITY_CACHEDIR=`pwd`/outputs/cache/.singularity/cache 
export MINIWDL__SINGULARITY__IMAGE_CACHE=`pwd`/outputs/cache/.cache/miniwdl 


toil-wdl-runner \
    --jobStore ./${wdl_name}_bigstore \
    --batchSystem slurm \
    --batchLogsDir ./${wdl_name}_logs \
    ${wdl_path} \
    ${input_json_dir}/${sample_id}_${wdl_name}.json \
    --outputDirectory analysis/${wdl_name} \
    --outputFile ${sample_id}_${wdl_name}_outputs.json \
    --runLocalJobsOnWorkers \
    --retryCount 0 \
    --disableProgress \
    2>&1 | tee log.txt

wait
echo "Done."
