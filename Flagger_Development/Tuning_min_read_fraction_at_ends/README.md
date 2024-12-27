# Tuning min read fraction for 3 sequencing platforms; HiFi, R9 and R10 (--minReadFractionAtEnds, -f)

## Overview

HMM-Flagger-v1.2 introduces a new parameter designed to reduce false positive calls at contig ends. This parameter specifies the minimum fraction of a read that must be mapped to a contig end for the mapper to consider it a valid mapping. By adjusting this fraction, users can better manage coverage expectations at contig ends. Most of the false positive calls addressed by this parameter belong to the false duplication category, arising from naturally lower read coverage at contig ends rather than actual misassemblies.

For simplicity I restricted this tuning experiment to the reads mapped with minimap2 to the validation contigs of the HG002-T2T assembly falsified with a misassembly rate of 2%. We already have that from the output of the workflow that I used for tuning the alpha matrix.

## Steps
For tuning `--minReadFractionAtEnds` I wrote a bash script that does these steps and can be executed on SLURM cluster:
- Takes a json file that contains the necessary files and parameters for each platform 
- Creates a binary file from the given coverage file to accelerate running HMM-Flagger multiple times
- Uses fasta index to make a bed file covering at most 50kb of contig ends. It will cover the whole contig if the contig was shorter than 100kb.
- Runs HMM-Flagger for different values of `--minReadFractionAtEnds` ranging from 0 to 1.0 with a step of 0.05 (21 different values)
- Counts the number of bases discordant with the truth bed file (given as an input) only on the 50kb of contig ends. The total number of FP and FN bases for each category will be reported in a tsv file

## Code and input jsons
The bash script is available here.
https://github.com/mobinasri/research_notes/blob/main/Flagger_Development/Tuning_min_read_fraction_at_ends/run_hmm_flagger_tune_f.sh 

The json files:
- [HiFi josn](https://github.com/mobinasri/research_notes/blob/main/Flagger_Development/Tuning_min_read_fraction_at_ends/HiFi_DC_1.2/run_hmm_flagger_tune_f.inputs.json)
- [ONT-R10 json](https://github.com/mobinasri/research_notes/blob/main/Flagger_Development/Tuning_min_read_fraction_at_ends/ONT_R1041_Dorado/run_hmm_flagger_tune_f.inputs.json)
- [ONT-R9 json](https://github.com/mobinasri/research_notes/blob/main/Flagger_Development/Tuning_min_read_fraction_at_ends/ONT_R941_Guppy6.3.7/run_hmm_flagger_tune_f.inputs.json)


## The output tsv files
- [HiFi tsv](https://github.com/mobinasri/research_notes/blob/main/Flagger_Development/Tuning_min_read_fraction_at_ends/HiFi_DC_1.2/tune_read_frac_table_HiFi_DC_1.2.tsv)
- [ONT-R10 tsv](https://github.com/mobinasri/research_notes/blob/main/Flagger_Development/Tuning_min_read_fraction_at_ends/ONT_R1041_Dorado/tune_read_frac_table_ONT_R1041_Dorado.tsv)
- [ONT-R9 tsv](https://github.com/mobinasri/research_notes/blob/main/Flagger_Development/Tuning_min_read_fraction_at_ends/ONT_R941_Guppy6.3.7/tune_read_frac_table_ONT_R941_Guppy6.3.7.tsv)
