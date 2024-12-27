# Tuning min read fraction for 3 sequencing platforms; HiFi, R9 and R10 (--minReadFractionAtEnds, -f)

## Overview

HMM-Flagger-v1.2 will have a new parameter for fixing false positive calls at contig ends. This parameter determines the minimum fraction of a read that 
should be mapped to a contig end so that the mapper will report it as a valid mapping. This fraction can be used for adjusting coverage expectations at contig ends. 
Most of the false positive calls that are expected to be removed with this parameter are from the false duplication category. They are called only because of the natural lower read coverage at contig ends but not because of the existence of a real misassembly.

For simplicity I restricted this tuning experiment to the reads mapped with minimap2 to the validation contigs of the HG002-T2T assembly falsified with a misassembly rate of 2%. We already have that from the output of the workflow that I used for tuning the alpha matrix.

For tuning `--minReadFractionAtEnds` I wrote a bash script that does these steps and can be executed on SLURM cluster:
- Takes a json file that contains the necessary files and parameters for each platform 
- Creates a binary file from the given coverage file to accelerate running HMM-Flagger multiple times
- Runs HMM-Flagger for different values of `--minReadFractionAtEnds` ranging from 0 to 1.0 with a step of 0.05 (21 different values)
- Counts the number of bases discordant with the truth bed file (given as an input). The number of FP and FN bases for each category will be reported in a tsv file

The bash script is available here.
