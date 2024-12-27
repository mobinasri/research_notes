# Tuning min read fraction for 3 sequencing platforms; HiFi, R9 and R10 (--minReadFractionAtEnds, -f)

## Overview

HMM-Flagger-v1.2 will have a new parameter for fixing false positive calls at contig ends. This parameter determines the minimum fraction of a read that 
should be mapped to a contig end so that the mapper will report it as a valid mapping. This fraction can be used for adjusting coverage expectations at contig ends. 
Most of the false positive calls that are expected to be removed with this parameter are ‘Dup’ blocks called only because of having lower coverage at contig ends 
not because of the existence of a real misassembly.

For simplicity I restrict this tuning analysis to the reads mapped with minimap2 to the validation contigs of the falsified assembly with a misassembly rate of 2%.
We already have that from the output of the workflow for tuning the alpha matrix.


```

```
