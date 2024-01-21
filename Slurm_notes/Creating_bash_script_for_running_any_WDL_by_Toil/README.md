Comment 1: 01/21/2024

I created a new bash script for running any WDL by Toil when we have a data table and we want run the related workflow on each row of the table.
Julian has already written a wdl for hifiasm [here](https://github.com/human-pangenomics/hprc_intermediate_assembly/blob/main/assembly/batch1/launch_hifiasm_array.sh) 
I just modified that script to make it work with any WDL and take the path to WDL as an argument.
