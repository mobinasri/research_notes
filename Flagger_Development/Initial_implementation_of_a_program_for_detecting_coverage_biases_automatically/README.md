## Initial implementation of a program for detecting coverage biases automatically

### Comment 1: 01/31/2023

#### Motivation
In evaluating CHM13-v1.0, HG002-T2T-v1.0.1 and HPRC-Year1 assemblies it has been observed clearly that a number of Human Satellites have biases in read depth of coverage.
Moreover these biases are varying across different sequencing platforms. For example in the recent investigation of the read alignments to HG002-T2T-v1.0.1 I noticed that 
for R10-UL data basecalled with Dorado the read coverage across HSat1B and HSat1A arrays are significantly more than its average over other parts of the genome (Look at
[this slide deck](https://docs.google.com/presentation/d/1vZEPkC3NyOxyIhlaQX26oltio8rALxeRvwmbvj5-9MY/edit#slide=id.g2b2d09284bb_0_185)).

#### Idea
Since each platform may have its specific biases in different HSats and it's time-consuming to identify them manually I want to write a program that does this automatically. This program should take a BAM file and also a list of BED files each of which contains the regions with a specific annotation. For example one BED file can point to all the regions annotated as HSat1A. It will basically calculate the read coverages along each annotation and find the coverage value with the highest frequency. This program also needs a BED file of all the remaining regions in the genome to use as a baseline. It will find the most frequent coverage value in the "baseline" regions. By comparing the most frequent coverage values in "baseline" vs other annotations we can find any significant deviation.

One key assumption here is that the most frequent coverage value is an estimation of the expected coverage in each annotation. It is a valid assumption as long as more than half of each annotation is assembled correctly. Therefore for assemblies that satellite arrays are two fragmented, erronoues or they are targets for evaluation there is not benefit in using this program.
