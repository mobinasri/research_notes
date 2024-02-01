## Initial implementation of a program for detecting coverage biases automatically

### Comment 1: 01/31/2023

In evaluating CHM13-v1.0, HG002-T2T-v1.0.1 and HPRC-Year1 assemblies it has been observed clearly that a number of Human Satellites have different types of biases in read depth of coverage.
Moreover these biases are varying across different sequencing platforms. For example in the recent investigation of the read alignments to HG002-T2T-v1.0.1 I noticed that 
for R10-UL data basecalled with Dorado the read coverage across HSat1B and HSat1A arrays are significantly more than its average over other parts of the genome (Look at
[this slide deck](https://docs.google.com/presentation/d/1vZEPkC3NyOxyIhlaQX26oltio8rALxeRvwmbvj5-9MY/edit#slide=id.g2b2d09284bb_0_185)).
