## Initial implementation of a program for detecting coverage biases automatically

### Comment 1: 01/31/2023

#### Motivation
In evaluating CHM13-v1.0, HG002-T2T-v1.0.1 and HPRC-Year1 assemblies it has been observed clearly that a number of Human Satellites have biases in read depth of coverage.
Moreover these biases are varying across different sequencing platforms. For example in the recent investigation of the read alignments to HG002-T2T-v1.0.1 I noticed that 
for R10-UL data basecalled with Dorado the read coverage across HSat1B and HSat1A arrays are significantly more than its average over other parts of the genome (Look at
[this slide deck](https://docs.google.com/presentation/d/1vZEPkC3NyOxyIhlaQX26oltio8rALxeRvwmbvj5-9MY/edit#slide=id.g2b2d09284bb_0_185)).

#### Idea
Since each platform may have its specific biases in different HSats and it's time-consuming to identify them manually I want to write a program that does this automatically. This program should take a BAM file and also a list of BED files each of which contains the regions with a specific annotation. For example one BED file can point to all the regions annotated as HSat1A. It will basically calculate the read coverages along each annotation and find the coverage value with the highest frequency. This program also needs a BED file of all the remaining regions in the genome to use as a baseline. It will find the most frequent coverage value in the "baseline" regions. By comparing the most frequent coverage values in "baseline" vs other annotations we can find any significant deviation.

One key assumption here is that the most frequent coverage value is an estimation of the expected coverage in each annotation. It is a valid assumption as long as more than half of each annotation is assembled correctly. Therefore for assemblies in which the satellite arrays are highly fragmented and erroneous this program might not be useful.

#### The initial implementation

Here is the link:
https://github.com/mobinasri/flagger/blob/dev-0.3.0/programs/src/bias_detector.c

I added `bias_detector` to the Flagger docker image `mobinasri/flagger:dev-v0.3.3`

#### Making small bam file for testing bias_detector
Using samtools I took a subset of the R10 UL data basecalled with Dorado and mapped to HG002-T2T-v1.0.1. I selected the p-arm of chr15_paternal since I know that this data has bias in HSAT1B. To make use of 
the multi-threading feature of `bias-detector` I had to update the header of the subset BAM file to only
contain the chr15_PATERNAL contig. I did this manually. 

(older versions of samtools could not pull data from s3 bucket so I had to update it to v0.19.2)

```
# pull data from s3 bucket and subset to chr15:1-20mb
./samtools view -hb https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/polishing/HG002/v1.0/mapping/ont_r10_ul_dorado/hg002v1.0_ont_r10_ul_dorado.pri.bam chr15_PATERNAL:1-20000000 > /private/groups/patenlab/masri/t2t/HG002_v1.0.1/hg002v1.0_ont_r10_ul_dorado.pri.chr15_pat_1_20mb.bam

# take the header and remove any contig other than chr15_PATERNAL
samtools view -H /private/groups/patenlab/masri/t2t/HG002_v1.0.1/hg002v1.0_ont_r10_ul_dorado.pri.chr15_pat_1_20mb.bam > /private/groups/patenlab/masri/t2t/HG002_v1.0.1/hg002v1.0_ont_r10_ul_dorado.pri.chr15_pat_1_20mb.header.sam

# open hg002v1.0_ont_r10_ul_dorado.pri.chr15_pat_1_20mb.header.sam and modify the header
vim /private/groups/patenlab/masri/t2t/HG002_v1.0.1/hg002v1.0_ont_r10_ul_dorado.pri.chr15_pat_1_20mb.header.sam

# take all records with no header
samtools view /private/groups/patenlab/masri/t2t/HG002_v1.0.1/hg002v1.0_ont_r10_ul_dorado.pri.chr15_pat_1_20mb.bam > /private/groups/patenlab/masri/t2t/HG002_v1.0.1/hg002v1.0_ont_r10_ul_dorado.pri.chr15_pat_1_20mb.no_header.sam


# make a new bam file with modified header
cat /private/groups/patenlab/masri/t2t/HG002_v1.0.1/hg002v1.0_ont_r10_ul_dorado.pri.chr15_pat_1_20mb.header.sam /private/groups/patenlab/masri/t2t/HG002_v1.0.1/hg002v1.0_ont_r10_ul_dorado.pri.chr15_pat_1_20mb.no_header.sam | samtools view -hb > /private/groups/patenlab/masri/t2t/HG002_v1.0.1/hg002v1.0_ont_r10_ul_dorado.pri.chr15_pat_1_20mb.bam

# index the new bam file
samtools index /private/groups/patenlab/masri/t2t/HG002_v1.0.1/hg002v1.0_ont_r10_ul_dorado.pri.chr15_pat_1_20mb.bam

```

So `/private/groups/patenlab/masri/t2t/HG002_v1.0.1/hg002v1.0_ont_r10_ul_dorado.pri.chr15_pat_1_20mb.bam` is the BAM file I want to test my program with.

Next I have to make annotation and baseline BED files to pass to `bias_detector`.

