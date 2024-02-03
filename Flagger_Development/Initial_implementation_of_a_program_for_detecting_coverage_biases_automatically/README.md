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

#### Making a small bam file for testing bias_detector
Using samtools I took a subset of the R10 UL data (basecalled with Dorado) that was mapped to HG002-T2T-v1.0.1. I selected the p-arm of chr15_paternal since I observed that this data has bias in HSAT1A/B. To make use of 
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

I made these 4 small BED files:
```
cat hg002v1.0.1.cenSatv1.0.chr15_pat_hsat1A.bed
chr15_PATERNAL	1298768	1850810
chr15_PATERNAL	7423750	7428399

cat hg002v1.0.1.cenSatv1.0.chr15_pat_hsat1B.bed
chr15_PATERNAL	6862343	7423750

cat hg002v1.0.1.cenSatv1.0.chr15_pat_asat.bed
chr15_PATERNAL	11654483	11660121	hor(S4/6C13/14/21/22H8)	100	.	11654483	11660121	255,146,0
chr15_PATERNAL	11675562	11676762	hor(S6C22H2-A)	100	.	11675562	11676762	255,146,0
chr15_PATERNAL	11825134	11826331	hor(S6C22H2-A)	100	.	11825134	11826331	255,146,0
chr15_PATERNAL	12654013	12657095	hor(S6C22H2-A)	100	.	12654013	12657095	255,146,0
chr15_PATERNAL	12950739	12967647	hor(S4C15H3-B)	100	.	12950739	12967647	255,146,0
chr15_PATERNAL	12967816	13021541	hor(S4C15H3-A)	100	.	12967816	13021541	255,146,0
chr15_PATERNAL	13021541	13067989	hor(S4C15H3-B)	100	.	13021541	13067989	255,146,0
chr15_PATERNAL	13068162	13135094	hor(S4C15H3-A)	100	.	13068162	13135094	255,146,0
chr15_PATERNAL	13201474	13202668	hor(S6C22H2-A)	100	.	13201474	13202668	255,146,0
chr15_PATERNAL	13324020	13325220	hor(S6C22H2-A)	100	.	13324020	13325220	255,146,0
chr15_PATERNAL	13438651	14073759	hor(S4C15H2)	100	.	13438651	14073759	255,146,0
chr15_PATERNAL	14182824	14186846	dhor(S02CMH2d)	100	.	14182824	14186846	153,0,0
chr15_PATERNAL	14193748	14196988	dhor(S2CMH4d)	100	.	14193748	14196988	153,0,0
chr15_PATERNAL	14196988	15332202	active_hor(S2C15H1L)	100	.	14196988	15332202	250,0,0


cat baseline.no_rDNA.no_hsat1.no_asat.chr15_pat_1_20mb.bed
chr15_PATERNAL	0	1298768
chr15_PATERNAL	1850810	2920268
chr15_PATERNAL	5889932	6862343
chr15_PATERNAL	7428399	11654483
chr15_PATERNAL	11660121	11675562
chr15_PATERNAL	11676762	11825134
chr15_PATERNAL	11826331	12654013
chr15_PATERNAL	12657095	12950739
chr15_PATERNAL	12967647	12967816
chr15_PATERNAL	13067989	13068162
chr15_PATERNAL	13135094	13201474
chr15_PATERNAL	13202668	13324020
chr15_PATERNAL	13325220	13438651
chr15_PATERNAL	14073759	14182824
chr15_PATERNAL	14186846	14193748
chr15_PATERNAL	15332202	20000000
```

Then I made a json file with annotation names and paths:
```
cat bed_files.json 
{
	"baseline":"/private/groups/patenlab/masri/t2t/HG002_v1.0.1/test_bias_detector/baseline.no_rDNA.no_hsat1.no_asat.chr15_pat_1_20mb.bed" ,
	"asat" : "/private/groups/patenlab/masri/t2t/HG002_v1.0.1/test_bias_detector/hg002v1.0.1.cenSatv1.0.chr15_pat_asat.bed",
	"hsat1A":"/private/groups/patenlab/masri/t2t/HG002_v1.0.1/test_bias_detector/hg002v1.0.1.cenSatv1.0.chr15_pat_hsat1A.bed" ,
	"hsat1B":"/private/groups/patenlab/masri/t2t/HG002_v1.0.1/test_bias_detector/hg002v1.0.1.cenSatv1.0.chr15_pat_hsat1B.bed"

}

```

Now I can run `bias_detector` with this command:
```
cd /private/groups/patenlab/masri/t2t/HG002_v1.0.1/test_bias_detector

bias_detector -i hg002v1.0_ont_r10_ul_dorado.pri.chr15_pat_1_20mb.bam -j bed_files.json -t4 -b "baseline" -d 0.2 > bias_table.tsv
```
Here is the output log
```
{
	"baseline":"/private/groups/patenlab/masri/t2t/HG002_v1.0.1/test_bias_detector/baseline.no_rDNA.no_hsat1.no_asat.chr15_pat_1_20mb.bed" ,
	"asat" : "/private/groups/patenlab/masri/t2t/HG002_v1.0.1/test_bias_detector/hg002v1.0.1.cenSatv1.0.chr15_pat_asat.bed",
	"hsat1A":"/private/groups/patenlab/masri/t2t/HG002_v1.0.1/test_bias_detector/hg002v1.0.1.cenSatv1.0.chr15_pat_hsat1A.bed" ,
	"hsat1B":"/private/groups/patenlab/masri/t2t/HG002_v1.0.1/test_bias_detector/hg002v1.0.1.cenSatv1.0.chr15_pat_hsat1B.bed"

}
[2024-02-03 07:07:58] Parsed  annotation baseline:/private/groups/patenlab/masri/t2t/HG002_v1.0.1/test_bias_detector/baseline.no_rDNA.no_hsat1.no_asat.chr15_pat_1_20mb.bed
[2024-02-03 07:07:58] Parsed  annotation asat:/private/groups/patenlab/masri/t2t/HG002_v1.0.1/test_bias_detector/hg002v1.0.1.cenSatv1.0.chr15_pat_asat.bed
[2024-02-03 07:07:58] Parsed  annotation hsat1A:/private/groups/patenlab/masri/t2t/HG002_v1.0.1/test_bias_detector/hg002v1.0.1.cenSatv1.0.chr15_pat_hsat1A.bed
[2024-02-03 07:07:58] Parsed  annotation hsat1B:/private/groups/patenlab/masri/t2t/HG002_v1.0.1/test_bias_detector/hg002v1.0.1.cenSatv1.0.chr15_pat_hsat1B.bed
[2024-02-03 07:07:58] Number of parsed annotations = 4
[2024-02-03 07:07:58] Created thread for parsing batch 0 (length=24064255)
[2024-02-03 07:07:58] Created thread for parsing batch 1 (length=24064255)
[2024-02-03 07:07:58] Created thread for parsing batch 2 (length=24064255)
[2024-02-03 07:07:58] Created thread for parsing batch 3 (length=24064252)
[2024-02-03 07:07:58] Started parsing reads in the block: chr15_PATERNAL	48128510	72192765
[2024-02-03 07:07:58] Started parsing reads in the block: chr15_PATERNAL	72192765	96257017
[2024-02-03 07:07:58] Started parsing reads in the block: chr15_PATERNAL	0	24064255
[2024-02-03 07:07:58] Started parsing reads in the block: chr15_PATERNAL	24064255	48128510
[2024-02-03 07:08:22] All batches are parsed.
[2024-02-03 07:08:22] Started sorting and merging blocks.
[2024-02-03 07:08:22] Merging coverage blocks is done.
[2024-02-03 07:08:22] Created block table with coverage data : tot_len=17589279, number=52421
[2024-02-03 07:08:22] Created block table for whole genome  : tot_len=96257017, number=1
[2024-02-03 07:08:22] Created block table for whole genome  : tot_len=96257017, number=1
[2024-02-03 07:08:22] Added 0-coverage whole genome blocks to coverage block tables : tot_len=113846296, number=52422
[2024-02-03 07:08:22] Added annotation blocks to coverage block tables: tot_len=130876632, number=52455
[2024-02-03 07:08:22] Started sorting and merging blocks
[2024-02-03 07:08:22] Merged blocks : tot_len=130876632, number=52455
[2024-02-03 07:08:22] Created final block table : tot_len=96257017, number=52457
{
	"baseline":"/private/groups/patenlab/masri/t2t/HG002_v1.0.1/test_bias_detector/baseline.no_rDNA.no_hsat1.no_asat.chr15_pat_1_20mb.bed" ,
	"asat" : "/private/groups/patenlab/masri/t2t/HG002_v1.0.1/test_bias_detector/hg002v1.0.1.cenSatv1.0.chr15_pat_asat.bed",
	"hsat1A":"/private/groups/patenlab/masri/t2t/HG002_v1.0.1/test_bias_detector/hg002v1.0.1.cenSatv1.0.chr15_pat_hsat1A.bed" ,
	"hsat1B":"/private/groups/patenlab/masri/t2t/HG002_v1.0.1/test_bias_detector/hg002v1.0.1.cenSatv1.0.chr15_pat_hsat1B.bed"

}
{
	"baseline":"/private/groups/patenlab/masri/t2t/HG002_v1.0.1/test_bias_detector/baseline.no_rDNA.no_hsat1.no_asat.chr15_pat_1_20mb.bed" ,
	"asat" : "/private/groups/patenlab/masri/t2t/HG002_v1.0.1/test_bias_detector/hg002v1.0.1.cenSatv1.0.chr15_pat_asat.bed",
	"hsat1A":"/private/groups/patenlab/masri/t2t/HG002_v1.0.1/test_bias_detector/hg002v1.0.1.cenSatv1.0.chr15_pat_hsat1A.bed" ,
	"hsat1B":"/private/groups/patenlab/masri/t2t/HG002_v1.0.1/test_bias_detector/hg002v1.0.1.cenSatv1.0.chr15_pat_hsat1B.bed"

}
Printing coverage table in stdout ...
[2024-02-03 07:08:22] Done.
```

#### Output table

```
cat bias_table.tsv

annotation	status	most_freq_cov	cov_diff_normalized	path
baseline	not_biased	55	+0.000	/private/groups/patenlab/masri/t2t/HG002_v1.0.1/test_bias_detector/baseline.no_rDNA.no_hsat1.no_asat.chr15_pat_1_20mb.bed
asat	not_biased	54	-0.018	/private/groups/patenlab/masri/t2t/HG002_v1.0.1/test_bias_detector/hg002v1.0.1.cenSatv1.0.chr15_pat_asat.bed
hsat1A	biased	98	+0.782	/private/groups/patenlab/masri/t2t/HG002_v1.0.1/test_bias_detector/hg002v1.0.1.cenSatv1.0.chr15_pat_hsat1A.bed
hsat1B	biased	135	+1.455	/private/groups/patenlab/masri/t2t/HG002_v1.0.1/test_bias_detector/hg002v1.0.1.cenSatv1.0.chr15_pat_hsat1B.bed
```

