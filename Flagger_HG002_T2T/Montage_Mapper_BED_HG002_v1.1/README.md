## Creating mapping BED file for montage_mapper.py

Here we want to create a 7-column BED file that specifies the mapper for different pairs of blocks in the HG002-T2T-v1.1 assembly. The mapper can be either minimap2 or centrolign.


Assembly v1.1 is taken from https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/HG002/assemblies/ 


### Split chromosomes into p-arm, peri/centromeric and q-arm to align them with different mappers.

- p-arm and q-arm of non-acro chromosomes with minimap2-2.26
- peri/centromeric (including HSat 2/3 and ASat) with centrolign
- P-arms of acrocentric chromosomes with centrolign

Censat annotation, hg002v1.1_v2.0/hg002v1.1.cenSatv2.0.bed, is taken from this area
https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/HG002/assemblies/annotation/centromere/

### Make whole genome bed

```
cd /private/groups/patenlab/masri/t2t/HG002_v1.1/assembly
cat hg002v1.1.fasta.fai | \
    awk '{print $1"\t0\t"$2}' | \
    bedtools sort -i - > hg002v1.1.bed
```


### Download annotation
```
cd /private/groups/patenlab/masri/t2t/HG002_v1.1/annotations/censat

wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/annotation/centromere/hg002v1.1_v2.0/hg002v1.1.cenSatv2.0.bed
```





### Get rDNA coordinates and merge them with “N” GAP coordinates
```
cat hg002v1.1.cenSatv2.0.bed | \
    grep -e rDNA -e GAP | \
    cut -f1-3  | \
    awk '$3-$2 > 10e3' | \
    bedtools merge -i - > hg002v1.1.cenSatv2.0.rDNA_plus_GAP.bed
```

### Get p-arm coordinates of the acrocentric chromosomes
```
cat hg002v1.1.cenSatv2.0.rDNA_plus_GAP.bed | \
    awk '{print $1"\t"0"\t"$2}' | \
    bedtools merge -i - | \
    grep -v -e chrY_ -e chr2_ > acro_parm.bed
```


### Keep stretches of HOR or HSat longer than 10kb, merge blocks closer than 1Mb and keep final merged blocks longer than 1Mb
```
cat hg002v1.1.cenSatv2.0.bed | \
    grep -e "hor" -e "hsat" -e "HSat" | \
    awk '10e3 < ($3-$2)' | \
    bedtools merge -i - -d 1e6 | \
    awk '1e6 < ($3-$2)' > hg002v1.1.cenSatv2.0.only_hor_hsat.gt_10k.merged_1mb.gt_1mb.bed
```






### Merge blocks with rDNA and surrounding N blocks and then subtract merged blocks from them. This is for keeping the bases upstream and downstream of rDNA but not rDNA themselves

```
cat hg002v1.1.cenSatv2.0.rDNA_plus_GAP.bed hg002v1.1.cenSatv2.0.only_hor_hsat.gt_10k.merged_1mb.gt_1mb.bed | \
    bedtools sort -i - | \
    bedtools merge -i - -d 1e6 | \
    bedtools subtract -a - -b hg002v1.1.cenSatv2.0.rDNA_plus_GAP.bed > hg002v1.1.cenSatv2.0.only_hor_hsat.gt_10k.merged_1mb.gt_1mb.adjusted_by_rdna.bed
```

### Added p-arms of acros

```
cat hg002v1.1.cenSatv2.0.only_hor_hsat.gt_10k.merged_1mb.gt_1mb.adjusted_by_rdna.bed  acro_parm.bed | \
    bedtools sort -i -  | \
    bedtools merge -i -  > hg002v1.1.cenSatv2.0.only_hor_hsat.gt_10k.merged_1mb.gt_1mb.adjusted_by_rdna.added_acro_parm.bed
```

### Remove sex chromosomes 

```
cat hg002v1.1.cenSatv2.0.only_hor_hsat.gt_10k.merged_1mb.gt_1mb.adjusted_by_rdna.added_acro_parm.bed | \
    grep -v "chrX" | \
    grep -v "chrY" > hg002v1.1.cenSatv2.0.only_hor_hsat.gt_10k.merged_1mb.gt_1mb.adjusted_by_rdna.added_acro_parm.no_sex.bed
```

### Separate pat and mat blocks

```
cat hg002v1.1.cenSatv2.0.only_hor_hsat.gt_10k.merged_1mb.gt_1mb.adjusted_by_rdna.added_acro_parm.no_sex.bed | grep "PAT" > centrolign_blocks.pat.bed

cat hg002v1.1.cenSatv2.0.only_hor_hsat.gt_10k.merged_1mb.gt_1mb.adjusted_by_rdna.added_acro_parm.no_sex.bed | grep "MAT" > centrolign_blocks.mat.bed
```



### Pair blocks that have to mapped with centrolign

```
paste centrolign_blocks.pat.bed centrolign_blocks.mat.bed | awk '{print $0"\tcentrolign"}' > centrolign_blocks.pairs.bed
```

```
cat centrolign_blocks.pairs.bed

chr10_PATERNAL	38578643	44232421	chr10_MATERNAL	38577130	43635884	centrolign
chr11_PATERNAL	50964481	53508439	chr11_MATERNAL	50938529	54665772	centrolign
chr12_PATERNAL	34580373	37782453	chr12_MATERNAL	34545968	37778542	centrolign
chr13_PATERNAL	0	2503151	chr13_MATERNAL	0	3646400	centrolign
chr13_PATERNAL	2779094	8564979	chr13_MATERNAL	9779362	12601132	centrolign
chr13_PATERNAL	9679403	13122311	chr13_MATERNAL	13765549	17236289	centrolign
chr14_PATERNAL	0	3424033	chr14_MATERNAL	0	3467031	centrolign
chr14_PATERNAL	5213545	10495661	chr14_MATERNAL	8664438	12419632	centrolign
chr14_PATERNAL	11578236	16762214	chr14_MATERNAL	13442573	18242469	centrolign
chr15_PATERNAL	0	2920268	chr15_MATERNAL	0	2240520	centrolign
chr15_PATERNAL	5895530	15332197	chr15_MATERNAL	5825891	18803957	centrolign
chr16_PATERNAL	34686993	37402212	chr16_MATERNAL	35888014	38204872	centrolign
chr16_PATERNAL	38566091	49744461	chr16_MATERNAL	39373481	46509338	centrolign
chr17_PATERNAL	21837775	28088572	chr17_MATERNAL	21673042	27285530	centrolign
chr18_PATERNAL	15476428	21348472	chr18_MATERNAL	15477187	19676097	centrolign
chr19_PATERNAL	24653975	29564703	chr19_MATERNAL	24693385	29421838	centrolign
chr1_PATERNAL	121897879	146932593	chr1_MATERNAL	121793663	137240935	centrolign
chr20_PATERNAL	26365907	33807407	chr20_MATERNAL	26299400	32960966	centrolign
chr21_PATERNAL	0	2558579	chr21_MATERNAL	0	4193474	centrolign
chr21_PATERNAL	5130853	7284755	chr21_MATERNAL	8143147	10929182	centrolign
chr21_PATERNAL	8402826	10490173	chr21_MATERNAL	12106868	13506641	centrolign
chr22_PATERNAL	0	2495861	chr22_MATERNAL	0	2681747	centrolign
chr22_PATERNAL	3848100	15001301	chr22_MATERNAL	8271869	18892213	centrolign
chr2_PATERNAL	88976293	94232987	chr2_MATERNAL	89056701	94465377	centrolign
chr3_PATERNAL	90555532	96708207	chr3_MATERNAL	90461937	96202232	centrolign
chr4_PATERNAL	49182798	54080686	chr4_MATERNAL	49229606	53502364	centrolign
chr5_PATERNAL	46686377	57738293	chr5_MATERNAL	46663888	51824561	centrolign
chr6_PATERNAL	58484657	63761325	chr6_MATERNAL	58406808	64147747	centrolign
chr7_PATERNAL	58365745	62990016	chr7_MATERNAL	58245110	63833413	centrolign
chr8_PATERNAL	44142979	46857663	chr8_MATERNAL	43875938	46736605	centrolign
chr9_PATERNAL	43097467	57492248	chr9_MATERNAL	45003882	69801493	centrolign
```


### Get blocks that have to mapped with minimap2
```
bedtools subtract -a ../../assembly/hg002v1.1.bed -b hg002v1.1.cenSatv2.0.only_hor_hsat.gt_10k.merged_1mb.gt_1mb.adjusted_by_rdna.added_acro_parm.no_sex.bed | \
    grep -v -e "chrEBV" -e "chrM" -e "chrY" -e "chrX" | \
    bedtools subtract -a - -b hg002v1.1.cenSatv2.0.rDNA_plus_GAP.bed > minimap2_blocks.bed
```


### Separate pat and mat blocks
```
cat minimap2_blocks.bed | grep PAT | bedtools sort -i - > minimap2_blocks.pat.bed

cat minimap2_blocks.bed | grep MAT | bedtools sort -i - > minimap2_blocks.mat.bed
```

### Pair  blocks that have to mapped with minimap2

```
paste minimap2_blocks.pat.bed minimap2_blocks.mat.bed | awk '{print $0"\tminimap2"}'  > minimap2_blocks.pairs.bed



cat minimap2_blocks.pairs.bed

chr10_PATERNAL	0	38578643	chr10_MATERNAL	0	38577130	minimap2
chr10_PATERNAL	44232421	135711605	chr10_MATERNAL	43635884	135899001	minimap2
chr11_PATERNAL	0	50964481	chr11_MATERNAL	0	50938529	minimap2
chr11_PATERNAL	53508439	133990155	chr11_MATERNAL	54665772	135235771	minimap2
chr12_PATERNAL	0	34580373	chr12_MATERNAL	0	34545968	minimap2
chr12_PATERNAL	37782453	133573713	chr12_MATERNAL	37778542	133580557	minimap2
chr13_PATERNAL	8564979	9679403	chr13_MATERNAL	12601132	13765549	minimap2
chr13_PATERNAL	13122311	109226264	chr13_MATERNAL	17236289	113336972	minimap2
chr14_PATERNAL	10495661	11578236	chr14_MATERNAL	12419632	13442573	minimap2
chr14_PATERNAL	16762214	105804016	chr14_MATERNAL	18242469	108656801	minimap2
chr15_PATERNAL	15332197	96256967	chr15_MATERNAL	18803957	100881320	minimap2
chr16_PATERNAL	0	34686993	chr16_MATERNAL	0	35888014	minimap2
chr16_PATERNAL	37402212	38566091	chr16_MATERNAL	38204872	39373481	minimap2
chr16_PATERNAL	49744461	93613879	chr16_MATERNAL	46509338	90398409	minimap2
chr17_PATERNAL	0	21837775	chr17_MATERNAL	0	21673042	minimap2
chr17_PATERNAL	28088572	84843407	chr17_MATERNAL	27285530	83770182	minimap2
chr18_PATERNAL	0	15476428	chr18_MATERNAL	0	15477187	minimap2
chr18_PATERNAL	21348472	80778407	chr18_MATERNAL	19676097	78908764	minimap2
chr19_PATERNAL	0	24653975	chr19_MATERNAL	0	24693385	minimap2
chr19_PATERNAL	29564703	61355724	chr19_MATERNAL	29421838	61317335	minimap2
chr1_PATERNAL	0	121897879	chr1_MATERNAL	0	121793663	minimap2
chr1_PATERNAL	146932593	252060741	chr1_MATERNAL	137240935	244022067	minimap2
chr20_PATERNAL	0	26365907	chr20_MATERNAL	0	26299400	minimap2
chr20_PATERNAL	33807407	67035379	chr20_MATERNAL	32960966	66071306	minimap2
chr21_PATERNAL	7284755	8402826	chr21_MATERNAL	10929182	12106868	minimap2
chr21_PATERNAL	10490173	44312598	chr21_MATERNAL	13506641	47311724	minimap2
chr22_PATERNAL	15001301	49430517	chr22_MATERNAL	18892213	53395367	minimap2
chr2_PATERNAL	0	88976293	chr2_MATERNAL	0	89056701	minimap2
chr2_PATERNAL	94232987	241873532	chr2_MATERNAL	94465377	242114397	minimap2
chr3_PATERNAL	0	90555532	chr3_MATERNAL	0	90461937	minimap2
chr3_PATERNAL	96708207	201513538	chr3_MATERNAL	96202232	201097962	minimap2
chr4_PATERNAL	0	49182798	chr4_MATERNAL	0	49229606	minimap2
chr4_PATERNAL	54080686	192384017	chr4_MATERNAL	53502364	191669995	minimap2
chr5_PATERNAL	0	46686377	chr5_MATERNAL	0	46663888	minimap2
chr5_PATERNAL	57738293	188875552	chr5_MATERNAL	51824561	183262120	minimap2
chr6_PATERNAL	0	58484657	chr6_MATERNAL	0	58406808	minimap2
chr6_PATERNAL	63761325	174411882	chr6_MATERNAL	64147747	174742601	minimap2
chr7_PATERNAL	0	58365745	chr7_MATERNAL	0	58245110	minimap2
chr7_PATERNAL	62990016	160104414	chr7_MATERNAL	63833413	160955030	minimap2
chr8_PATERNAL	0	44142979	chr8_MATERNAL	0	43875938	minimap2
chr8_PATERNAL	46857663	146786823	chr8_MATERNAL	46736605	146740677	minimap2
chr9_PATERNAL	0	43097467	chr9_MATERNAL	0	45003882	minimap2
chr9_PATERNAL	57492248	131519435	chr9_MATERNAL	69801493	143803229	minimap2
```

### Put all block pairs in a single file

```
cat centrolign_blocks.pairs.bed minimap2_blocks.pairs.bed | bedtools sort -i - > all_block_pairs.bed
```
