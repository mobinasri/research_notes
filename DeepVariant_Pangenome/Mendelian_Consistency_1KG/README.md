# Mendelian consistency analysis on 1KG trios

To measure the improvement of pangenome-aware DV over linear-reference-based DV on non-GIAB samples, we ran both methods on 20 trios from the 1000 Genomes Project (1KG). 
Trios that shared any sample with the HPRC-v1.1 panel were excluded to avoid evaluating with the haplotypes that are present in the pangenome. 
These trios were selected randomly from all 1KG trios; 10 trios have male child and 10 have female. 
We then counted variants inconsistent with Mendelian inheritance between the parents and the child. 
For trios with female child we excluded chromosome Y and for the male ones we excluded chromosomes X and Y. 
All variants with at least one non-reference allele in a sample per trio were included for the analysis. 
For example it might be a hom-ref variant (gt=0/0) in child but a non-hom-ref call in a parent, which would be counted as an error.

For all 20 trios, number of Mendelian errors decreased when pangenome-aware DV was used with an average of 5.24% and 7.9% reduction in SNP and InDel errors over linear-ref-based DV. 
The largest improvement was observed the HG02584's trio with 8.79% and 12.69% reduction in the number of SNP and InDel errors.
We also computed Medelian error rates by diving error counts by the total number of variants.
The average improvement for Mendelian error rates were 6.17% and 9.53% for SNPs and InDels respectively and the largest rates were observed in HG02584's trio with 9.55% and 13.93% reduction in SNP and InDel error rates.

## Analysis steps

Below we provide a summary of all the steps we performed in [the jupyter notebook that was used for this analysis.](https://github.com/mobinasri/research_notes/blob/main/DeepVariant_Pangenome/Mendelian_Consistency_1KG/scripts/1KG_Trio_Analysis.ipynb)

#### Selecting 20 random 1KG trios

The information for all pedigrees in 1KG are available through this link:
https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/1kGP.3202_samples.pedigree_info.txt

We downloaded the table and excluded the trios with missing parents or that have any shared sample with the HPRC-v1.1 panel. Next we randomly selected 10 trios with male and 10 trios with female child.
In total there are 60 1KG samples for which we need linear-ref-based and pangenome-aware DV calls. 

##### Male trios
| index | sampleID | fatherID | motherID | sex |
| ----- | -------- | -------- | -------- | --- |
| 0     | HG02686  | HG02684  | HG02685  | 1   |
| 1     | HG01609  | HG01608  | HG01607  | 1   |
| 2     | NA10865  | NA11891  | NA11892  | 1   |
| 3     | HG02584  | HG02582  | HG02583  | 1   |
| 4     | HG02838  | HG02836  | HG02837  | 1   |
| 5     | NA12707  | NA12716  | NA12717  | 1   |
| 6     | HG04215  | HG03679  | HG03642  | 1   |
| 7     | HG02815  | HG02813  | HG02814  | 1   |
| 8     | HG03942  | HG03941  | HG03940  | 1   |
| 9     | HG01579  | HG01577  | HG01578  | 1   |

##### Female trios
| index | sampleID | fatherID | motherID | sex |
| ----- | -------- | -------- | -------- | --- |
| 0     | HG01343  | HG01341  | HG01342  | 2   |
| 1     | NA19659  | NA19658  | NA19657  | 2   |
| 2     | HG03230  | HG03228  | HG03229  | 2   |
| 3     | HG02495  | HG02493  | HG02494  | 2   |
| 4     | NA18484  | NA18486  | NA18488  | 2   |
| 5     | HG02785  | HG02783  | HG02784  | 2   |
| 6     | HG03834  | HG03833  | HG03832  | 2   |
| 7     | HG00700  | HG00698  | HG00699  | 2   |
| 8     | HG01949  | HG01947  | HG01948  | 2   |
| 9     | HG03719  | HG03725  | HG03722  | 2   |


### Download DV call sets for 60 1KG samples (20 trios) 

The DV call sets for all 3,202 1KG samples are available in the following links. For each chromosome there is a separate VCF file.

```
# for pangenome-aware DV
gs://brain-genomics-public/research/cohort/1KGP/vg/pangenome_aware_to_grch38/glnexus_output

# for linear-ref-based DV
gs://brain-genomics-public/research/cohort/1KGP/vg/dv_grch38
```

Since the VCF files are large (e.g., ~100 GB for chr1), we downloaded them in small batches and extracted only the 60 samples required for each chromosome.
We used the scripts `scripts/run_filter_hom_ref.linear_ref_based.bash` and `scripts/run_filter_hom_ref.pangenome_aware.bash`, which download four VCF files at a time, 
create subsetted VCF files containing the 60 samples, and then remove the original large VCF files. 
Subsetting was performed using the Python script `scripts/filter_hom_ref.py`, which retains only variants with at least one non-reference allele among the 60 samples.


#### Note:
The names of the samples for different chromosomes are different for pangenome-aware DV VCF files though they are equivalent. 
For example `ERR3239276_NA06985` for chr22 is equivalent to `NA06985` for chr1. 
For chr1 to 11 samples don't have that ERR prefix but for the rest of the chromosomes some samples like NA06985 have that. 
As mentioned in the jupyter notebook the names are standardized by removing that "ERRXXXX_" before doing the analysis.


### Computing Mendelian consistency
We parsed the 60-sample vcf files in the notebook using cyvcf2 library and iterated over the variants with at least 
one non-ref allele among the samples in each trio. The code is available in this notebook:
https://github.com/mobinasri/research_notes/blob/main/DeepVariant_Pangenome/Mendelian_Consistency_1KG/scripts/1KG_Trio_Analysis.ipynb

We created two tables with the Menedelian stats; one for SNPs and one for InDels, which are available in the following link:

https://docs.google.com/spreadsheets/d/1KDIC3qUzqsCc368-qfyyFdjaAKZPQCAFmJYXZ98SuGg/edit?gid=367147034#gid=367147034
