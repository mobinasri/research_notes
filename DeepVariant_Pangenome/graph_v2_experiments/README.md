## Running pangenome-aware DeepVariant with HPRC graph-v2.0

VG team has been working on generating new graphs using HPRC-Release2 assemblies. Based on the latest internal tests the accuracy of variant calling with the latest and greatest graph-v2.0 is 
at least comparable with graph-v1.1. I'm going to use the eval graph (without HG002 sample) to test pangenome-aware DeepVariant. We have to run haplotype sampling since the number of haplotypes 
in graph-v2.0 is about 460 and we cannot put all of them in the DeepVariant pileup tensors. I generated haplotype-sampled graphs with 32 and 16 haplotypes and used the resulting gbz files for 
mapping HG002 Illumina reads with vg-1.66 and variant calling with pangenome-aware DV-v1.9.0. 
The graphs are available in this directory on the UCSC Phoenix cluster `/private/groups/cgl/hprc-graphs/hprc-v2.0-feb28/snarl-filtering/`. I'm using the GRCh38-based graph for this analysis.

### Removig CHM13 from graph-v2.0

During haplotype sampling the CHM13 reference will remain intact in the graph regardless of the reads we use for sampling. It is not possible to skip it while generating pileup tensors in DeepVariant
so I removed the CHM13 reference from the original graph using the wdl [remove_sample_from_gbz.wdl](https://github.com/mobinasri/vg_wdl/blob/master/workflows/remove_sample_from_gbz.wdl). 

The related SLURM job was submitted on the UCSC Phoenix cluster using the tables and bash script available in the `remove_chm13/` directory.

```
cd /private/groups/patenlab/masri/haplotype_sampling/graph_v2.0/remove_chm13_vg_1.66
bash run_remove_sample_from_gbz.sh 
```

### Haplotype sampling

I ran haplotype sampling with 32 and 16 haplotypes using the wdl [haplotype_sampling_customized.wdl](https://github.com/mobinasri/vg_wdl/blob/master/workflows/haplotype_sampling_customized.wdl).

The related SLURM job was submitted on the UCSC Phoenix cluster using the tables and bash script available in the `remove_chm13/` directory.
```
cd /private/groups/patenlab/masri/haplotype_sampling/graph_v2.0/test_samples_vg_1.66
bash run_haplotype_sampling.sh
```

### Mapping reads with vg-1.66

I used the haplotype-sampled gbz file with 32 haplotypes to map HG002 Illumina reads to GRCh38. I didn't pass the hapl file and kmer kff file since if you pass these files to vg 
it will run haplotype sampling in the diploid mode as an internal step. I wanted to use exactly the same graph (with 32 haplotypes) both for mapping and running pangenome-aware DV.

The related SLURM job was submitted on the UCSC Phoenix cluster using the json file and bash script available in the `map_vg/` directory.
```
cd /private/groups/patenlab/masri/haplotype_sampling/map_short_reads_graph_v2_ec1M.vg_1.66
sbatch map_vg_bam.sh HG002_illumina/inputs.json
```

### Running Pangenome-aware DV 

I used the DeepVariant model trained with graph-v1.1 and 32 haplotypes. The checkpoint files were available in these links:
```
https://storage.googleapis.com/brain-genomics-public/research/pangenome_aware_dv_paper/may_2025/models/pangenome_aware_dv_32_haps/checkpoint-179200-0.98890-1.data-00000-of-00001
https://storage.googleapis.com/brain-genomics-public/research/pangenome_aware_dv_paper/may_2025/models/pangenome_aware_dv_32_haps/checkpoint-179200-0.98890-1.index
https://storage.googleapis.com/brain-genomics-public/research/pangenome_aware_dv_paper/may_2025/models/pangenome_aware_dv_32_haps/example_info.json
```

The related SLURM job was submitted on the UCSC Phoenix cluster using the json file and bash script available in the `pang_aware_dv/hap_32` directory.
```
cd /private/groups/patenlab/masri/haplotype_sampling/map_short_reads_graph_v2_ec1M.vg_1.66/HG002_illumina/run_png_aware_dv/hap_32
sbatch /private/groups/patenlab/masri/apps/bash_scripts/pangenome_aware_deepvariant.bash pangenome_aware_deepvariant.inputs.json
```

### Run Happy

The output vcf file was benchmarked against the T2T-Q100-v1.1 truth set using two high-confidence bed files; HG002-T2T-v1.1 and the intersection of HG002-GIAB-v4.2.1 and HG002-T2T-v1.1
The intersected BED file was created to measure the performance in easier regions.

The related SLURM job was submitted on the UCSC Phoenix cluster using the json file and bash script available in the 
`pang_aware_dv/hap_32/happy/t2t_and_giab_conf` and  `pang_aware_dv/hap_32/happy/t2t_conf` directories.
```
# for T2T-v1.1 bed file
cd /private/groups/patenlab/masri/haplotype_sampling/map_short_reads_graph_v2_ec1M.vg_1.66/HG002_illumina/run_png_aware_dv/hap_32/happy/t2t_conf
sbatch /private/groups/patenlab/masri/apps/bash_scripts/happy.bash happy.inputs.json

# for intersected bed file
cd /private/groups/patenlab/masri/haplotype_sampling/map_short_reads_graph_v2_ec1M.vg_1.66/HG002_illumina/run_png_aware_dv/hap_32/happy/t2t_and_giab_conf
sbatch /private/groups/patenlab/masri/apps/bash_scripts/happy.bash happy.inputs.json
```


### Upload output files to the GCP bucket

#### Haplotype-sampled GBZ files
```
gs://pepper-deepvariant/mobinasri/pangenome_paper/graph_v2_ec1M_vg_1.66/HG002/illumina/gbz_files/HG002.novaseq.pcr-free.chm13_removed.hap_num_16.gbz
gs://pepper-deepvariant/mobinasri/pangenome_paper/graph_v2_ec1M_vg_1.66/HG002/illumina/gbz_files/HG002.novaseq.pcr-free.chm13_removed.hap_num_32.gbz
```

#### Read BAM files (Created with hap32 gbz files)
```
gs://pepper-deepvariant/mobinasri/pangenome_paper/graph_v2_ec1M_vg_1.66/HG002/illumina/read_mapping/HG002.novaseq.pcr-free.graph_v2.0_eval.ec1M.vg_1.66.hap32.bam
gs://pepper-deepvariant/mobinasri/pangenome_paper/graph_v2_ec1M_vg_1.66/HG002/illumina/read_mapping/HG002.novaseq.pcr-free.graph_v2.0_eval.ec1M.vg_1.66.hap32.bam.bai
```

#### VCF files
```
gs://pepper-deepvariant/mobinasri/pangenome_paper/graph_v2_ec1M_vg_1.66/HG002/illumina/pang_aware_dv_hap32/HG002.novaseq.pcr-free.graph_v2.0_eval.ec1M.vg_1.66.hap_32_head737001992.vcf.gz
gs://pepper-deepvariant/mobinasri/pangenome_paper/graph_v2_ec1M_vg_1.66/HG002/illumina/pang_aware_dv_hap32/HG002.novaseq.pcr-free.graph_v2.0_eval.ec1M.vg_1.66.hap_32_head737001992.vcf.gz.tbi
```

#### Happy outputs (T2T and GIAB intersection BED file)
```
gs://pepper-deepvariant/mobinasri/pangenome_paper/graph_v2_ec1M_vg_1.66/HG002/illumina/pang_aware_dv_hap32/t2t_and_giab_conf_happy/HG002.novaseq.pcr-free.graph_v2.0_eval.ec1M.vg_1.66.hap_32_head737001992.summary.csv
gs://pepper-deepvariant/mobinasri/pangenome_paper/graph_v2_ec1M_vg_1.66/HG002/illumina/pang_aware_dv_hap32/t2t_and_giab_conf_happy/HG002.novaseq.pcr-free.graph_v2.0_eval.ec1M.vg_1.66.hap_32_head737001992.vcf.gz
gs://pepper-deepvariant/mobinasri/pangenome_paper/graph_v2_ec1M_vg_1.66/HG002/illumina/pang_aware_dv_hap32/t2t_and_giab_conf_happy/HG002.novaseq.pcr-free.graph_v2.0_eval.ec1M.vg_1.66.hap_32_head737001992.vcf.gz.tbi
```

#### Happy outputs (T2T BED file)
```
gs://pepper-deepvariant/mobinasri/pangenome_paper/graph_v2_ec1M_vg_1.66/HG002/illumina/pang_aware_dv_hap32/t2t_conf_happy/HG002.novaseq.pcr-free.graph_v2.0_eval.ec1M.vg_1.66.hap_32_head737001992.summary.csv
gs://pepper-deepvariant/mobinasri/pangenome_paper/graph_v2_ec1M_vg_1.66/HG002/illumina/pang_aware_dv_hap32/t2t_conf_happy/HG002.novaseq.pcr-free.graph_v2.0_eval.ec1M.vg_1.66.hap_32_head737001992.vcf.gz
gs://pepper-deepvariant/mobinasri/pangenome_paper/graph_v2_ec1M_vg_1.66/HG002/illumina/pang_aware_dv_hap32/t2t_conf_happy/HG002.novaseq.pcr-free.graph_v2.0_eval.ec1M.vg_1.66.hap_32_head737001992.vcf.gz.tbi
```

