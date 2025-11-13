# Using HPRC pangenome v2 for mapping and variant calling with pangenome-aware DV

## Using Adam's BAM file made with `sep8` graph
### The `v2-sep8` graph
Glenn Hickey fixed some of the big snarl issues in the initial HPRC graph-v2 (with the name `ec1M`) and made a new one called `sep8`. We should use the `-eval` version of the graph since it does not have HG002 sample. 
Here is the path where this graph is available in the UCSC Pheonix cluster:

```
/private/home/ghickey/dev/work/chrom-ordering/hprc-sep8-mc-grch38-eval/hprc-sep8-mc-grch38-eval.gbz
/private/home/ghickey/dev/work/chrom-ordering/hprc-sep8-mc-grch38-eval/hprc-sep8-mc-grch38-eval.hapl
```

They are also available in the HPRC S3 bucket:
```
https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/scratch/2025_09_08_minigraph_cactus/hprc-sep8-mc-grch38-eval.gbz
https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/scratch/2025_09_08_minigraph_cactus/hprc-sep8-mc-grch38-eval.hapl
```

### Removing CHM13 and haplotype-sampling

To use the graph for pangenome-aware DV we have to remove the CHM13 sample from the haplotypes and then haplotype-sample it to 32 haplotypes. I added this line to `remove_chm13_vg_1.66/data_table.csv`
```
sample_id,gbz
hprc-sep8-mc-grch38-eval,/private/home/ghickey/dev/work/chrom-ordering/hprc-sep8-mc-grch38-eval/hprc-sep8-mc-grch38-eval.gbz
```
and used this `input_mapping.csv`:
```
input,type,value
RemoveSampleFromGraph.IN_GBZ_FILE,scalar,$input.gbz
RemoveSampleFromGraph.CREATE_INDEX_OPTIONS,scalar,"--snarl-limit 1"
RemoveSampleFromGraph.SAMPLE_NAME_TO_REMOVE,scalar,"CHM13"
RemoveSampleFromGraph.CORES,scalar,16
RemoveSampleFromGraph.DOCKER_IMAGE,scalar,"quay.io/vgteam/vg:v1.66.0"
```

I made the inputs json file with `remove_chm13_vg_1.66/run_remove_sample_from_gbz.sh` and saved it in `remove_sample_from_gbz_input_jsons/hprc-sep8-mc-grch38-eval_remove_sample_from_gbz.json`
```
{
  "RemoveSampleFromGraph.IN_GBZ_FILE": "/private/home/ghickey/dev/work/chrom-ordering/hprc-sep8-mc-grch38-eval/hprc-sep8-mc-grch38-eval.gbz",
  "RemoveSampleFromGraph.CREATE_INDEX_OPTIONS": "--snarl-limit 1",
  "RemoveSampleFromGraph.SAMPLE_NAME_TO_REMOVE": "CHM13",
  "RemoveSampleFromGraph.CORES": 16,
  "RemoveSampleFromGraph.DOCKER_IMAGE": "quay.io/vgteam/vg:v1.66.0"
}
```

Next I did haplotype sampling for which I made `test_samples_vg_1.66/data_table.csv`
```
sample_id,reads
HG002.novaseq.pcr-free,"['https://storage.googleapis.com/brain-genomics-public/research/pangenome_aware_dv_paper/jan_2025/fastq/HG002.novaseq.pcr-free.35x.R1.fastq.gz','https://storage.googleapis.com/brain-genomics-public/research/pangenome_aware_dv_paper/jan_2025/fastq/HG002.novaseq.pcr-free.35x.R2.fastq.gz']"
```
and `test_samples_vg_1.66/input_mapping_graph_v2_sep28_GRCh38_eval.csv`
```
input,type,value
HaplotypeSampling.IN_OUTPUT_NAME_PREFIX,scalar,$input.sample_id
HaplotypeSampling.DIPLOID,scalar,"false"
HaplotypeSampling.INPUT_READ_FILE_ARRAY,array,$input.reads
HaplotypeSampling.HAPLOTYPE_NUMBER_ARRAY,array,"[16, 32]"
HaplotypeSampling.REFERENCE_FASTA,scalar,"/private/groups/patenlab/masri/common/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"
HaplotypeSampling.HAPL_FILE,scalar,"/private/groups/patenlab/masri/haplotype_sampling/graph_v2.0/remove_chm13_vg_1.66/runs_toil_slurm/hprc-sep8-mc-grch38-eval/analysis/remove_sample_from_gbz_outputs/hprc-sep8-mc-grch38-eval.CHM13_removed.hapl"
HaplotypeSampling.IN_DIST_FILE,scalar,"/private/groups/patenlab/masri/haplotype_sampling/graph_v2.0/remove_chm13_vg_1.66/runs_toil_slurm/hprc-sep8-mc-grch38-eval/analysis/remove_sample_from_gbz_outputs/hprc-sep8-mc-grch38-eval.CHM13_removed.dist"
HaplotypeSampling.IN_GBZ_FILE,scalar,"/private/groups/patenlab/masri/haplotype_sampling/graph_v2.0/remove_chm13_vg_1.66/runs_toil_slurm/hprc-sep8-mc-grch38-eval/analysis/remove_sample_from_gbz_outputs/hprc-sep8-mc-grch38-eval.CHM13_removed.gbz"
HaplotypeSampling.CORES,scalar,8
HaplotypeSampling.DOCKER_IMAGE,scalar,"quay.io/vgteam/vg:v1.66.0"
```

The output GBZ files are available on the UCSC cluster (one with hap=32 and one with hap=16):
```
/private/groups/patenlab/masri/haplotype_sampling/graph_v2.0/test_samples_vg_1.66/runs_toil_slurm_graph_v2_sep28_GRCh38_eval/HG002.novaseq.pcr-free/analysis/haplotype_sampling_customized_outputs/HG002.novaseq.pcr-free.hap_num_32.gbz
/private/groups/patenlab/masri/haplotype_sampling/graph_v2.0/test_samples_vg_1.66/runs_toil_slurm_graph_v2_sep28_GRCh38_eval/HG002.novaseq.pcr-free/analysis/haplotype_sampling_customized_outputs/HG002.novaseq.pcr-free.hap_num_16.gbz
```

### Running pangenome-aware DV with Adam's BAM file and hap-sampled graph-sep8

Adam Novak mapped HG002 NovaSeq reads to `graph-v2-sep8` with `vg-v1.68` using this json file:
```
at >/private/groups/patenlab/anovak/trash/for-mobin/mobin-alignment-input.json <<EOF
{
  "Giraffe.INPUT_READ_FILE_1": "https://storage.googleapis.com/brain-genomics-public/research/pangenome_aware_dv_paper/jan_2025/fastq/HG002.novaseq.pcr-free.35x.R1.fastq.gz",
  "Giraffe.INPUT_READ_FILE_2": "https://storage.googleapis.com/brain-genomics-public/research/pangenome_aware_dv_paper/jan_2025/fastq/HG002.novaseq.pcr-free.35x.R2.fastq.gz",
  "Giraffe.GBZ_FILE": "/private/home/ghickey/dev/work/chrom-ordering/hprc-sep8-mc-grch38-eval/hprc-sep8-mc-grch38-eval.gbz",
  "Giraffe.HAPL_FILE": "/private/home/ghickey/dev/work/chrom-ordering/hprc-sep8-mc-grch38-eval/hprc-sep8-mc-grch38-eval.hapl",
  "Giraffe.SAMPLE_NAME": "HG002",
  "Giraffe.OUTPUT_SINGLE_BAM": true,
  "Giraffe.HAPLOTYPE_SAMPLING": true,
  "Giraffe.HAPLOTYPE_NUMBER": 32,
  "Giraffe.DIPLOID": true,
  "Giraffe.SET_REFERENCE": "GRCh38",
  "Giraffe.VG_DOCKER": "quay.io/vgteam/vg:v1.68.0"
}
EOF
```

He also shared the Toil commands he ran to submit the job to SLURM (commit: `5fc3c2fe5241046f43bd503232fd624169a3ef11`):
```
toil-wdl-runner '#workflow/github.com/vgteam/vg_wdl/Giraffe:master' \
    --input /private/groups/patenlab/anovak/trash/for-mobin/mobin-alignment-input.json \
    --batchLogsDir /private/groups/patenlab/anovak/trash/for-mobin/logs \
    --jobStore /private/groups/patenlab/anovak/trash/for-mobin/tree \
    --outputFile /private/groups/patenlab/anovak/trash/for-mobin/out.json \
    --outputDirectory /private/groups/patenlab/anovak/trash/for-mobin/out \
    --logFile /private/groups/patenlab/anovak/trash/for-mobin/log2.txt \
    --batchSystem slurm \
    --slurmTime 11:59:59 \
    --caching=False
```

Note: I used `vg-v1.66` for haplotype sampling above, which is different from the version Adam used for mapping.


Next I ran pangenome-aware DV with this json file
```
cd /private/groups/patenlab/masri/haplotype_sampling/map_short_reads_graph_v2_sep8.vg_1.68/HG002_illumina_AdamNovak/run_png_aware_dv/hap_32

cat pangenome_aware_deepvariant.inputs.json
{
	"bam" : "/private/groups/patenlab/masri/haplotype_sampling/map_short_reads_graph_v2_sep8.vg_1.68/HG002_illumina/HG003.novaseq.pcr-free.35x.graph_v2_sep8.hap32_to_dip.vg-1.68.0.bam", 
	"pangenome" : "/private/groups/patenlab/masri/haplotype_sampling/graph_v2.0/test_samples_vg_1.66/runs_toil_slurm_graph_v2_sep28_GRCh38_eval/HG002.novaseq.pcr-free/analysis/haplotype_sampling_customized_outputs/HG002.novaseq.pcr-free.hap_num_32.gbz",
	"model_ckpt" : "/private/groups/patenlab/masri/apps/bash_scripts/png_aware_dv_models/hap_32/checkpoint-179200-0.98890-1",
	"threads" : "64",
	"bin_version" : "pangenome_aware_deepvariant-1.9.0",
	"reference_fasta" : "/private/groups/patenlab/masri/apps/reference/GRCh38_no_alt_analysis_set.fasta",
	"output_prefix" : "HG002.novaseq.pcr-free.graph_v2.0_eval.sep8.vg_1.68.hap_32_head737001992",
	"output_dir" : "/private/groups/patenlab/masri/haplotype_sampling/map_short_reads_graph_v2_sep8.vg_1.68/HG002_illumina/run_png_aware_dv/hap_32",
	"mount_dir" : "/private/groups/patenlab/masri",
	"pangenome_height" : "37"
}
```
