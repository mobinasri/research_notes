## Create assembly-based truth sets using HPRC-R2 assemblies

The samples HG001-007 are being used frequently for training models for variant calling therefore the models might overfit to these samples and testing on these samples
cannot reveal this bias. To measure how much the models are overfit to these samples we should create some high-quality truth sets using data from samples other than HG001-007. 
I follow the steps mentioned below to create such truth sets:
- Find 5 HPRC-R2 samples with phased diploid assemblies polished with DeepPolisher (using PacBio-HiFi and ONT data).
  Those samples should be absent from HPRC-R1 since we want to use graph-v1.1 for testing pangenome-aware DeepVariant.
- Map assemblies to GRCh38 or CHM13-v2.0 with dipcall (I prefer CHM13-v2.0 but I have to check if DRAGEN works with it)
- Get output VCF and high-confidence BED file
