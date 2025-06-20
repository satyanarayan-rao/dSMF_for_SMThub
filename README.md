# dSMF footprint generation pipeline for SMTHub

SMTHub hosts FootPrint BigBed files generated from NOMe-seq and dSMF data for now. User may generate the bigbed file using this pipeline for dSMF data. Please note that if you have NOMe-seq data, please use the [SMF_for_SMTHub](https://github.com/satyanarayan-rao/SMF_for_SMThub) github repo.

## Create virtual environment and install required tools

**NOTE:** You can skip this step if you have already created this virtual environment while working on [SMF_for_SMTHub](https://github.com/satyanarayan-rao/SMF_for_SMThub) 

```
mamba create -n smthub_production python=3.12
mamba activate smthub_production

mamba install -c bioconda snakemake=7.26 scanf trim-galore bedtools bismark samtools bamtools pandas
```

### Generate the BigBed file

```
snakemake  --snakefile nome_seq_data_to_smf_bigbed.smk suppressed_merged/suppressed_merged_demo_s2_to_dm3_with_wobble_1_min_fp_10_and_mvec.bb --configfile configs/config.yaml -j4

```
