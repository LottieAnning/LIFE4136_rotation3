#!/bin/bash
## This script creates an .args file and then filters a vcf by the populations you specify
## Whenever you are changing the populations you use you will change the populations in the --pop argument of the retrieve_IDs_updated_FIX.py python script

# activate conda env outside of script
conda activate /shared/apps/conda/bio2

# get sample names to include in the filtered vcf
python3 retrieve_IDs_updated_FIX.py \
	-i Chrom_1_noSnakemake.lyrata.bipassed.dp.m.bt.1pct.ld_pruned.vcf \
	--pop 'BZD','PEK','SCT','TEM','GYE','JOH','KAG','LIC','LOI','MAU','MOD','PIL','SCB','SWB','HAB','ROK','FRE','OCH','KEH' \
	--same_file 'yes' \
	-odir ./ \
	-opre tetraploid.args \

# prodce dictionary file for reference fasta
gatk CreateSequenceDictionary -R lyrata.fasta

# create index file for the fasta file
samtools faidx lyrata.fasta

# index the vcf file
gatk IndexFeatureFile -I Chrom_1_noSnakemake.lyrata.bipassed.dp.m.bt.1pct.ld_pruned.vcf

# filter the vcf to only contain biallelic variants we are interested in
gatk SelectVariants \
	-R lyrata.fasta \
	-V Chrom_1_noSnakemake.lyrata.bipassed.dp.m.bt.1pct.ld_pruned.vcf \
	-sn tetraploid.args \
	--output tetraploids.vcf

# deactivate conda env
conda deactivate