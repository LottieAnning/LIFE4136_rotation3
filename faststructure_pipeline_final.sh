#!/bin/bash

# This script is going to run the whole faststructure pipeline. Running it will:
# 1. Produce a filtered vcf with only specific populations
# 2. Create VCFs for each population in the filtered vcf you created (i.e OCH.vcf,PEK.vcf). This becomes useful later on when trying to find signatures of allopolyploidy
# 3. Convert your filtered vcf into a .str file so it can work with faststructure, and reorders the .str file so the faststructure output is in the order we want.
# 4. Run faststructure on our data and produce csv files for K values 2-4 in a format acceptable by the ommicsspeaks website 

# Before running you will need:
# 1. A directory  (~/faststructure), and within that directory there should be:
    # A. The lyrata.fasta file
    # B. The Chrom_1_noSnakemake.lyrata.bipassed.dp.m.bt.1pct.ld_pruned.vcf file containing all the populations. If you have the compressed .gz file make sure you use the unzip command  
    # C. A file called samples_to_exclude.args that contains all the outlier individuals for each population, so we can exclude them from the filtered vcf

# 2. Scripts:
    # A. retrieve_IDs_updated_FIX.py
    # B. reorder_str_file.py
    # C. Cochlearia_create_structure_file.py
    # D. faststructure_pipeline_final.sh

# There are quite a lot of files and directories created from this pipeline but only a few are important. Important output from this pipeline:
# 1. There should be a directory (~/final_omicsspeaks_output) and it should contain the csv for K values 2-4 which can be put into the ommicsspeaks website
# 2. A directory (~/final_populations) that contains the filtered vcf with all the populations you chose. This can be used to put into the R scripts to make PCA and splitstree output
# 3. A directory (~/individual_population_files) that contains the vcf files for each individual population from the group of populations we settled on. This is useful for calculating allele frequencies


# NOTE
# If you want to change the populations used you will have to change the pops used in line 61, 100 and 161

# Actual code 

### Prepare your environment:

#make a directory for the omicsspeaks csv output
mkdir final_omicsspeaks_output
#make a directory for population specific .str files we produce when reordering the faststructure input
mkdir individual_str_output
#make a directory for all filtered vcf output data
mkdir final_populations
#make directory to store all the output for the individuals
mkdir individual_population_files
#make directories to store faststructure output
mkdir -p faststructure_output/vcf_dir
mkdir faststructure_output/final_svg_files

#get the 'attributes' of your environment 
source $HOME/.bash_profile

## Part1 Getting a filtered vcf. The populations should be given in the order you want the final faststructure output to be in. 
## Whenever you are changing the populations you use you will change the populations in the --pop argument of the retrieve_IDs_updated_FIX.py python script

# activate conda env outside of script
conda activate /shared/apps/conda/bio2

# make directory for all filtered vcf output data
mkdir final_populations

# get sample names to include in the filtered vcf
python3 retrieve_IDs_updated_FIX.py \
	-i Chrom_1_noSnakemake.lyrata.bipassed.dp.m.bt.1pct.ld_pruned.vcf \
	--pop 'KEH','BZD','OCH','FRE','ROK','HAB','KAG','MAU','JOH','MOD','LIC','PEK' \
	--same_file 'yes' \
	-odir final_populations \
	-opre diploids_hybird_1arenosa_1lyrata \
	--fast_struc 'yes' \
	-xcl samples_to_exclude.args 

# prodce dictionary file for reference fasta
gatk CreateSequenceDictionary \
	-R lyrata.fasta

# create index file for the fasta file
samtools faidx lyrata.fasta

# index the vcf file
gatk IndexFeatureFile \
	-I Chrom_1_noSnakemake.lyrata.bipassed.dp.m.bt.1pct.ld_pruned.vcf

# filter the vcf to only contain biallelic variants we are interested in
gatk SelectVariants \
	-R lyrata.fasta \
	-V Chrom_1_noSnakemake.lyrata.bipassed.dp.m.bt.1pct.ld_pruned.vcf \
	-sn final_populations/gatk_args_output/*.args \
	--output final_populations/diploids_hybird_1arenosa_1lyrata.vcf

# deactivate conda env
conda deactivate

## Part2 Producing individual vcf files for each of the populations you have narrowed down to 
# activate conda env outside of script
conda activate /shared/apps/conda/bio2

# make directory to store all the output for the individuals
mkdir individual_population_files

# get sample names to include in the filtered vcf
# Now do everything again but produce files for the individual populations
python3 retrieve_IDs_updated_FIX.py \
        -i Chrom_1_noSnakemake.lyrata.bipassed.dp.m.bt.1pct.ld_pruned.vcf \
        --pop 'KEH','BZD','OCH','FRE','ROK','HAB','KAG','MAU','JOH','MOD','LIC','PEK' \
        --same_file 'no' \
        -odir individual_population_files \
	-xcl samples_to_exclude.args 

# prodce dictionary file for reference fasta
gatk CreateSequenceDictionary \
	-R ~/lyrata_arenosa_data/lyrata.fasta

# create index file for the fasta file
samtools faidx ~/lyrata_arenosa_data/lyrata.fasta

# go through all the arg files to produce a vcf for each individual population
for file in individual_population_files/gatk_args_output/*.args ; do
	
	# get name of just file and not the full path
	individual_file=$(basename "$file")
	
	## filter the vcf to only contain biallelic variants we are interested in
	gatk SelectVariants \
		-R lyrata.fasta \
		-V Chrom_1_noSnakemake.lyrata.bipassed.dp.m.bt.1pct.ld_pruned.vcf \
		-sn "$file" \
		--output individual_population_files/"$individual_file.vcf"
done

# deactivate conda environment
conda deactivate

## Part3 Convert filtered vcf file to a str file
# use python3 in our cloud hpc environment
# has numpy version 1.19.5

# make directories to store faststructure output
mkdir -p faststructure_output/vcf_dir
mkdir faststructure_output/final_svg_files

# remove everything from the directory with vcfs
# rm faststructure_output/vcf_dir/*.vcf

# remove everything from the str files
# rm faststructure_output/vcf_dir/vcf_to_str/*.str

# move file with filtered vcf to the vcf directory
cp final_populations/*.vcf faststructure_output/vcf_dir

# convert polyploids data to format acceptable to fastSTRUCTURE
python3 Cochlearia_create_structure_file.py \
        -v faststructure_output/vcf_dir/ \
        -o filtered_vcf_converted_to_str \
        -s true


## remove first and last line from .str file
sed -i '1d;$d' faststructure_output/vcf_dir/vcf_to_str/filtered_vcf_converted_to_str.StructureInputDiploidized.str

## Part4 reorder the .str file based on what order we've done
 
python3 reorder_str_file.py \
	-i faststructure_output/vcf_dir/vcf_to_str/filtered_vcf_converted_to_str.StructureInputDiploidized.str \
	-o faststructure_output/vcf_dir/vcf_to_str/filtered_vcf_converted_to_str.StructureInputDiploidized2.str \
	-p 'KEH','BZD','OCH','FRE','ROK','HAB','KAG','MAU','JOH','MOD','LIC','PEK' \
	--str_dir individual_str_output

# actiavte conda environment for faststructure
conda activate /shared/conda/faststructure

# move into final svg directory so all output gets put into there
cd faststructure_output/final_svg_files

## Part5 run fasstructure command on our reordered .str file with different K values
for i in {1..4}; do 
	echo $i
	python /shared/conda/faststructure/bin/structure.py \
		-K $i \
		--input ~/fast_structure/faststructure_output/vcf_dir/vcf_to_str/filtered_vcf_converted_to_str.StructureInputDiploidized2 \
		--output ~/fast_structure/faststructure_output/final_svg_files/filtered_populations \
		--format str \
		--full

done

# run choose K
python /shared/conda/faststructure/bin/chooseK.py --input ~/faststructure_output/final_svg_files/filtered_populations

# We do not run distruct because the output from distruct is not good enough to be put into our reports/presentations. We will instead us omicsspeaks

# deactivate conda environment
conda deactivate

## Part6 Run code to create csv file necessary for ommics speaks. 
# for K3
paste -d '\t' fast_structure/final_populations/faststructure_files/diploids_hybird_1arenosa_1lyrata_faststructure_popnames.txt fast_structure/faststructure_output/final_svg_files/filtered_populations.3.meanQ > omics_speaks_K3.tsv

cat omics_speaks_K3.tsv | tr '\t' ',' | tr -s '[:blank:]' ',' > omics_speaks_K3.csv

# for K2
paste -d '\t' fast_structure/final_populations/faststructure_files/diploids_hybird_1arenosa_1lyrata_faststructure_popnames.txt fast_structure/faststructure_output/final_svg_files/filtered_populations.2.meanQ > omics_speaks_K2.tsv

cat omics_speaks_K2.tsv | tr '\t' ',' | tr -s '[:blank:]' ',' > omics_speaks_K2.csv

# for K4
paste -d '\t' fast_structure/final_populations/faststructure_files/diploids_hybird_1arenosa_1lyrata_faststructure_popnames.txt fast_structure/faststructure_output/final_svg_files/filtered_populations.4.meanQ > omics_speaks_K4.tsv

cat omics_speaks_K4.tsv | tr '\t' ',' | tr -s '[:blank:]' ',' > omics_speaks_K4.csv
