#!/bin/bash
# This script runs the whole genome scan pipeline. For certain scripts you will have to change the populations in the arguments depending on what population contrast you are using 

# source home profile
source $HOME/.bash_profile

# activate environment with all the required tools
conda activate selection_scan_environment

# get AC tables for one suspected allopolyploid population and one pure lyrata population
gatk VariantsToTable \
        -V individual_population_files/BZD_individuals.args.vcf \
        -F CHROM \
        -F POS \
        -F AC \
        -F AN \
        -F DP \
        -O BZD_raw.table

gatk VariantsToTable \
        -V individual_population_files/LIC_individuals.args.vcf \
        -F CHROM \
        -F POS \
        -F AC \
        -F AN \
        -F DP \
        -O LIC_raw.table

echo 'raw tables produced!'

# run G1 script. populations need to be changed in the arguments
# use the allele count tables generated in the previous step to create plots of different metrics like histogram of Fst. An example directory which these will be stored in is the plots are in a folder stores in the directory e.g LICBZD/graphs
# we set the number of SNPs in each analysed window to be 2. It is a crude way of a doing things because there is a lot of random noise that will be picked up using only 2 SNPs but we couldnt find genes to do with meiosis using reasonable genes 
# the -cut sets it so that we are looking for top 0.5% of empirical outliers. This might help a bit with the increased noise caused by the small SNPs per window
# In terms of window size for now the theory is that a window is considered as the base pairs (bp) between two SNPs, but the window size cant exceed 26560 bp. 
# So if we have a SNP at pos 20 and a SNP at pos 40, the window size is 20 and it is less than 26560 so that window is considered in the metrics calculations
python3 G1_outliers.py \
        -i ./ \
        -coh1 LIC \
        -coh2 BZD \
        -snps 2 \
        -cut 0.5

echo 'done with G1, starting G2'

# edit gff to have NW scaffold name instead of scaffold_1
# we need this step because the G2 and G3 script code recognises the scaffold name to be NW_003302555.1 and not scaffold_1
sed 's/scaffold_1/NW_003302555.1/' LyV2.gff > LyV2_edited.gff

echo 'edited gff with success'

# run G2 script. populations need to be changed in the argyments
# this script takes the directory with the files produced from the G1 script as input. The -an argument is the lyrata gff we have been given.
# the gf script is a list of lyrata genes and their homologs in thaliana (with their functions as well)
# the most important file is a file with ALL in the name like this 'LICBZD_2SNPs_5000ppm_ALL_0ol_GF.txt'. It contains all the genes caught in the selection scan name. 
# each line has the name of the gene in lyrata, the homo/ortholog in thaliana, and its function
python3 G2_genes.py -i LICBZD/ -an LyV2_edited.gff -gf LyV2_TAIR10orth_des_20150927.txt

echo 'Done with G2 script'

# create outfile to be used by R script
# the outfile is essentially just a filtered version of the gff that is compatible with the R script in G3 genes 
# it contains the lyrate gene names (AL....), scaffold number, start position, stop position, and strand direction
python3 gff_to_out_file.py -g LyV2.gff -o final_lyrata_genes.out 

echo 'new out file created'

# run the G3 script. Populations need to be changed in the arguments 
# this is the final script in the pipeline. The -suffix, -out, -snps, and -ovlp arguments. The -snps argument has to be the same as the SNPs per window set in the G1 outliers script
# it takes the .out file we generated in the last command as input, as well as the names of the cohorts we are using
# the expected output from this is a allele frequence difference per SNP (AFDsnp),FST,DD, and Dxy plot for candidate genes that have to do with meiosis (it has been hard coded that if a meiosis gene is not caught in the selection scan then no plots will be produced)
# all the plots have a red arrow highlighting the candidate gene, but the most important of the plots is is the AFDsnp plot, why?
# the plots are in a folder stores in the directory e.g LICBZD/genes/directory_with_graphs_as_pdf
python3 G3_graphs.py -coh1 LIC -coh2 BZD -geneor final_lyrata_genes.out -suffix _Aa -out ALL -snps 2 -ovlp 0 -win 30000

echo 'done with everything'
