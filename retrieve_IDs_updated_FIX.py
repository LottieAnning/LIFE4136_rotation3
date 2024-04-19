# Script to extract individuals and populations from vcf file to be used in faststructure, as well as just individuals to be used in GATK by producing a .args file
import re
import argparse
import os

# Example run command
#python3 ~/lyrata_arenosa_scripts/retrieve_IDs_updated_FIX.py \
#	-i ~/lyrata_arenosa_data/Chrom_1_noSnakemake.lyrata.bipassed.dp.m.bt.1pct.ld_pruned.vcf \
#	--pop 'KEH','BZD','OCH','FRE','ROK','HAB','KAG','KAR','HLY','NKM','PIZ' \
#	--same_file 'yes' \
#	-odir ~/lyrata_arenosa_data \
#	-opre final_tetraploids \
#	--fast_struc 'yes' \
#	-xcl ~/lyrata_arenosa_data/final_tetraploids_exclusions.args

my_vars = argparse.ArgumentParser(description='Taking a vcf file')

my_vars.add_argument('-i', type=str, required=True, help='give the path to tsv')
my_vars.add_argument('-p','--pop',dest='p', type=str, required=True, help='List format is pop1,pop2,pop3')
my_vars.add_argument('--same_file', type=str,dest='same_file', required=True, help='Determines if all output is written to the same file',default='no')
my_vars.add_argument('-odir',type=str, required=True, help='Give the output file directory')
my_vars.add_argument('-opre',type=str, required=False, help='Give the output file prefix when writing output. Needs to be given if same file flag or fast struc file is yes')
my_vars.add_argument('--fast_struc',type=str, required=False,default = 'no', help='If flagged as yes gives an output file that is compatible with the output of faststructure')
my_vars.add_argument('-xcl',type=str, required=False, help='Samples from a given population to exclude',default='no_exclusions')

vcf_file = my_vars.parse_args()

# if they have provided a file to exclude stuff
if vcf_file.xcl != 'no_exclusions':

# open file with samples to exclude
    xcl_file = open(vcf_file.xcl,'r')

    # put samples to exclude into a list
    samples_to_exclude = xcl_file.readlines()
    samples_to_exclude = [samp.rstrip() for samp in samples_to_exclude]

# set sample to exlucde to be an empty string if no xcl file is given
elif vcf_file.xcl == 'no_exclusions':
    samples_to_exclude = []

# Open main file
with open(vcf_file.i,'r+') as f:

	# Get the positions for the line in the vcf with headers
    vcf_headers_pos = re.search('#CHROM.*',f.read())
	
	# Go back to top of file
    f.seek(0)

	# Get list of vcf headers
    vcf_headers = f.read()[vcf_headers_pos.start():vcf_headers_pos.end()]
	
    # Do something different depending on if they want the output written to the same file or not
    if vcf_file.same_file == 'no':

        
        # create a directory to store all .args output files 
        os.makedirs(os.path.join(vcf_file.odir,'gatk_args_output'),exist_ok=True)
	
        # Go through each unique population name
        for pop_name in vcf_file.p.split(','):

            # Make a new directory for each population
            #os.mkdir(os.path.join(vcf_file.odir,pop_name))

            # Open a new file using the population name
            new_file = open(f'{vcf_file.odir}/gatk_args_output/{pop_name}_individuals.args','w+')

            # Open a new file for faststructure compatible output as well
            individual_population_file = open(f'{vcf_file.odir}/gatk_args_output/{pop_name}_individuals_populations.txt','w+')

            # Look in vcf header line for populations with that name
            for header in vcf_headers.split():

                # Write out the IDs into separate text files	
                if header.startswith(pop_name) and header.rstrip() not in samples_to_exclude:
                    
                    new_file.write(f'{header.rstrip()}\n')

                    individual_population_file.write(f'{header.rstrip()}\t{pop_name}\n')
                    
            new_file.close()
            individual_population_file.close()

    # If they want all names in the same file then use a different extension in the output file and write output one line at a time
    elif vcf_file.same_file == 'yes':
        
        # Make directory to store the .args file needed for gatk
        os.makedirs(os.path.join(vcf_file.odir,'gatk_args_output'),exist_ok=True)
        
        # Open a new file using the population name
        new_file = open(f'{vcf_file.odir}/gatk_args_output/{vcf_file.opre}.args','w+')

        # Go through each unique population name
        for pop_name in vcf_file.p.split(','):
                
            # Look in vcf header line for populations with that name
            for header in vcf_headers.split():
                
                # Write out the IDs into separate text files	
                if header.startswith(pop_name) and header.rstrip() not in samples_to_exclude:
                    
                    new_file.write(f'{header.rstrip()}\n')

        new_file.close()

    # create a popname file to use		
    if vcf_file.fast_struc == 'yes':

        # Make a directory to store the files taken by faststructure
        os.makedirs(os.path.join(vcf_file.odir,'faststructure_files'),exist_ok=True)        

        # Open a new file using the population name
        faststruc_file = open(f'{vcf_file.odir}/faststructure_files/{vcf_file.opre}_faststructure_popnames.txt','w+')

        # Go through each unique population name
        for pop_name in vcf_file.p.split(','):

            # Look in vcf header line for populations with that name
            for header in vcf_headers.split():

                # Write out the IDs into separate text files
                if header.startswith(pop_name) and header.rstrip() not in samples_to_exclude:

                    faststruc_file.write(f'{header.rstrip()}\t{pop_name.rstrip()}\n')

        faststruc_file.close()


