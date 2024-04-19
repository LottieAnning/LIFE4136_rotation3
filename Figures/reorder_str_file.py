# This script reorders the population in the .str file so the faststructure output is better to look at

# import modules
import os
import argparse

my_vars = argparse.ArgumentParser(description='Taking input parameters')

my_vars.add_argument('-i', '--input_file',dest='file_one', type=str, metavar='input file 1', required=True, help='give the path to input .str file')
my_vars.add_argument('-o', '--output_file',dest='out_file', type=str, metavar='output file', required=True, help='give the path to the output file')
my_vars.add_argument('-p', '--pops',dest='populations',type=str, metavar='ordered populations', required=False, help='give the populations as a string separated by spaces and in an order you want your faststructure output to be in')
my_vars.add_argument('--str_dir',dest='individual_str_directory', type=str, metavar='output directory', required=True, help='give the path to the output directory for all the individual .str files')




input_args = my_vars.parse_args()

# populations taken in as input 'KEH','OCH','FRE','ROK','HAB','BZD','GYE','KAG','MAU','JOH','MOD','LIC','PEK'
# also need the str file to reorder it

# Take input as a list of popuations in the order you want ['KEH','BZD','OCH','FRE','ROK','HAB','KAG']
pop_order = input_args.populations.split(',')

# store path to input str file in a variable '~/fastSTRUCTURE_scripts_and_guidance-20240314/faststructure_output/vcf_dir/vcf_to_str/final_tetraploids.StructureInputDiploidized.str'
str_file = input_args.file_one

# store the value for the output file '~/fastSTRUCTURE_scripts_and_guidance-20240314/faststructure_output/vcf_dir/vcf_to_str/final_tetraploids.StructureInputDiploidized2.str'
output_file = input_args.out_file

for pop in pop_order:
    os.system(f"grep {pop} {str_file} > {input_args.individual_str_directory}/{pop}.str")
    os.system(f"cat {input_args.individual_str_directory}/{pop}.str >> {output_file}")



