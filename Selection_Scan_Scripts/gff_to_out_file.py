# This script is meant to produce a .out file compatible with the G3 graphs python script in the selection scan pipeline.
# The .out file has the gene name, scaffold number, start position, stop position, and strand direction
# Absolute paths must be given to the file


# import modules
import argparse
import os

my_vars = argparse.ArgumentParser(description='Taking an input and output file')

my_vars.add_argument('-g','--gff',dest='gff_file', type=str, required=True, help='give the absolute path to the gff file query')

my_vars.add_argument('-o','--out',dest='output_file', type=str, metavar='output file', required=True,
                        help='give the absolute path to the output file name with .out at the end')

# my_vars.add_argument('-d','--dir',dest='home_dir', type=str, required=True,
#                         help='give a ~ only')

my_files = my_vars.parse_args()

print(f"Users home directory is {os.path.expanduser('~')}")

# get all the lines that have gene in them but remove miRNA and tRNA genes
os.system(f"grep 'scaffold.*gene' {my_files.gff_file} | grep -v 'miRNA_gene' | grep -v 'tRNA_gene' > {os.path.expanduser('~')}/intermediate_gff_genes.gff")

# open the files we will use
with open(f"{os.path.expanduser('~')}/intermediate_gff_genes.gff",'r+') as gff_file, open(my_files.output_file,'w+') as out_file:
    
    
    # create var for gene names
    my_genes = []
    
    # store every line of data into a variable
    query_data = gff_file.readlines()

    # Go through every line of data and first get all the info we need before writing output
    for line in query_data:
        
        # separate each column of data on the line
        line_data = line.split('\t')

        # view what the line data looks like
        # print(line_data)

        # find the gene name
        gene_info = line_data[8].split(';')

        actual_gene_name = gene_info[0].replace('ID=','')

        # get start and stop positions
        gene_start = line_data[3]

        gene_end = line_data[4]

        # get chromsome number
        scaffold = line_data[0].replace('scaffold_','')

        # get strand direction
        strand_direciton = line_data[6]

        # write output
        out_file.write(f"{actual_gene_name}\t{scaffold}\t{gene_start}\t{gene_end}\t{strand_direciton}\n")

# remove intermediate file
os.system(f"rm {os.path.expanduser('~')}/intermediate_gff_genes.gff")
        
 
