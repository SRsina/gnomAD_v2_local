#!/usr/bin/env python3
from pysam import VariantFile
import pandas as pd
import os

# Importing the dataframe as it was defined in the previous pipline
genes = pd.read_excel("./data/df.xlsx")

# Cleaning the data
a = len(genes)
print(f'Your data frame has {a} genes')
# Make a unique list of your genes
genes = genes.drop_duplicates(subset=['Gene_symbol'])
b = len(genes)
print(f'Your data frame has {b} unique genes')
# Excluding the mitocondrian chromosome as they are not present in gnomad v2.1.1 lift_over hg38
genes = genes[genes['Chromosomes'] != 'chrMT']
c = len(genes)
print(
    f'Your data frame has {b - c} genes which their position is in chrMT which our data set does not cover so they are excluded ')
# Sort the genes to run faster
genes = genes.sort_values(by=['Chromosomes', 'Annotation Genomic Range Start', 'Annotation Genomic Range Stop'],
                          ignore_index=True)


# open the input_vcf file with VariantFile

# If you are working with exome you need to un coment the following code
# vcf_file = VariantFile("./data/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz")

# If you are working with genome you need to un coment the following code
# vcf_file = VariantFile("./data/gnomad.genomes.r2.1.1.sites.liftover_grch38.vcf.bgz)


for i in range(len(genes)):
    gene_symbol = genes.loc[i, 'Gene_symbol']
    chromosome = genes.loc[i, 'Chromosomes']
    start_position = genes.loc[i, 'Annotation Genomic Range Start']
    end_position = genes.loc[i, 'Annotation Genomic Range Stop']

    # If you are working with exome you need to un coment the following code
    # output_file = "./results/exome/" + gene_symbol + '_filtered_variants.vcf
    
    # If you are working with genome you need to un coment the following code
    # output_file = "./results/genome/" + gene_symbol + '_filtered_variants.vcf'

    # test if outfile is present and skip if it is
    if os.path.exists(output_file) or os.path.exists(output_file + '.gz'):
        print(f" {gennum}/{len(genes)} Skipping {gene_symbol} on {chromosome} {start_position} .. {end_position} ")
        continue

    # fetch position from vcf
    records = vcf_file.fetch(chromosome, start_position, end_position)
    # if something found save to vcf outfile
    if records:
        # If you are working with exome you need to un coment the following code
        # output_file = "./results/exome/" + gene_symbol + '_filtered_variants.vcf'
        
        # If you are working with exome you need to un coment the following code
        # output_file = "./results/genome/" + gene_symbol + '_filtered_variants.vcf'
        
        # produce the output
        vcf_extracted = VariantFile(output_file, 'w', header=vcf_file.header)
        for record in records:
            vcf_extracted.write(record)
        vcf_extracted.close()
        print(f" {gennum}/{len(genes)} Found {gene_symbol} on {chromosome} {start_position} .. {end_position}")
# close the vcf file
vcf_file.close()
