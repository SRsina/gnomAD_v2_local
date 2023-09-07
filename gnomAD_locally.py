#!/usr/bin/env python3


from pysam import VariantFile
import pandas as pd
import os

# genes = pd.read_excel("../../../Downloads/variant_annotation/input/df.xlsx")
genes = pd.read_excel("./input/df.xlsx")
a = len(genes)
print(f'Your data frame has {a} genes')
genes = genes.drop_duplicates(subset=['Gene_symbol'])
b = len(genes)
print(f'Your data frame has {b} unique genes')
genes = genes[genes['Chromosomes'] != 'chrMT']
c = len(genes)
print(
    f'Your data frame has {b - c} genes which their position is in chrMT which our data set does not cover so they are excluded ')
genes = genes.sort_values(by=['Chromosomes', 'Annotation Genomic Range Start', 'Annotation Genomic Range Stop'],
                          ignore_index=True)

# open the in_vcf file with VariantFile
# chr1 test run
# vcf_file = VariantFile(
#     "../../../Downloads/variant_annotation/gnomad_v211_liftover/exome/gnomad.exomes.r2.1.1.sites.liftover_grch38_chr1.vcf.gz")
# for exom files
# vcf_file = VariantFile("../../../Downloads/variant_annotation/gnomad_v211_liftover/exome/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz")
# for genome files
# vcf_file = VariantFile(
#    "../../../Downloads/variant_annotation/gnomad_v211_liftover/genome/gnomad.genomes.r2.1.1.sites.liftover_grch38.vcf.bgz)
# vcf_file = VariantFile("./gnomad/genome/gnomad.genomes.r2.1.1.sites.liftover_grch38.vcf.bgz")
vcf_file = VariantFile("./gnomad/exome/gnomad.exomes.r2.1.1.sites.liftover_grch38_chr1.vcf.gz")
# iterate on the genes
gennum = 0
# chr1 test run
# for i in range(10):
for i in range(len(genes)):
    gennum += 1
    gene_symbol = genes.loc[i, 'Gene_symbol']
    chromosome = genes.loc[i, 'Chromosomes']
    start_position = genes.loc[i, 'Annotation Genomic Range Start']
    end_position = genes.loc[i, 'Annotation Genomic Range Stop']

    # test if outfile is present and skip if it is
    output_file = "./output/exome/" + gene_symbol + '_filtered_variants.vcf'
    if os.path.exists(output_file) or os.path.exists(output_file + '.gz'):
        print(f" {gennum}/{len(genes)} Skipping {gene_symbol} on {chromosome} {start_position} .. {end_position} ")
        continue

    # fetch position from vcf
    records = vcf_file.fetch(chromosome, start_position, end_position)
    # if something found save to vcf outfile
    if records:
        # output_file = "../../../Downloads/variant_annotation/output_vcf/exome/" + gene_symbol + '_filtered_variants.vcf'
        # for genome files
        # output_file = "../../../Downloads/variant_annotation/output_vcf/genome/" + gene_symbol + '_filtered_variants.vcf'
        vcf_extracted = VariantFile(output_file, 'w', header=vcf_file.header)
        for record in records:
            vcf_extracted.write(record)
        vcf_extracted.close()
        print(f" {gennum}/{len(genes)} Found {gene_symbol} on {chromosome} {start_position} .. {end_position}")
# close the vcf file
vcf_file.close()