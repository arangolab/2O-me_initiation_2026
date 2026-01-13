import os
from collections import defaultdict

# Set the directory containing the BED files
bed_file = "/projects/b1042/Arangolab/2ometh/5utrSequence/nmPositionsGenome/nearCognateAllNmBoxplot/allGene_1.bed"
NUG_gene_file = "nugNMgenes.txt"
CUG_gene_file = "cugNMgenes.txt"
GUG_gene_file = "gugNMgenes.txt"

# Use a set to store unique gene names
unique_genes = set()
NUG_genes = set()
CUG_genes = set()
GUG_genes = set()

nug_genes_list = ["ATG", "GTG", "CTG"]

# Open the output .bed file for writing
# Loop through all .bed files in the directory
with open(bed_file, 'r') as f:
    for line in f:
        line = line.strip()
        fields = line.split('\t')
        if len(fields) == 6:
            gene = fields[-1]
            unique_genes.add(gene)
            codon = fields[-3]
            if codon in nug_genes_list:
                NUG_genes.add(gene)
            if codon == "GTG":
                GUG_genes.add(gene)
            if codon == "CTG":
                CUG_genes.add(gene)

# Write unique gene names
with open(NUG_gene_file, 'w') as out:
    for gene in sorted(NUG_genes):
        out.write(gene + '\n')

with open(GUG_gene_file, 'w') as out:
    for gene in sorted(GUG_genes):
        out.write(gene + '\n')

with open(CUG_gene_file, 'w') as out:
    for gene in sorted(CUG_genes):
        out.write(gene + '\n')

print("Finished processing.")