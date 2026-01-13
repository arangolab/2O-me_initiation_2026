import os
from collections import defaultdict

# Set the directory containing the BED files
bed_dir = "/projects/b1042/Arangolab/2ometh/5utrSequence/nmPositionsGenome/allCodons"
output_gene_file = "pos1AnalyzedGenes.txt"
output_codon_file = "pos1CodonFreq.txt"

# Use a set to store unique gene names
unique_genes = set()
# Dictionary to store codon frequencies
codon_counts = defaultdict(int)

# Loop through all .bed files in the directory
for filename in os.listdir(bed_dir):
    if filename.endswith("Codons.bed"):
        filepath = os.path.join(bed_dir, filename)
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                # Skip section headers like "+1 Nm"
                if line.startswith('+'):
                    continue
                fields = line.split('\t')
                if len(fields) == 6:
                    gene_name = fields[2]
                    codon = fields[5]
                    unique_genes.add(gene_name)
                    codon_counts[codon] += 1

# Write unique gene names
with open(output_gene_file, 'w') as out:
    for gene in sorted(unique_genes):
        out.write(gene + '\n')

# Write codon frequencies
with open(output_codon_file, 'w') as out:
    for codon, count in sorted(codon_counts.items()):
        out.write(f"{codon}\t{count}\n")

print(f"Finished processing.")
print(f"Unique genes written to: {output_gene_file}")
print(f"Codon frequencies written to: {output_codon_file}")
