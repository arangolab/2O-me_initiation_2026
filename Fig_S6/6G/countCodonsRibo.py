from collections import Counter
from Bio import SeqIO
from itertools import product
import csv

# Generate all 64 possible codons
bases = ['A', 'T', 'C', 'G']
all_codons = [''.join(c) for c in product(bases, repeat=3)]

# Initialize codon counter with all codons set to zero
codon_counts = Counter({codon: 0 for codon in all_codons})

# Read FASTA and count codons in a sliding window
for record in SeqIO.parse("/projects/b1042/Arangolab/2ometh/5utrSequence/nmPositionsGenome/allCodons/filteredRiboSeq5utr.fasta", "fasta"):
    seq = str(record.seq).upper()
    for i in range(len(seq) - 2):  # sliding window
        codon = seq[i:i+3]
        if codon in codon_counts:
            codon_counts[codon] += 1

with open("5utrCodonCountsRiboSeq.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["Codon", "Count"])
    for codon in sorted(all_codons):
        writer.writerow([codon, codon_counts[codon]])

