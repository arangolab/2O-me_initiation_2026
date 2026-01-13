from collections import defaultdict
from Bio import SeqIO

# Dictionary to store the longest sequence per gene
gene_to_longest_seq = {}

# Read FASTA
for record in SeqIO.parse("riboseq5utrSequences.fasta", "fasta"):
    if record.seq.startswith("Sequence"): # skip entries without sequences present
         continue
    gene_name = record.description.split("|")[-4]
    transcript_length = int(record.description.split("|")[-2])
    if gene_name not in gene_to_longest_seq:
            gene_to_longest_seq[gene_name] = record
    else:
        existing_length = int(gene_to_longest_seq[gene_name].description.split("|")[-2])
        if transcript_length > existing_length:
            gene_to_longest_seq[gene_name] = record

# Write output FASTA
with open("filteredRiboSeq5utr.fasta", "w") as out_f:
    SeqIO.write(gene_to_longest_seq.values(), out_f, "fasta")
