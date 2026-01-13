from Bio import SeqIO
from Bio.Seq import Seq
import argparse
import sys

def check_arg(args=None):
    parser = argparse.ArgumentParser(description='Extract flanking genomic sequence around Nm sites')
    parser.add_argument('-i', '--input', help='Input BED file', required=True)
    parser.add_argument('-o', '--output', help='Output file', required=True)
    parser.add_argument('--genome', help='Path to genome FASTA', required=True)
    parser.add_argument('--flank', type=int, default=5, help='Number of bases to include upstream and downstream (default: 5)')
    return parser.parse_args(args)

args = check_arg(sys.argv[1:])

# Load genome FASTA
genome = SeqIO.to_dict(SeqIO.parse(args.genome, "fasta"))

with open(args.output, "w") as out, open(args.input) as bed:
    for line in bed:
        fields = line.strip().split()
        if len(fields) < 6:
            print(f"Skipping malformed BED line: {line.strip()}")
            continue

        chrom, start, end, consequence, gene, strand = fields[:6]
        if chrom not in genome:
            print(f"Chromosome {chrom} not found in genome FASTA.")
            continue
        if consequence != "5UTR":
            continue

        start = int(start)
        center = int(start)
        flank_start = max(0, (center - args.flank - 1))
        flank_end = center + args.flank


        seq = genome[chrom].seq[flank_start:flank_end]
        if strand == '-':
            seq = seq.reverse_complement()

        out.write(f"{chrom}\t{start}\t{gene}\t{strand}\t{seq}\n")
