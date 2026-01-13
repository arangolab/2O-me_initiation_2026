import argparse
import sys

# Function to parse command line arguments
def check_arg(args=None):
    parser = argparse.ArgumentParser(description='Write pattern count given pattern and text')
    parser.add_argument('-i', '--input',
                        help='path to input file',
                        required=True
                        )
    parser.add_argument('-o', '--output',
                        help='output file name',
                        required=True
                        )
    return parser.parse_args(args)

# Retrieve command line arguments and assign to variables
args = check_arg(sys.argv[1:])
infile = args.input
outfile = args.output

# Define the near-cognate codons
near_cognate_codons = {"ATG", "GTG", "TTG", "CTG", "ATA", "ATC", "ATT", "AGG", "ACG", "AAG"}

# Function to extract rows with near-cognate codons and write with matched codon
def extract_near_cognate_codons(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        sections = {
            "-3 Nm": (8, 11, []),
            "-2 Nm": (7, 10, []),
            "-1 Nm": (6, 9, []),
            "+1 Nm": (5, 8, []),
            "+2 Nm": (4, 7, []),
            "+3 Nm": (3, 6, []),
            "+4 Nm": (2, 5, [])
        }

        for line in infile:
            columns = line.strip().split('\t')
            if len(columns) < 5:
                continue  # Skip malformed lines
            seq = columns[4]

            for label, (start, end, collector) in sections.items():
                if end <= len(seq):
                    codon = seq[start:end]
                    if codon in near_cognate_codons:
                        new_line = line.strip() + '\t' + codon + '\n'
                        collector.append(new_line)

        for label, (_, _, lines) in sections.items():
            outfile.write(label + '\n')
            for item in lines:
                outfile.write(item)

# Run the function
extract_near_cognate_codons(infile, outfile)

print(f"Rows with near-cognate codons have been extracted to {outfile}")