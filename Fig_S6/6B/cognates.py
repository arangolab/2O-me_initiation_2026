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

# Function to extract rows with near-cognate codons
def extract_near_cognate_codons(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        one = [] # initialize list to store seq's with Nm at +1 position
        two = []
        three = []
        four = []
        minusone = []
        minustwo = []
        minusthree = []
        for line in infile:
            columns = line.strip().split('\t')
            seq = columns[7]
            # Check if any near-cognate codon is a substring of seq and subset them based on Nm location
            if any(codon in seq[5:8] for codon in near_cognate_codons):
                one.append(line)
            if any(codon in seq[4:7] for codon in near_cognate_codons):
                two.append(line)
            if any(codon in seq[3:6] for codon in near_cognate_codons):
                three.append(line)
            if any(codon in seq[2:5] for codon in near_cognate_codons):
                four.append(line)
            if any(codon in seq[6:9] for codon in near_cognate_codons):
                minusone.append(line)
            if any(codon in seq[7:10] for codon in near_cognate_codons):
                minustwo.append(line)
            if any(codon in seq[8:11] for codon in near_cognate_codons):
                minusthree.append(line)
        # Write the contents of the lists to the output file, each item separated by a newline
        outfile.write("-3 Nm" + '\n')
        for item in minusthree:
            outfile.write(item)
        outfile.write("-2 Nm" + '\n')
        for item in minustwo:
            outfile.write(item)
        outfile.write("-1 Nm" + '\n')
        for item in minusone:
            outfile.write(item)
        outfile.write("+1 Nm" + '\n')
        for item in one:
            outfile.write(item)
        outfile.write("+2 Nm" + '\n')
        for item in two:
            outfile.write(item)
        outfile.write("+3 Nm" + '\n')
        for item in three:
            outfile.write(item)
        outfile.write("+4 Nm" + '\n')
        for item in four:
            outfile.write(item)

# Run the function
extract_near_cognate_codons(infile, outfile)

print(f"Rows with near-cognate codons have been extracted to {outfile}")