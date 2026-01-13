import pandas as pd
from collections import Counter
import argparse
import sys

# Function to parse command line arguments
def check_arg(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input',
                        help='path to input file',
                        required=True
                        )
    parser.add_argument('-o', '--output',
                        help='output file name',
                        required=True
                        )
    return parser.parse_args(args)

args = check_arg(sys.argv[1:])
infile = args.input
outfile = args.output

# Adjust column names based on your input structure
column_names = [
    'chrom', 'start', 'end', 'gene',
    'full_chrom', 'full_start', 'full_end',
    'consequence', 'full_gene', 'strand'
]

# Load the data
df = pd.read_csv(infile, sep='\t', header=None, names=column_names)

# Group by the key fields
group_cols = ['chrom', 'start', 'end', 'gene']
grouped = df.groupby(group_cols)

condensed_rows = []

for name, group in grouped:
    most_common_consequence = Counter(group['consequence']).most_common(1)[0][0]
    representative_row = group.iloc[0].copy()

    # Set most common consequence
    representative_row['consequence'] = most_common_consequence

    condensed_rows.append(representative_row)

# Create a condensed DataFrame
condensed_df = pd.DataFrame(condensed_rows)

# Save to output file
condensed_df.to_csv(outfile, sep='\t', index=False)
