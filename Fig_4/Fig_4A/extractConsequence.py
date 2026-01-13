import pandas as pd
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

def condense_rows(input_file, output_file):
    # Read the tab-delimited file with no header
    df = pd.read_csv(input_file, sep='\t', header=None)

    # Assign meaningful column names
    df.columns = ['chrom', 'start', 'end', 'id1', 'id2', 'strand', 'chrom2', 'region_start', 'region_end', 'consequence', 'gene', 'strand2']

    # Sanity check: ensure strand columns match
    assert (df['strand'] == df['strand2']).all(), "Strand mismatch between columns"

    # Group by chromosome, start, end, and strand
    grouped = df.groupby(['chrom', 'start', 'end', 'strand'])

    condensed_rows = []

    for group_key, group_df in grouped:
        chrom, start, end, strand = group_key

        # Most common consequence
        consequence_mode = group_df['consequence'].mode()
        top_consequence = consequence_mode.iloc[0] if not consequence_mode.empty else 'NA'

        # Most common gene name
        gene_mode = group_df['gene'].mode()
        top_gene = gene_mode.iloc[0] if not gene_mode.empty else 'NA'

        # Output: chrom, start, end, consequence, gene, strand
        condensed_rows.append([chrom, start, end, top_consequence, top_gene, strand])

    # Create output DataFrame
    condensed_df = pd.DataFrame(condensed_rows, columns=['chrom', 'start', 'end', 'consequence', 'gene', 'strand'])

    # Save to output file
    condensed_df.to_csv(output_file, sep='\t', index=False, header=False)


condense_rows(infile, outfile)