#!/usr/bin/env python3
"""
Convert CSV file to BED format.
"""

def read_bed2(filename):
    """Read bed2 file (chr_Position, gene, etc....)"""
    records = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('chr_Position'):
                continue  # Skip header line
            parts = line.split(',')
            if len(parts) >= 6:
                chrom_pos = parts[0].split('_')
                if len(chrom_pos) != 2:
                    continue  # Skip malformed lines
                chrom, start = chrom_pos
                # Skip if start is not a valid number (e.g., header)
                try:
                    start = int(start)
                except ValueError:
                    continue
                chrom = chrom.replace('chr', '')
                end = start + 1  # end position
                records.append({
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'gene': parts[1]
                })
    return records


def write_bed_format(records, output_file):
    """Write records to BED format file"""
    with open(output_file, 'w') as out:
        for record in records:
            out.write(f"{record['chrom']}\t{record['start']}\t{record['end']}\t{record['gene']}\n")
    
    print(f"Wrote {len(records)} records to {output_file}")


def main():
    # Input file
    C42_csv = '../supplementaryTableS4.csv'
    output_file = 'formattedAllC42sites.bed'
    
    print(f"Reading {C42_csv}...")
    bed2_records = read_bed2(C42_csv)
    print(f"  Found {len(bed2_records)} records")
    
    print("Writing BED format...")
    write_bed_format(bed2_records, output_file)


if __name__ == '__main__':
    main()