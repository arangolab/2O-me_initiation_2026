import os
import glob

def parse_block_header(header):
    # Convert block headers like "-3 Nm" or "+2 Nm" into output filenames
    prefix = "allGene_"
    value = header.strip().split()[0]
    if value.startswith("-"):
        return f"{prefix}m{value[1:]}.bed"
    elif value.startswith("+"):
        return f"{prefix}{value[1:]}.bed"
    else:
        raise ValueError(f"Unexpected block header: {header}")

def collect_nm_files(input_pattern):
    output_data = {}

    for filename in glob.glob(input_pattern):
        with open(filename, 'r') as f:
            current_file = None
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.endswith("Nm"):
                    current_file = parse_block_header(line)
                    if current_file not in output_data:
                        output_data[current_file] = []
                else:
                    parts = line.split('\t')
                    if len(parts) < 6:
                        continue  # Skip malformed lines
                    chrom = parts[0]
                    pos = parts[1]
                    strand = parts[3]
                    codon = parts[-1]
                    gene = parts[2]
                    output_data[current_file].append(f"{chrom}\t{pos}\t{pos}\t{codon}\t{strand}\t{gene}")

    # Write each output file
    for out_file, lines in output_data.items():
        with open(out_file, 'a') as f:
            for line in lines:
                f.write(line + '\n')

if __name__ == "__main__":
    # Change the pattern below if your files are in a subdirectory or have a unique suffix
    collect_nm_files("*Cognates.bed")
