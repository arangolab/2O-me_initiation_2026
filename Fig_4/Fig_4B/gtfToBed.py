def gtf_to_bed_with_consequence(gtf_file, bed_file):
    valid_features = {
        'five_prime_utr': "5UTR",
        'three_prime_utr': "3UTR",
        'CDS': 'CDS'
    }

    with open(gtf_file, 'r') as infile, open(bed_file, 'w') as outfile:
        for line in infile:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            chrom, source, feature, start, end, score, strand, frame, attributes = fields

            # Only keep UTR and CDS features
            if feature not in valid_features:
                continue

            # Only keep protein_coding entries
            if 'gene_type "protein_coding"' not in attributes and 'transcript_type "protein_coding"' not in attributes:
                continue

            # Extract gene_name
            gene_name = 'NA'
            for attr in attributes.strip().split(';'):
                attr = attr.strip()
                if attr.startswith('gene_name'):
                    gene_name = attr.split('"')[1]
                    break

            consequence = valid_features[feature]

            # BED is 0-based, half-open
            bed_start = int(start) - 1
            bed_end = int(end)

            # Output format: chrom, start, end, consequence, gene_name, strand
            outfile.write(f'{chrom}\t{bed_start}\t{bed_end}\t{consequence}\t{gene_name}\t{strand}\n')

gtf_to_bed_with_consequence('/projects/b1042/Arangolab/2ometh/canonicalStartSites/Homo_sapiens.GRCh38.113.chr.gtf', 'genomeConsequence0based.bed')