# Figure 4C: Codon Biases

## File preparation per study
Since we are using 2'O-methylation maps from different publications, file preparation differs slightly based on the original data provided.
### Nanopore-DRS (Li et. al. 2024)
Table S5 (excel file) contains a list of Nm sites detected in HEK293T, HeLa, and HepG2 cell lines. Each sheet contains the results for each cell line. The dataset was missing strand and mRNA consequence (5UTR, CDS, or 3UTR) information. 

We used this awk code to remove any Nm sites with a ratio below 10%. 
```
awk -F'\t' 'NR==1 || $3 >= 0.10 { print $1, $2 }' OFS='\t' HEKallData.txt > HEK_rawFiltered.txt
awk -F'\t' 'NR==1 || $3 >= 0.10 { print $1, $2 }' OFS='\t' HeLaallData.txt > HeLa_rawFiltered.txt
awk -F'\t' 'NR==1 || $3 >= 0.10 { print $1, $2 }' OFS='\t' HepG2allData.txt > HepG2_rawFiltered.txt
```
Awk was used to reformat and sort BED files to prepare for downstream analysis.
```
awk -F '\t' '{
  split($1, a, "_"); # original format is "chr1_12345"
  chrom = substr(a[1], 4);  # remove "chr"
  start = a[2];
  end = start + 1;
  gene = $2; # extract gene name
  print chrom "\t" start "\t" end "\t" gene;
}' HEK_rawFiltered.txt | sort -k1,1n -k2,2n > HEKconsequence10.bed

awk -F '\t' '{
  split($1, a, "_"); # original format is "chr1_12345"
  chrom = substr(a[1], 4);  # remove "chr"
  start = a[2];
  end = start + 1;
  gene = $2; # extract gene name
  print chrom "\t" start "\t" end "\t" gene;
}' HepG2_rawFiltered.txt | sort -k1,1n -k2,2n > HepG2consequence10.bed

awk -F '\t' '{
  split($1, a, "_"); # original format is "chr1_12345"
  chrom = substr(a[1], 4);  # remove "chr"
  start = a[2];
  end = start + 1;
  gene = $2; # extract gene name
  print chrom "\t" start "\t" end "\t" gene;
}' HeLa_rawFiltered.txt | sort -k1,1n -k2,2n > HeLaconsequence10.bed
```
Since the raw data did not include strand, we need to merge the data with a GTF file that contains the strand.

A python script was used to convert Ensembl's Homo_sapiens.GRCh38.113.chr.gtf file into a BED file with genomic coordinates per mRNA consequence for each gene, genomeConsequence0based.bed
```
python gtfToBed.py
```
Then, bedtools intersect was used to overlay Nm site BED files and genomeConsequence0based.bed.
```
bedtools intersect -a HEKconsequence10.bed -b genomeConsequence0based.bed -wa -wb > HEKconsequencetmp0based.bed
bedtools intersect -a HepG2consequence10.bed -b genomeConsequence0based.bed -wa -wb > HepG2consequencetmp0based.bed
bedtools intersect -a HeLaconsequence10.bed -b genomeConsequence0based.bed -wa -wb > HeLaconsequencetmp0based.bed
```
After bedtools, we added the strand to Nm site files by retaining Nm sites where the gene reported at the Nm site matched the GTF file gene. 
```
awk '$4 == $9' HEKconsequencetmp0based.bed > HEKlocation0based.bed
awk '$4 == $9' HepG2consequencetmp0based.bed > HepG2location0based.bed
awk '$4 == $9' HeLaconsequencetmp0based.bed > HeLalocation0based.bed
```
Since different transcripts from the same gene may have different mRNA consequences at an Nm site, this script only includes the most common consequence per Nm site. 
```
python addConsequence.py -i HEKlocation0based.bed -o HEKlocation0based_filtered.bed
python addConsequence.py -i HepG2location0based.bed -o HepG2location0based_filtered.bed
python addConsequence.py -i HeLalocation0based.bed -o HeLalocation0based_filtered.bed
```
We have added the strand and mRNA consequence information, but still need to reformat the BED file for downstream analysis.
```
awk 'NR > 1 { print $1, $2, $3 - 1, $4, $8, $10 }' OFS="\t" HEKlocation0based_filtered.bed > HEK0based_allSites.bed
awk 'NR > 1 { print $1, $2, $3 - 1, $4, $8, $10 }' OFS="\t" HeLalocation0based_filtered.bed > HeLa0based_allSites.bed
awk 'NR > 1 { print $1, $2, $3 - 1, $4, $8, $10 }' OFS="\t" HepG2location0based_filtered.bed > HepG20based_allSites.bed
```
These BED files can now be analyzed in parallel with the datasets from the other publications.

### NJU-seq (Tang et. al. 2024)
Table S2 (excel file) contains Nm sites detected in HEK293T, HeLa, and A549 cell lines. Each sheet contains the results for each cell line. Columns 'Chr', 'Position', 'Strand', 'Nm', 'ID', 'Name', and 'Distribution' were extracted from the original excel file into ```A549_Filtered.txt, HEK_Filtered.txt, and HeLa_Filtered.txt```. These files needed to be reformatted for downstream analysis.
```
awk 'BEGIN{OFS="\t"} {if ($6=="Forward") $6="+"; else if ($6=="Reverse") $6="-"; print}' A549_Filtered.txt | sort -k1,1 -k2,2n > A549_allSites.bed
awk 'BEGIN{OFS="\t"} {if ($6=="Forward") $6="+"; else if ($6=="Reverse") $6="-"; print}' HEK_Filtered.txt | sort -k1,1 -k2,2n > HEK_allSites.bed
awk 'BEGIN{OFS="\t"} {if ($6=="Forward") $6="+"; else if ($6=="Reverse") $6="-"; print}' HeLa_Filtered.txt | sort -k1,1 -k2,2n > HeLa_allSites.bed
```
The data can now be analyzed in parallel with the datasets from the other publications.

### Nm-mut-seq (Chen et. al. 2023)
All supplementary tables were combined into one excel file. Each sheet contains the results for each cell line. Sheets "S5_HeLa_mRNA_WT" and "S6_HepG2_mRNA_WT" contain all necessary genomic information. Columns 'chr', 'position', 'position2', 'mRNA segment', 'name', and 'strand' were extracted from the original excel file into ```HepG2_allSites.txt and HeLa_allSites.txt```. Files were then sorted for downstream analysis.
```
sort -k1,1 -k2,2n HepG2_allSites.txt > HepG2_allSites.bed
sort -k1,1 -k2,2n HeLa_allSites.txt > HeLa_allSites.bed
```
Now these files can be analyzed in parallel with the datasets from the other publications.

### Nm-seq (Dai et. al. 2018)
Nm site maps were downloaded from GSE90164. These BED files contained strand information, but not mRNA consequence. Each BED file contains the results for each cell line (HeLa and HEK293T). 

Bedtools intersect was used to add all mRNA consequences per Nm site. (see Nanopore-DRS section for the generation of genomeConsequence0based.bed)
```
bedtools intersect -a GSE90164_HeLamRNA.Nm.genome_sorted.bed -b genomeConsequence0based.bed -wa -wb -s > HeLaNm.hg38noLiftover_annotated.bed
bedtools intersect -a GSE90164_HEKmRNA.Nm.genome_sorted.bed -b genomeConsequence0based.bed -wa -wb -s > HEKNm.hg38noLiftover_annotated.bed
```
Since different transcripts from the same gene may have different mRNA consequences at an Nm site, this script only includes the most common consequence per Nm site. 
```
python extractConsequence.py -i HeLaNm.hg38noLiftover_annotated.bed -o HeLaNm.hg38noLiftover_filtered.bed
python extractConsequence.py -i HEKNm.hg38noLiftover_annotated.bed -o HEKNm.hg38noLiftover_filtered.bed
```
Now that we have added mRNA consequence, the data from this study can be analyzed in parallel with the datasets from the other publications.

## Extracting the sequence flanking Nm sites
Nm maps from all studies and cell lines were renamed to distinguish studies and cell lines as ```mlmHEKSites.bed, mlmHeLaSites.bed, mlmHepG2Sites.bed, njuA549Sites.bed, njuHeLaSites.bed, njuHEKSites.bed, nmHEKSites.bed, nmHeLaSites.bed, nmmtHeLaSites.bed, and nmmtHepG2Sites.bed```. Then, the sequence flanking Nm sites in the 5'UTR were extracted using a python script.
```
python extractSeq.py -i "$BED" -o "$BED_SEQ" --genome "Homo_sapiens.GRCh38.113"
```
Once the sequence was extracted, another script was used to identify Nm sites around near-cognate codons and categorize them by Nm position relative to the near-cognate codon.
```
python cognates.py -i "$BED_SEQ" -o "$BED_COGNATES"
```
## Extracting Nm[+1] Sites from other Nm positions
We extracted Nm[+1] sites into ```mlmHEKCodons.bed, mlmHeLaCodons.bed, mlmHepG2Codons.bed, njuA549Codons.bed, njuHeLaCodons.bed, njuHEKCodons.bed, nmHEKCodons.bed, nmHeLaCodons.bed, nmmtHeLaCodons.bed, and nmmtHepG2Codons.bed```. From these files, we extracted a list of all genes containing Nm[+1] and counts the occurance of each codon with the script below. 
```
python analyzedGenes.py
```
## Identifying near-cognate codon frequency in 5'UTRs
Using a list of genes that were known to be expressed in HR Ribo-seq, we extracted a FASTA file containing 5'UTR sequences for genes of interest from Ensembl's biomart. 
```
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" formatter = "FASTA" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
				
	<Dataset name = "hsapiens_gene_ensembl" interface = "default" >
		<Filter name = "external_gene_name" value = riboSeqGenes.txt />
		<Attribute name = "5utr" />
		<Attribute name = "chromosome_name" />
		<Attribute name = "start_position" />
		<Attribute name = "end_position" />
		<Attribute name = "external_gene_name" />
		<Attribute name = "ensembl_transcript_id" />
		<Attribute name = "transcript_length" />
		<Attribute name = "strand" />
	</Dataset>
</Query>
```
Using this FASTA file, we selected for the longest recorded 5'UTR per gene. The following script outputs a count of each codon present.
```
python parseFastaRibo.py
```
The following script outputs a count of every occurance of a codon.
```
python countCodonsRibo.py
```
## Analyzing codon biases of near-cognate codons for Nm[+1] vs 5'UTRs of expressed genes
```codonBias.Rmd``` generates a barplot comparing codon biases in near-cognate codons. The markdown is codonBias.html.
The data sets needed to reproduce the plots are:
* 5utrCodonCountsRiboSeq.csv
* pos1CodonFreq.txt