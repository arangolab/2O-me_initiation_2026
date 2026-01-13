# Supplementary Data 6G

## File preparation per study
Since we are using 2'O-methylation maps from different publications, file preparation differs slightly based on the original data provided.
### Nanopore-DRS (Li et. al. 2024)
Table S5 (excel file) contains a list of Nm sites in HeLa cells. The data was missing each Nm site's strand and mRNA consequence (5UTR, CDS, or 3UTR) information. 

#### Preparing HeLa Dataset
We used this awk code to remove any Nm sites with a ratio below 10%. 
```
awk -F'\t' 'NR==1 || $3 >= 0.10 { print $1, $2 }' OFS='\t' HeLaallData.txt > HeLa_rawFiltered.txt
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
}' HeLa_rawFiltered.txt | sort -k1,1n -k2,2n > HeLaconsequence10.bed
```
#### Appending strand & mRNA consequence to HeLa
Since the raw data did not include strand, we need to merge the data with a GTF file that contains the strand.

A python script was used to convert Ensembl's Homo_sapiens.GRCh38.113.chr.gtf file into a BED file with genomic coordinates per mRNA consequence for each gene, genomeConsequence0based.bed
```
python gtfToBed.py
```
*script in S6A directory*

Then, bedtools intersect was used to overlay Nm site BED files and genomeConsequence0based.bed.
```
bedtools intersect -a HeLaconsequence10.bed -b genomeConsequence0based.bed -wa -wb > HeLaconsequencetmp0based.bed
```
After bedtools, we added the strand to Nm site files by retaining Nm sites where the gene reported at the Nm site matched the GTF file gene. 
```
awk '$4 == $9' HeLaconsequencetmp0based.bed > HeLalocation0based.bed
```
Since different transcripts from the same gene may have different mRNA consequences at an Nm site, this script only includes the most common consequence per Nm site. 
```
python addConsequence.py -i HeLalocation0based.bed -o HeLalocation0based_filtered.bed
```
*script in S6A directory*

Now we have added the strand and mRNA consequence!
## Extracting the sequence flanking Nm sites
Nm maps were renamed to distinguish cell lines as ```mlmHeLaSites.bed```. Then, the sequence flanking Nm sites in the 5'UTR were extracted using a python script.
```
python extractSeq.py -i "$BED" -o "$BED_SEQ" --genome "Homo_sapiens.GRCh38.113"
```
*script in S7B directory*

Once the sequence was extracted, another script was used to identify Nm sites around near-cognate codons and categorize them by Nm position relative to the near-cognate codon.
```
python cognates.py -i "$BED_SEQ" -o "$BED_COGNATES"
```
*script in S7B directory*

## Extracting Nm[+1] Sites from other Nm positions
We extracted Nm[+1] sites into ```mlmHeLaCodons.bed```. From these files, we extracted a list of all genes containing Nm[+1] and counts the occurance of each codon with the script below. 
```
python analyzedGenesHeLaNano.py
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
```codonBiasHeLaNanoBiasesOnly.Rmd``` generates a barplot comparing codon biases in near-cognate codons. The markdown is codonBiasHeLaNanoBiasesOnly.html.
The data sets needed to reproduce the plots are:
* 5utrCodonCounts.csv
* pos1CodonFreqHELANANO.txt
