# Supplementary Data 6D
## File preparation
### Nanopore-DRS (Li et. al. 2024)
Table S4 (excel file) contains a list of Nm sites in C4-2 cells and Table S5 (excel file) contains a list of Nm sites in HeLa cells. The data was missing each Nm site's strand and mRNA consequence (5UTR, CDS, or 3UTR) information. 

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
#### Preparing the C4-2 Dataset
Data is already filtered for Nm sites with ratio > 10%. The python script below reformats the data for downstream analysis with bedtools.
```
python formatForBedtools_NOFBL.py
```
*script in S6A directory*
#### Appending strand & mRNA consequence to HeLa & C4-2 datasets
Since the raw data did not include strand, we need to merge the data with a GTF file that contains the strand.

A python script was used to convert Ensembl's Homo_sapiens.GRCh38.113.chr.gtf file into a BED file with genomic coordinates per mRNA consequence for each gene, genomeConsequence0based.bed
```
python gtfToBed.py
```
*script in S6A directory*

Then, bedtools intersect was used to overlay Nm site BED files and genomeConsequence0based.bed.
```
bedtools intersect -a HeLaconsequence10.bed -b genomeConsequence0based.bed -wa -wb > HeLaconsequencetmp0based.bed
bedtools intersect -a formattedAllC42sites.bed -b genomeConsequence0based.bed -wa -wb > formattedAllC42sitestmp0based.bed
```
After bedtools, we added the strand to Nm site files by retaining Nm sites where the gene reported at the Nm site matched the GTF file gene. 
```
awk '$4 == $9' HeLaconsequencetmp0based.bed > HeLalocation0based.bed
awk '$4 == $9' formattedAllC42sitestmp0based.bed > formattedAllC42sites0based.bed
```
Since different transcripts from the same gene may have different mRNA consequences at an Nm site, this script only includes the most common consequence per Nm site. 
```
python addConsequence.py -i HeLalocation0based.bed -o HeLalocation0based_filtered.bed
python addConsequence.py -i formattedAllC42sites0based.bed -o AllC42sites0based_filtered.bed
```
*script in S6A directory*

Now we have added the strand and mRNA consequence!
## Extracting the sequence flanking Nm sites
Nm maps were renamed to distinguish cell lines as ```mlmHeLaSites.bed, AllC42sites0based_allSites.bed```. Then, the sequence flanking Nm sites in the 5'UTR were extracted using a python script.
```
python extractSeq.py -i "$BED" -o "$BED_SEQ" --genome "Homo_sapiens.GRCh38.113"
```
*script in S6B directory*

Once the sequence was extracted, another script was used to identify Nm sites around near-cognate codons and categorize them by Nm position relative to the near-cognate codon.
```
python cognates.py -i "$BED_SEQ" -o "$BED_COGNATES"
```
*script in S6B directory*

## Upset plot for 5'UTR Nm[+1] sites
```NmPlus1inHeLaAllC42UpsetPlot.Rmd``` generates an upset plot showing common 5'UTR Nm[+1] sites between HeLa and C4-2 cell lines. The markdown is NmPlus1inHeLaAllC42UpsetPlot.html.
The data sets needed to reproduce the plots are:
* AllC42Cognates.bed (file in S6B)
* mlmHeLaCognates.bed (file in S6B)
