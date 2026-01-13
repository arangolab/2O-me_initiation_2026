# Extended Data 4A

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
Now that we have added the strand and mRNA consequence information, the data from this study can be analyzed in parallel with the datasets from the other publications.

### NJU-seq (Tang et. al. 2024)
Table S2 (excel file) contains Nm sites detected in HEK293T, HeLa, and A549 cell lines. Each sheet contains the results for each cell line. The data from this study already contained all necessary genomic information, therefore it can be analyzed in parallel with the datasets from the other publications.

### Nm-mut-seq (Chen et. al. 2023)
All supplementary tables were combined into one excel file. Each sheet contains the results for each cell line. Sheets "S5_HeLa_mRNA_WT" and "S6_HepG2_mRNA_WT" contain all necessary genomic information, therefore it can be analyzed in parallel with the datasets from the other publications.

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

## Analyzing mRNA consequence per study and cell line
```every2OSITEdotplot.Rmd``` creates dotplots showing the distribution of Nm sites for each mRNA consequence per study and cell line. The markdown is every2OSITEdotplot.html.
The data sets needed to reproduce the plots are:
* NJUSeq_Data.xlsx
* Nm-Mut-seq Supplementary Tables.xlsx
* HeLalocation0based_filtered.bed
* HepG2location0based_filtered.bed
* HEKlocation0based_filtered.bed
* HeLaNm.hg38noLiftover_filtered.bed
* HEKNm.hg38noLiftover_filtered.bed