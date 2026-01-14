# Supplementary Data 7C

## File preparation
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
Now we have added the strand and mRNA consequence!
## Extracting the sequence flanking Nm sites
Nm maps were renamed as ```mlmHeLaSites.bed```. Then, the sequence flanking Nm sites in the 5'UTR were extracted using a python script.
```
python extractSeq.py -i "$BED" -o "$BED_SEQ" --genome "Homo_sapiens.GRCh38.113"
```
*script in S6B directory*

Once the sequence was extracted, another script was used to identify Nm sites around near-cognate codons and categorize them by Nm position relative to the near-cognate codon.
```
python cognates.py -i "$BED_SEQ" -o "$BED_COGNATES"
```
*script in S6B directory*

We now have a list of 5'UTR Nm[+1] sites, including gene name for downstream analysis of translation efficiency (TE).
## Deriving Translation Efficiency
### Data extraction
We extracted publically available data from GSE105248, a study that performed RNA-seq with matched Ribo-seq in HeLa cells. These cells are modified for doxycylin induced FBL knockdown. The raw RNA-seq counts were provided by the author, but we need to derive RPFs to calculate TE.
### Deriving RPFs
Raw Ribo-seq FASTQ sequence files were trimmed to remove adapters using TrimGalore/0.6.10
```
trim_galore --fastqc -o "$TRIMMED_FASTQ" "$FASTQ"
```
Trimmed FASTQ files were aligned to GRCh38 using STAR/2.7.9a
```
STAR --runMode alignReads \
--genomeDir /genome_index/GRCh38 \
--readFilesIn $TRIMMED_FASTQ \
--outFilterMultimapNmax 3 \
--outFilterMultimapScoreRange 1 \
--outSAMprimaryFlag AllBestScore \
--outFileNamePrefix $OUTPUT_PREFIX \
--outSAMattributes All \
--outSAMtype BAM SortedByCoordinate \
--outFilterMismatchNmax 2 \
--outFilterType BySJout \
--outReadsUnmapped Fastx \
--outFilterScoreMin 10 \
--outSAMattrRGline ID:foo SM:$SAMPLE_NAME \
--alignEndsType EndToEnd \
--quantMode TranscriptomeSAM \
--readFilesCommand zcat
```
Generate RPF counts with featureCounts in subread/2.0.3
```
featureCounts \
    -a "$GTF" \
    -o counts.txt \
    -g gene_id \
    -s 0 \
    "${FILES[@]}"
```
Now we have RPF counts we can use to calculate translation efficiency!
## Plotting TE in FBL-KD vs WT depending on Nm status
```NmPlus1HeLaTELog2Fold.Rmd``` generates a boxplot comparing TE in 5'UTR Nm[+1] genes vs genes without 5'UTR Nm. The script normalizes the RNA-seq and Ribo-seq counts, calculates TE, and analyzes fold change between FBL-KD and WT cells.
The data sets needed to reproduce the plots are:
* Download raw RNA-seq counts from GSE105248
* HeLalocation0based_filtered.bed (file in S6A)
* mlmHeLaCognates.bed (file in S6B)
