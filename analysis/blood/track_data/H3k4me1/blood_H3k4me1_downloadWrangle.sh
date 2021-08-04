# Read me for H3k4me1

# **new data source**
# - downloaded from: 
#     - UCSC: http://hgdownload.soe.ucsc.edu/goldenPath/hg18/encodeDCC/wgEncodeBroadChipSeq/
# - download date: 
#     - lol it depends when i finally do the thing 
# - genome build 
#     - hg18 
# - tissue 
#     - gm12878: 
#     - a lymphoblastoid cell line: lympthoblastoid is a typy of/ modified lymphocyte, which is a white blood cell
#     - read more: https://www.genome.gov/encode-project-common-cell-types
# - lab method: 
#     - H3k4me1 Chip-seq. read more on general histone Chip seq:: http://genome.ucsc.edu/cgi-bin/hgFileUi?db=hg18&g=wgEncodeBroadChipSeq
# - file format & head: 
#     - wig format 
#     - fixedStep chrom=chr1 start=1 step=25
#     - 19
#     - 21
#     - 33.5
#     - 39.5
#     - 36.5
#     - 33.5
#     - 34.5
#     - 31.5
#     - 32.25
# - what the data represents: 
#     - General histone Chip seq: Density graph (wiggle) of signal enrichment based on processed data.
#     - General histone Chip seq: Fragment densities were computed by counting the number of reads overlapping each position in the genome (counted as 1 from the start of each alignment to 200 bp downstream and by 0.25 from 200 to 300 bp, and displayed at 25bp resolution). Discrete intervals of ChIP-seq fragment enrichment were identified using a scan statistics approach, under the assumption of uniform background signal


#BASH SCRIPT 

#download the files 
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg18/encodeDCC/wgEncodeBroadChipSeq/wgEncodeBroadChipSeqSignalGm12878H3k4me1.wig.gz -P data/blood/track_data/H3k4me1/

#convert to bed 
#The wig2bed script converts both variable - and fixed -step, 1-based, closed [start, end] UCSC Wiggle format (WIG) to sorted, 0-based, half-open [start-1, end) extended BED data.
#In the case where WIG data are sourced from bigWigToWig or other tools that generate 0-based, half-open [start-1, end) WIG, a --zero-indexed option is provided to generate coordinate output without any re-indexing. 
gunzip data/blood/track_data/H3k4me1/wgEncodeBroadChipSeqSignalGm12878H3k4me1.wig.gz
wig2bed --multisplit < data/blood/track_data/H3k4me1/wgEncodeBroadChipSeqSignalGm12878H3k4me1.wig > data/blood/track_data/H3k4me1/wgEncodeBroadChipSeqSignalGm12878H3k4me1.bed

#tabix big file 
bgzip data/blood/track_data/H3k4me1/wgEncodeBroadChipSeqSignalGm12878H3k4me1.bed
tabix -s 1 -b 2 -e 3 data/blood/track_data/H3k4me1/wgEncodeBroadChipSeqSignalGm12878H3k4me1.bed.gz

#for loop to sep into chroms 
for f in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22
do 
tabix data/blood/track_data/H3k4me1/wgEncodeBroadChipSeqSignalGm12878H3k4me1.bed.gz $f > data/blood/track_data/H3k4me1/H3k4me1_$f.bed
bgzip data/blood/track_data/H3k4me1/H3k4me1_$f.bed
tabix -s 1 -b 2 -e 3 data/blood/track_data/H3k4me1/H3k4me1_$f.bed.gz
done