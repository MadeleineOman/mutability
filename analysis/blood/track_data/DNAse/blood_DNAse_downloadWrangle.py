# **ReadME for the DNAse data** 
# - downloaded from: 
#     - UCSC: http://hgdownload.cse.ucsc.edu/goldenpath/hg18/encodeDCC/wgEncodeUwDnaseSeq/
# - download date: 
#     - lol it depends when i finally do the thing 
# - genome build 
#     - hg18 
# - tissue 
#     - gm12878: 
#     - a lymphoblastoid cell line: lympthoblastoid is a typy of/ modified lymphocyte, which is a white blood cell
#     - read more: https://www.genome.gov/encode-project-common-cell-types
# - lab method: 
#     - Digital DNaseI methodology. read more : http://genome.ucsc.edu/cgi-bin/hgFileUi?db=hg18&g=wgEncodeUwDnaseSeq
# - file format & head: 
#     - .wig file, variable step, whole genome 
#     - file head 
#         - variableStep chrom=chr1 span=20
#         - 105    1
#         - 125    1
#         - 145    1
#         - 165    1
#     - there are 2 available replicates, we are only using rep1 here 
# - what the data represents: 
#     - DNaseI sensitivity is directly reflected in raw tag density (Signal), which is shown in the track as density of tags mapping within a 150 bp sliding window (at a 20 bp step across the genome)


#DOWNLOAD THE FILES 
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg18/encodeDCC/wgEncodeUwDnaseSeq/wgEncodeUwDnaseSeqRawSignalRep1Gm12878.wig.gz -P ../../../../data/blood/track_data/DNAse
    
#CONVERt TO .bed USING BEDOPS wig2bed
#The wig2bed script converts both variable - and fixed -step, 1-based, closed [start, end] UCSC Wiggle format (WIG) to sorted, 0-based, half-open [start-1, end) extended BED data.
#In the case where WIG data are sourced from bigWigToWig or other tools that generate 0-based, half-open [start-1, end) WIG, a --zero-indexed option is provided to generate coordinate output without any re-indexing.
gunzip ../../../../data/blood/track_data/DNAse/wgEncodeUwDnaseSeqRawSignalRep1Gm12878.wig.gz 
wig2bed --multisplit < ../../../../data/blood/track_data/DNAse/wgEncodeUwDnaseSeqRawSignalRep1Gm12878.wig > ../../../../data/blood/track_data/DNAse/wgEncodeUwDnaseSeqRawSignalRep1Gm12878.bed #don;t ask me why but the gunzip piped into wig2bed doesn;t work. need the steps seperately 

#TABIX THE WHOLE BED FILE
bgzip ../../../../data/blood/track_data/DNAse/wgEncodeUwDnaseSeqRawSignalRep1Gm12878.bed
tabix -s 1 -b 2 -e 3 ../../../../data/blood/track_data/DNAse/wgEncodeUwDnaseSeqRawSignalRep1Gm12878.bed.gz

#SEPERATE INTO SEPERATE CHROMOSOME FILES 
for f in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22; do 
    tabix ../../../../data/blood/track_data/DNAse/wgEncodeUwDnaseSeqRawSignalRep1Gm12878.bed.gz $f > ../../../../data/blood/track_data/DNAse/DNAse_$f.bed
    bgzip ../../../../data/blood/track_data/DNAse/DNAse_$f.bed
    tabix -s 1 -b 2 -e 3 ../../../../data/blood/track_data/DNAse/DNAse_$f.bed.gz
done