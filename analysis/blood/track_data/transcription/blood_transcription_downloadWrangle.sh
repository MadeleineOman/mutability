# - downloaded from: 
#     - https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg18&g=wgEncodeCaltechRnaSeq
#     - http://hgdownload.soe.ucsc.edu/goldenPath/hg18/encodeDCC/wgEncodeCaltechRnaSeq/
# - download date: 
#     - lol it depends when i finally do the thing 
# - genome build 
#     - hg18
# - tissue 
#     - gm12878 
#     - a lymphoblastoid cell line: lympthoblastoid is a typy of/ modified lymphocyte, which is a white blood cell
#     - read more: https://www.genome.gov/encode-project-common-cell-types
# - lab method: 
#     - rna-seq of seelcted poly-a tail RNa (out of total cellular RNA) https://genome.ucsc.edu/cgi-bin/hgTrackUi?g=wgEncodeCaltechRnaSeq
#         -  The 2x75 n.t. reads were mapped serially, first with the Bowtie program (Langmead et al., 2009) against the genome and UCSC known-gene splice junctions (Splice Sites). Bowtie-unmapped reads were then mapped using BLAT to find evidence of novel splicing, by requiring at least 10 bp on the short-side of the splice. obtained as pairs from both ends cDNAs resulting from random priming. https://genome.ucsc.edu/cgi-bin/hgEncodeVocab?type=readType
#         - 200 bp insertion 
# - file format & head: 
#     - track type=bedGraph name="GM_2x75"  priority=0.010 visibility=full color=0,0,255
#     - chr1 35 110 0.0277
#     - chr1 1297 1317 0.0139
#     - chr1 1317 1372 0.0277
#     - chr1 1372 1392 0.0139
#     - chr1 1509 1584 0.0092
# - what the data represents: 
#     - Raw Signal: Density graph (wiggle) of signal enrichment based on a normalized aligned read density (RPKM) for non strand-specific reads. The RPKM measure assists in visualizing the relative amount of a given transcript across multiple samples. 


#download 
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg18/encodeDCC/wgEncodeCaltechRnaSeq/wgEncodeCaltechRnaSeqRawSignalRep1Gm12878CellLongpolyaBb12x75.wig.gz -P data/blood/track_data/transcription/

#remove the header and convert from space delimited to tab delimited 
zcat data/blood/track_data/transcription/wgEncodeCaltechRnaSeqRawSignalRep1Gm12878CellLongpolyaBb12x75.wig.gz | sed 's/ /\t/g' | tail -n +2  > data/blood/track_data/transcription/transcription.bed  

#bgzip 
bgzip data/blood/track_data/transcription/transcription.bed

#tabix 
tabix -s 1 -b 2 -e 3 data/blood/track_data/transcription/transcription.bed.gz