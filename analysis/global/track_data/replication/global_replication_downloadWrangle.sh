### **download and wrangle the koren replication timing file**

# - downloaded from: 
#     - https://www.encodeproject.org/experiments/ENCSR098AZD/
# - download date: 
# - genome build 
#     - hg 19, converted to hg18 
# - tissue 
#     - HeLa-S3 G2 phase
# - lab method: 
#     - https://www.encodeproject.org/documents/50ccff70-1305-4312-8b09-0311f7681881/@@download/attachment/wgEncodeUwRepliSeq.html.pdf
# - file format & head: 
# - what the data represents: 
#     - Percentage-normalized Signal: Replication signal at 1 kb intervals as a percentage of normalized +/ -25 kb tag densities for all cell cycle fractions (G1/G1b, S1, S2, S3, S4, G2).

#wget 
wget https://www.encodeproject.org/files/ENCFF001GOJ/@@download/ENCFF001GOJ.bigWig -P data/global/track_data/replication/

#convert to bed file. below is the description for the UCSC utilities readme file in data/applications/UCSC_genomeBrowser_Blat/aboutExecutibles.txt
# ================================================================
# ========   bigWigToBedGraph   ====================================
# ================================================================
# ### kent source version 416 ###
# bigWigToBedGraph - Convert from bigWig to bedGraph format.
# usage:
#    bigWigToBedGraph in.bigWig out.bedGraph
# options:
#    -chrom=chr1 - if set restrict output to given chromosome
#    -start=N - if set, restrict output to only that over start
#    -end=N - if set, restict output to only that under end
#    -udcDir=/dir/to/cache - place to put cache for remote bigBed/bigWigs
data/applications/UCSC_genomeBrowser_Blat/bigWigToBedGraph data/global/track_data/replication/ENCFF001GOJ.bigWig data/global/track_data/replication/ENCFF001GOJ.bed

### **using liftover**
# - usage: liftOver oldFile map.chain newFile
#     - http://hgdownload.soe.ucsc.edu/downloads.html#liftover
#         - The links to liftOver over.chain files can be found in the corresponding assembly sections above. For example, the link for the mm5-to-mm6 over.chain file is located in the mm5 downloads section
#     - found the chain file in the"liftover": files in UCSC doiwnloads page under hg19 http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/

#get the chain over file for liftover 
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg18.over.chain.gz -P data/global/track_data/replication/

#use liftover 
data/applications/UCSC_genomeBrowser_Blat/liftOver data/global/track_data/replication/ENCFF001GOJ.bed data/global/track_data/replication/hg19ToHg18.over.chain.gz data/global/track_data/replication/replication_unsorted.bed  data/global/track_data/replication/unliftged.bed 

#i double checked it manually against UCSC web application using entire chr3 
#0) create test files (smaller) to use on the ucsc web application 
#grep "chr3" ../../../../data/global/track_data/replication/ENCFF001GOJ.bed > ../../../../data/global/track_data/replication/liftoverTest_hg19_chr3.bed
#grep "chr3" ../../../../data/global/track_data/replication/replication.bed > ../../../../data/global/track_data/replication/liftoverTest_hg18_chr3.bed
#1) download the liftoverTest_hg19_chr3.bed file
#2) convert to hg18 on https://genome.ucsc.edu/cgi-bin/hgLiftOver
#3) reupload and then compare differences with diff command: 
#diff ../../../../data/global/track_data/replication/liftoverTest_hg18_chr3.bed ../../../../data/global/track_data/replication/hglft_genome_54b4c_c8e0f0.bed

#sort the file 
sort -k1,1 -k2,2n data/global/track_data/replication/replication_unsorted.bed > data/global/track_data/replication/replication.bed 

#bgzip and tabix 
bgzip data/global/track_data/replication/replication.bed
tabix -p bed data/global/track_data/replication/replication.bed.gz