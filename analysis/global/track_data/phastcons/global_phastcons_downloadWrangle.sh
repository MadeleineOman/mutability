### phastcons download and wrangle 

# - downloaded from: 
#     - http://hgdownload.soe.ucsc.edu/goldenPath/hg18/phastCons44way/primates/
# - download date: 
#     - lol it depends when i finally do the thing 
# - genome build 
#     - hg18 
# - tissue 
#     - NA
# - lab method: 
#     - This track shows multiple alignments of 44 vertebrate species and measurements of evolutionary conservation using two methods (phastCons and phyloP) from the PHAST package, for all species (vertebrate) and two subsets (primate and placental mammal). The multiple alignments were generated using multiz and other tools in the UCSC/Penn State Bioinformatics comparative genomics alignment pipeline. Conserved elements identified by phastCons are also displayed in this track. 
#     - read more: https://genome.ucsc.edu/cgi-bin/hgTables (UCSC table browser: hg18 > comparitive genomics > conservation > primatecons > describe table schema ) 
# - file format & head: 
#     - wig format, turned into bed
#     - seperate chroms as the file are HUGE 
# - what the data represents: 
#     - PhastCons (which has been used in previous Conservation tracks) is a hidden Markov model-based method that estimates the probability that each nucleotide belongs to a conserved element, based on the multiple alignment. It considers not just each individual alignment column, but also its flanking columns. By contrast, phyloP separately measures conservation at individual columns, ignoring the effects of their neighbors. As a consequence, the phyloP plots have a less smooth appearance than the phastCons plots, with more "texture" at individual sites. The two methods have different strengths and weaknesses. PhastCons is sensitive to "runs" of conserved sites, and is therefore effective for picking out conserved elements. PhyloP, on the other hand, is more appropriate for evaluating signatures of selection at particular nucleotides or classes of nucleotides (e.g., third codon positions, or first positions of miRNA target sites). 


#wget 
for f in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22; do 
    wget http://hgdownload.soe.ucsc.edu/goldenPath/hg18/phastCons44way/primates/$f.phastCons44way.primates.wigFix.gz -P data/global/track_data/phastcons/
done    

#gunzip 
for f in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22; do 
    gunzip data/global/track_data/phastcons/$f.phastCons44way.primates.wigFix.gz 
done 

#conver to bed, then bgzip right away (the files get HUGE)
for f in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22; do 
    convert2bed -i wig < data/global/track_data/phastcons/$f.phastCons44way.primates.wigFix >  data/global/track_data/phastcons/phastcons_$f.bed
    bgzip data/global/track_data/phastcons/phastcons_$f.bed
done   

#tabix 
for f in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22; do 
    tabix -p bed data/global/track_data/phastcons/phastcons_$f.bed.gz
done 

#rm uneccessary files to get back more space 
rm  data/global/track_data/phastcons/*.wigFix