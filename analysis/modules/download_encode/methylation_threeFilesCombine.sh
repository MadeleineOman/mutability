### blood methylation download wrangle 

# - genome build 
#     - GrCh38
# - lab method:  
#     - WG bisulfite seqeuncing --> turns cytosines into uracil but not methylated cytosines. then WGS to see where there are still Cs 
#     - read more : https://www.encodeproject.org/documents/964e2676-d0be-4b5d-aeec-f4f02310b221/@@download/attachment/WGBS%20pipeline%20overview.pdf
# - file format & head: 
#     - OG: bigbed format (type of condesned binary) 
# - what the data represents: 
#     - Description of bedMethyl file: The bedMethyl file is a bed9+2 file containing the number of reads and the percent methylation. Each column represents the following: 
#         1. Reference chromosome or scaffold 
#         2. Start position in chromosome 
#         3. End position in chromosome 
#         4. Name of item 
#         5. Score from 0-1000. Capped number of reads 
#         6. Strandedness, plus (+), minus (-), or unknown (.) 
#         7. Start of where display should be thick (start codon) 
#         8. End of where display should be thick (stop codon) 
#         9. Color value (RGB) 
#         10. Coverage, or number of reads 
#         11. Percentage of reads that show methylation at this position in the genome 


#wget 
wget https://www.encodeproject.org/files/$2/@@download/$2.bigBed -P data/$1/track_data/methylation/ -o data/$1/track_data/methylation/wget_output_CHG.txt
wget https://www.encodeproject.org/files/$3/@@download/$3.bigBed -P data/$1/track_data/methylation/ -o data/$1/track_data/methylation/wget_output_CHH.txt
wget https://www.encodeproject.org/files/$4/@@download/$4.bigBed -P data/$1/track_data/methylation/ -o data/$1/track_data/methylation/wget_output_CpG.txt

#convert to bed 
echo $1 "converting to bed" 
data/applications/UCSC_genomeBrowser_Blat/bigBedToBed data/$1/track_data/methylation/$2.bigBed data/$1/track_data/methylation/methylation_CHG_hg38.bed
data/applications/UCSC_genomeBrowser_Blat/bigBedToBed data/$1/track_data/methylation/$3.bigBed data/$1/track_data/methylation/methylation_CHH_hg38.bed
data/applications/UCSC_genomeBrowser_Blat/bigBedToBed data/$1/track_data/methylation/$4.bigBed data/$1/track_data/methylation/methylation_CpG_hg38.bed
rm data/$1/track_data/methylation/$2.bigBed
rm data/$1/track_data/methylation/$3.bigBed
rm data/$1/track_data/methylation/$4.bigBed

#prepare for liftover 
awk '{print $1"\t"$2"\t"$3"\t"$4"|"$5"|"$6"|"$7"|"$8"|"$9"|"$10"|"$11"|"$12"|"$13"|"$14}' data/$1/track_data/methylation/methylation_CHG_hg38.bed > data/$1/track_data/methylation/methylation_CHG_hg38_preppedLiftover.bed
awk '{print $1"\t"$2"\t"$3"\t"$4"|"$5"|"$6"|"$7"|"$8"|"$9"|"$10"|"$11"|"$12"|"$13"|"$14}' data/$1/track_data/methylation/methylation_CHH_hg38.bed > data/$1/track_data/methylation/methylation_CHH_hg38_preppedLiftover.bed
awk '{print $1"\t"$2"\t"$3"\t"$4"|"$5"|"$6"|"$7"|"$8"|"$9"|"$10"|"$11"|"$12"|"$13"|"$14}' data/$1/track_data/methylation/methylation_CpG_hg38.bed > data/$1/track_data/methylation/methylation_CpG_hg38_preppedLiftover.bed
rm data/$1/track_data/methylation/methylation_CHG_hg38.bed
rm data/$1/track_data/methylation/methylation_CHH_hg38.bed
rm data/$1/track_data/methylation/methylation_CpG_hg38.bed

# #liftvoer hg38-->19 
echo $1 "lifting over hg38-->hg19"
data/applications/UCSC_genomeBrowser_Blat/liftOver data/$1/track_data/methylation/methylation_CHG_hg38_preppedLiftover.bed analysis/modules/download_encode/hg38ToHg19.over.chain.gz data/$1/track_data/methylation/methylation_CHG_hg19.bed  data/$1/track_data/methylation/methylation_CHG_unliftgedhg38.bed
data/applications/UCSC_genomeBrowser_Blat/liftOver data/$1/track_data/methylation/methylation_CHH_hg38_preppedLiftover.bed analysis/modules/download_encode/hg38ToHg19.over.chain.gz data/$1/track_data/methylation/methylation_CHH_hg19.bed  data/$1/track_data/methylation/methylation_CHH_unliftgedhg38.bed
data/applications/UCSC_genomeBrowser_Blat/liftOver data/$1/track_data/methylation/methylation_CpG_hg38_preppedLiftover.bed analysis/modules/download_encode/hg38ToHg19.over.chain.gz data/$1/track_data/methylation/methylation_CpG_hg19.bed  data/$1/track_data/methylation/methylation_CpG_unliftgedhg38.bed
rm data/$1/track_data/methylation/methylation_CHG_hg38_preppedLiftover.bed 
rm data/$1/track_data/methylation/methylation_CHH_hg38_preppedLiftover.bed
rm data/$1/track_data/methylation/methylation_CpG_hg38_preppedLiftover.bed

# #liftvoer hg19-->18 
echo $1 "lifting over hg19-->hg18"
data/applications/UCSC_genomeBrowser_Blat/liftOver data/$1/track_data/methylation/methylation_CHG_hg19.bed analysis/modules/download_encode/hg19ToHg18.over.chain.gz data/$1/track_data/methylation/methylation_CHG_hg18.bed  data/$1/track_data/methylation/methylation_CHG_unliftgedhg19.bed
data/applications/UCSC_genomeBrowser_Blat/liftOver data/$1/track_data/methylation/methylation_CHH_hg19.bed analysis/modules/download_encode/hg19ToHg18.over.chain.gz data/$1/track_data/methylation/methylation_CHH_hg18.bed  data/$1/track_data/methylation/methylation_CHH_unliftgedhg19.bed
data/applications/UCSC_genomeBrowser_Blat/liftOver data/$1/track_data/methylation/methylation_CpG_hg19.bed analysis/modules/download_encode/hg19ToHg18.over.chain.gz data/$1/track_data/methylation/methylation_CpG_hg18.bed  data/$1/track_data/methylation/methylation_CpG_unliftgedhg19.bed
rm data/$1/track_data/methylation/methylation_CHG_hg19.bed 
rm data/$1/track_data/methylation/methylation_CHH_hg19.bed
rm data/$1/track_data/methylation/methylation_CpG_hg19.bed

#bcombine and sort 
echo $1 "combining "
cat data/$1/track_data/methylation/methylation_CHG_hg18.bed > data/$1/track_data/methylation/methylation.bed 
cat data/$1/track_data/methylation/methylation_CHH_hg18.bed >> data/$1/track_data/methylation/methylation.bed 
cat data/$1/track_data/methylation/methylation_CpG_hg18.bed >> data/$1/track_data/methylation/methylation.bed 
rm data/$1/track_data/methylation/methylation_CHG_hg18.bed 
rm data/$1/track_data/methylation/methylation_CHH_hg18.bed 
rm data/$1/track_data/methylation/methylation_CpG_hg18.bed 

echo $1 "sorting"
sort -k1,1V -k2,2n -o data/$1/track_data/methylation/methylation.bed data/$1/track_data/methylation/methylation.bed 
#the sort resulted in "chrom bloacks not continuous" . apparently a probablem with sorting the chrom name. can add -V (-k1,1V) to enforce good sort 

#tabix 
echo $1 "tabixing"
bgzip data/$1/track_data/methylation/methylation.bed
tabix -s 1 -b 2 -e 3 data/$1/track_data/methylation/methylation.bed.gz
