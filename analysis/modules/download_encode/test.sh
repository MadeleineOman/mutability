#first agrument = tissue 
#second argument = predictor 
#thir argument = filename 
#fourht argument = web link to file 

touch data/$1/track_data/$2/wget_output.txt

wget $4 -P data/$1/track_data/$2/ -o data/$1/track_data/$2/wget_output.txt

#convert from bigwig to bed 
data/applications/UCSC_genomeBrowser_Blat/bigWigToBedGraph data/$1/track_data/$2/$3 data/$1/track_data/$2/$2_hg38.bed
rm data/$1/track_data/$2/$3
# wc -l data/$1/track_data/$2/$2.bed

#prepr for liftovber: get the liftover chain files (n eed to do it in 2 steps, no hg38-->hg18 chain file exists....)
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz -P data/$1/track_data/$2/ -o data/$1/track_data/$2/hg38ToHg19_wget_output.txt
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg18.over.chain.gz -P data/$1/track_data/$2/ -o data/$1/track_data/$2/hg19ToHg18_wget_output.txt

# #liftvoer hg38-->19 
data/applications/UCSC_genomeBrowser_Blat/liftOver data/$1/track_data/$2/$2_hg38.bed data/$1/track_data/$2/hg38ToHg19.over.chain.gz data/$1/track_data/$2/$2_hg19.bed  data/$1/track_data/$2/$2_unliftgedhg38.bed
# # wc -l data/$1/track_data/$2/$2* 
rm data/$1/track_data/$2/$2_hg38.bed #removing the unused file
# # echo "19M rows lost in hg38-->hg19 liftover out of 218M"

# #liftover hg19-->18 
data/applications/UCSC_genomeBrowser_Blat/liftOver data/$1/track_data/$2/$2_hg19.bed data/$1/track_data/$2/hg19ToHg18.over.chain.gz data/$1/track_data/$2/$2_unsorted.bed  data/$1/track_data/$2/$2_unliftgedhg19.bed
rm data/$1/track_data/$2/$2_hg19.bed #removing the unused file 
# # wc -l data/$1/track_data/$2/$2* # 
# # echo "450k rows lost in hg19-->hg18 liftover out of 218M"

# #sort
sort -k1,1 -k2,2n data/$1/track_data/$2/$2_unsorted.bed > data/$1/track_data/$2/$2.bed 
rm data/$1/track_data/$2/$2_unsorted.bed

# bgzip and tabix 
bgzip data/$1/track_data/$2/$2.bed
tabix -p bed data/$1/track_data/$2/$2.bed.gz
