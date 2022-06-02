#first agrument = tissue 
#second argument = predictor 
#thir argument = filename1 
#fourht argument =  filename 2
#fifth argument = web link to file1
#sixth argument = web link to file 2 


#download and convert file 1 
wget $5 -P data/$1/track_data/$2/ -o data/$1/track_data/$2/$3_wget_output.txt
data/applications/UCSC_genomeBrowser_Blat/bigWigToBedGraph data/$1/track_data/$2/$3.bigWig data/$1/track_data/$2/$2_hg38_$3.bed
rm data/$1/track_data/$2/$3.bigWig
echo "wget and convert file 1 \n"

#download and convert file 2 
wget $6 -P data/$1/track_data/$2/ -o data/$1/track_data/$2/$4_wget_output.txt
data/applications/UCSC_genomeBrowser_Blat/bigWigToBedGraph data/$1/track_data/$2/$4.bigWig data/$1/track_data/$2/$2_hg38_$4.bed
rm data/$1/track_data/$2/$4.bigWig
echo "wget and convert file 2 \n"

#merge the two files 
cat data/$1/track_data/$2/$2_hg38_$3.bed > data/$1/track_data/$2/$2_hg38.bed
wc -l data/$1/track_data/$2/$2_hg38.bed
echo "lines in the first file"
cat data/$1/track_data/$2/$2_hg38_$4.bed >> data/$1/track_data/$2/$2_hg38.bed
wc -l  data/$1/track_data/$2/$2_hg38.bed
echo "lines in the merged file after adding second file"

# wc -l data/$1/track_data/$2/$2.bed

# prepr for liftovber: get the liftover chain files (n eed to do it in 2 steps, no hg38-->hg18 chain file exists....)
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz -P data/$1/track_data/$2/ -o data/$1/track_data/$2/hg38ToHg19_wget_output.txt
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg18.over.chain.gz -P data/$1/track_data/$2/ -o data/$1/track_data/$2/hg19ToHg18_wget_output.txt

# #liftvoer hg38-->19 
data/applications/UCSC_genomeBrowser_Blat/liftOver data/$1/track_data/$2/$2_hg38.bed data/$1/track_data/$2/hg38ToHg19.over.chain.gz data/$1/track_data/$2/$2_hg19.bed  data/$1/track_data/$2/$2_unliftgedhg38.bed
wc -l data/$1/track_data/$2/$2_hg19.bed
rm data/$1/track_data/$2/$2_hg38.bed #removing the unused file
# # echo "19M rows lost in hg38-->hg19 liftover out of 218M"

# #liftover hg19-->18 
data/applications/UCSC_genomeBrowser_Blat/liftOver data/$1/track_data/$2/$2_hg19.bed data/$1/track_data/$2/hg19ToHg18.over.chain.gz data/$1/track_data/$2/$2_unsorted.bed  data/$1/track_data/$2/$2_unliftgedhg19.bed
rm data/$1/track_data/$2/$2_hg19.bed #removing the unused file 
wc -l data/$1/track_data/$2/$2_unsorted.bed
# # echo "450k rows lost in hg19-->hg18 liftover out of 218M"

# #sort
sort -k1,1 -k2,2n data/$1/track_data/$2/$2_unsorted.bed > data/$1/track_data/$2/$2.bed 
wc -l data/$1/track_data/$2/$2.bed
rm data/$1/track_data/$2/$2_unsorted.bed

# bgzip and tabix 
bgzip data/$1/track_data/$2/$2.bed
tabix -p bed data/$1/track_data/$2/$2.bed.gz
