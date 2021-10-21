### download and wranlge the male and femlae germline 

# - male germline ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# - downloaded from: 
#     - https://www.encodeproject.org/experiments/ENCSR755LFM/
#     - about the file: https://www.encodeproject.org/files/ENCFF193HUI/
# - download date: 
#     - lol it depends when i finally do the thing 
# - genome build 
#     - GRCh38 --> hg18 liftover 
# - tissue 
#     - testis 
#     - male embryo (maybe not best... they only start producing sperm later)
# - lab method:
#     - polyA plus RNA-seq
#     - general overview: https://www.encodeproject.org/documents/6354169f-86f6-4b59-8322-141005ea44eb/@@download/attachment/Long%20RNA-seq%20pipeline%20overview.pdf
# - file format & head: 
# - what the data represents: 
#     - signal of unique reads

#wget 
wget https://www.encodeproject.org/files/ENCFF193HUI/@@download/ENCFF193HUI.bigWig -P data/germline/track_data/transcription/

#convert from bigwig to bed 
data/applications/UCSC_genomeBrowser_Blat/bigWigToBedGraph data/germline/track_data/transcription/ENCFF193HUI.bigWig data/germline/track_data/transcription/transcription_male.bed

#prepr for liftovber: get the liftover chain files (n eed to do it in 2 steps, no hg38-->hg18 chain file exists....)
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz -P data/germline/track_data/transcription/
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg18.over.chain.gz -P data/germline/track_data/transcription/ 

#liftvoer hg38-->19 
data/applications/UCSC_genomeBrowser_Blat/liftOver data/germline/track_data/transcription/transcription_male.bed data/germline/track_data/transcription/hg38ToHg19.over.chain.gz data/germline/track_data/transcription/transcription_male_hg19.bed  data/germline/track_data/transcription/transcription_male_unliftgedhg38.bed
# wc -l ../../../../data/germline/track_data/transcription/transcription_male* 
echo "35k rows lost in male hg38-->hg19 liftover out of 26M"

#liftover hg19-->18 
data/applications/UCSC_genomeBrowser_Blat/liftOver data/germline/track_data/transcription/transcription_male_hg19.bed data/germline/track_data/transcription/hg19ToHg18.over.chain.gz data/germline/track_data/transcription/transcription_male_hg18.bed  data/germline/track_data/transcription/transcription_male_unliftgedhg19.bed
# wc -l data/germline/track_data/transcription/transcription_male* # 
echo "21k rows lost in male hg19-->hg18 liftover out of 26M"

#sort
sort -k1,1 -k2,2n data/germline/track_data/transcription/transcription_male_hg18.bed > data/germline/track_data/transcription/transcription_male_hg18_sorted.bed 

#bgzip and tabix 
bgzip data/germline/track_data/transcription/transcription_male_hg18_sorted.bed
tabix -p bed data/germline/track_data/transcription/transcription_male_hg18_sorted.bed.gz


# - female germline ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# - downloaded from: 
#     - https://www.encodeproject.org/experiments/ENCSR727VTD/
#     - https://www.encodeproject.org/files/ENCFF753LZJ/
# - download date: 
#     - lol it depends when i finally do the thing 
# - genome build 
#     - hg38 --> hg18 
# - tissue 
#     - Homo sapiens ovary tissue female embryo
# - lab method:
#     - polyA plus RNA-seq
#     - general details here: https://www.encodeproject.org/documents/6354169f-86f6-4b59-8322-141005ea44eb/@@download/attachment/Long%20RNA-seq%20pipeline%20overview.pdf
# - file format & head: 
#     - chr1	5080	5091	0.19661
#     - chr1	5091	5093	0.15729
#     - chr1	5093	5094	0.14418
#     - chr1	5094	5099	0.13108
# - what the data represents: 
#     - signal of unique reads

#wget 
wget https://www.encodeproject.org/files/ENCFF753LZJ/@@download/ENCFF753LZJ.bigWig -P  data/germline/track_data/transcription/

#convert from bigwig to bed 
data/applications/UCSC_genomeBrowser_Blat/bigWigToBedGraph /data/germline/track_data/transcription/ENCFF753LZJ.bigWig data/germline/track_data/transcription/transcription_female.bed

#liftvoer hg38-->19
data/applications/UCSC_genomeBrowser_Blat/liftOver data/germline/track_data/transcription/transcription_female.bed data/germline/track_data/transcription/hg38ToHg19.over.chain.gz data/germline/track_data/transcription/transcription_female_hg19.bed  data/germline/track_data/transcription/transcription_female_unliftgedhg38.bed
#### wc -l data/germline/track_data/transcription/transcription_female* 
echo "88k rows lost in female hg38-->hg19 liftover out of 42M"

#liftvoer hg19-->18 
data/applications/UCSC_genomeBrowser_Blat/liftOver data/germline/track_data/transcription/transcription_female_hg19.bed data/germline/track_data/transcription/hg19ToHg18.over.chain.gz data/germline/track_data/transcription/transcription_female_hg18.bed  data/germline/track_data/transcription/transcription_female_unliftgedhg19.bed
# wc -l data/germline/track_data/transcription/transcription_female* 
echo "43k rows lost in female hg19-->hg18 liftover out of 42M"

#sort 
sort -k1,1 -k2,2n data/germline/track_data/transcription/transcription_female_hg18.bed > data/germline/track_data/transcription/transcription_female_hg18_sorted.bed

#bzip and tabix 
bgzip data/germline/track_data/transcription/transcription_female_hg18_sorted.bed
tabix -p bed data/germline/track_data/transcription/transcription_female_hg18_sorted.bed.gz