#liftover and tabix blood mutations 
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg18.over.chain.gz -P data/blood/mutations/ #het the over chain file for liftover 
data/applications/UCSC_genomeBrowser_Blat/liftOver data/blood/mutations/all_blood_mutations_rearrangedFroLiftover.bed data/blood/mutations/hg19ToHg18.over.chain.gz data/blood/mutations/blood_mutations_hg18.bed data/blood/mutations/liftover_unlifted.bed #do the liftover 

#sort 
sort -k1,1 -k2,2n data/blood/mutations/blood_mutations_hg18.bed > data/blood/mutations/blood_mutations_hg18_sorted.bed

#return to tab delimited 
sed 's/|/\t/g' data/blood/mutations/blood_mutations_hg18_sorted.bed > data/blood/mutations/blood_mutations_hg18_sorted_tabdelim.bed
