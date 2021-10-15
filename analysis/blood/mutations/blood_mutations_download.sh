###download and combine the mutation files for blood SC muts 

#wrangle hte manula file created by copying the DSMNC blood files 
awk -F, '{print $1}' analysis/blood/mutations/DSMNC_list_all_files_blood.csv > data/blood/mutations/blood_filesnamesOnly.txt #get only the col with the name . "," as the sep
tail -c +4 data/blood/mutations/blood_filesnamesOnly.txt > data/blood/mutations/blood_filesnamesOnly_noAnnoyingStratChar.txt #remove the invisible start character 

#wget 
cat data/blood/mutations/blood_filesnamesOnly_noAnnoyingStratChar.txt | while read line || [[ -n $line ]];
do 
    wget --no-check-certificate -q "https://dsmnc.big.ac.cn/DBtables/$line.txt" -P data/blood/mutations/
done 

#combine all files into 1 
touch data/blood/mutations/all_blood_mutations.txt #create the empty file 
head -n 1 data/blood/mutations/ID_1_individual_1_single-cell_2.txt > data/blood/mutations/all_blood_mutations.txt #insert the header first 
cat data/blood/mutations/blood_filesnamesOnly_noAnnoyingStratChar.txt | while read line || [[ -n $line ]]; #read and append each file into the vlank one we just made 
do 
    tail -n +2 "data/blood/mutations/$line.txt" >> data/blood/mutations/all_blood_mutations.txt #use tail -n +2 to remove the header of each subsequent file we append
done 


#data exploration 
# wc -l blood_filesnamesOnly_noAnnoyingStratChar.txt #checks the numberof fioles downloaded. should be 24 
# wc -l ../../../data/blood/mutations/all_blood_mutations.txt #checks the total # mutations in the combined file. should be 9171 
# cat ../../../data/blood/mutations/all_blood_mutations.txt | awk -F "\t" '{print $14}' | sort | uniq -c #checking the types of tissues used for single cell seq 
    # 2328 blood(bone marrow mononuclear cells,CD34+)
    # 3314 blood(hematopoietic stem/progenitor cells,HSPCs)
    # 3529 blood(mesenchymal stem cells,MSCs)
# cat ../../../data/blood/mutations/all_blood_mutations.txt | awk -F "\t" '{print $16}' | sort | uniq -c # checking the types of tissues used as controls 
    # 2328 blood(bone marrow mononuclear cells,CD34+)(bulk)
    # 3529 blood(mesenchymal stem cells,MSCs)(bulk)
    # 3232 bone marrow mononuclear cells
    # 82 peripheral blood cells
# cat ../../../data/blood/mutations/all_blood_mutations.txt | awk -F "\t" '{print $15}' | sort | uniq -c # checking the SC method 
    # 82 whole exome sequencing(single-stem-cell clonal culture)
    # 5857 whole genome sequencing(iPSC based single-cell clonal culture)
    # 3232 whole genome sequencing(single-stem-cell clonal culture)