# imports ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import random
import pandas as pd
import numpy as np 
from numpy.random import choice
import collections
from Bio import AlignIO
import pysam 
from datetime import datetime
import gzip
import multiprocessing
import sys 
import json 
import time # for timing the loop 
#home-made modules : 
sys.path.append('/research/projects/hsapiens/mutability/analysis/global/track_data/annotation/') 
import annotation_handling
sys.path.append('/research/projects/hsapiens/mutability/analysis/modules/closest_value/') 
import closest_val


tmp_file_dir =""

# command line input ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

tissue = sys.argv[1]
# tissue = "liver"
model_desc = sys.argv[2]
# model_desc = "model4"
list_of_surrounding_contexts = json.loads(sys.argv[3]) #note if you chnagee/ increase this, then you oncrese the buffer zone (not using sites in the buffer, len(dna)-max_distance )
# list_of_surrounding_contexts = [1,100,10000]

if tissue == "germline": 
    mutations_lines = open(tmp_file_dir+'data/germline/mutation_data/mutations_hg18_final.bed').readlines()
    mutations_df = pd.read_table(tmp_file_dir+'data/germline/mutation_data/mutations_hg18_final.bed',sep="\t",header = None)
    mutations_df.columns = ["chromosome","start","fake_end","ref","alt","Fathers_age_at_conception","Mothers_age_at_conception"]
elif tissue in ["blood","liver","skin"]:
    mutations_lines = open(tmp_file_dir+'data/{t}/mutations/mutations.bed'.format(t=tissue)).readlines()
    mutations_df = pd.read_table(tmp_file_dir+'data/{t}/mutations/mutations.bed'.format(t=tissue),sep="\t",header = None)
    mutations_df.columns = ["chromosome","start","fake_end","ref","alt","ID","VAF","Gene name", "Region", "AA", "COSMIC", "Species", "Gender", "Age_in_years",           
                            "Tissue/Cell type","Single-cell_genomics_biotechnology_or_Method","Control_sample_or_tissue"]
else: 
    print(tissue," tissue specified not yet supported  !!~~~!!!~~~~!!")

    

#dictionry where i specify which col contains the information in the datafile , 0 indexed 
tracksColFile_dict = json.load(open(tmp_file_dir+"data/{t}/objects/{m}/tracksColDict.txt".format(m=model_desc,t=tissue)))#  


#mutant sites ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#create a dictionary where each chrom key will have a n empty list 
duplicate_lines = []
muts_bychrom_dict = {}
for x in range(1,23): 
    key_string = 'chr{n}'.format(n=x)
    muts_bychrom_dict[key_string] = []
    
print(tissue,": fill each chromosome's empty list  with the sites for that chrom")
non_chrnMuts = []#create list of chrom names that dont belong to chrN format --> disgnostic 
for line in (mutations_lines[1:]): 
    if line[0]=="c":                                               #aking sure the line is a chrN (lots of weird junk..) 
        chrom_mut = line.split("\t")[0]
        mut_startSite = line.split("\t")[1]                  #getting rid of the weird double(hgopefully) 
        if chrom_mut in muts_bychrom_dict.keys():                  #controlling for chrX/chrY
            if mut_startSite not in muts_bychrom_dict[chrom_mut]: 
                muts_bychrom_dict[chrom_mut].append(mut_startSite)
            else: duplicate_lines.append(line)
        else: 
            non_chrnMuts.append(chrom_mut)

#testing making usre the only sites that dont make it are sex chromosome mutations 
timestamp = datetime.now().strftime("%Y_%m_%d")
error_log = str("df created on "+timestamp)
error_log+=(str(len(non_chrnMuts))+"  non chrN muts (ommited) from these lables: "+str(list(np.unique(non_chrnMuts)))+"\n")

#add the sites infro from file 
sites = []#sites = list of sites 
for chrom_key in muts_bychrom_dict.keys(): 
    for mutation_element in muts_bychrom_dict[chrom_key]: 
        sites.append([chrom_key, int(mutation_element),1]) #the 1 is for mutation status column. 1 = yes 


#non -mutant sites ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#get chrom length information so I can perform weighted choice for non-mut site selection"
ChromLengths = pd.read_csv(tmp_file_dir+'data/global/sequence/hg38_chromosomelengths.csv') #read in the csv file of hg38 chrom lengths I found on the internets 
total_length=0 #lets sum (get the total length) 
for length in list(ChromLengths.Length): 
    total_length+=int(length.replace(",",""))

#build dictionary to store porbability 
dict_lengths = {}#creat emepty dictionary 
for x in range (0,22): 
        tmp_index = x +1
        length = str(ChromLengths[x:x+1]).split()[4]
        length = length.replace(",", "")
        length = int(length)
        dict_lengths["chr"+str(tmp_index)] = length

#make the porbability of choosing a chrom based on length 
list_chroms = ['chr' + str(i) for i in range(1, 23)]
list_chrom_probabilities = []
for chrom in list_chroms: 
    list_chrom_probabilities.append(dict_lengths[chrom]/total_length)
list_chrom_probabilities[0] = list_chrom_probabilities[0]+1-sum(list_chrom_probabilities) # adds the 0.00000001 left from rounding errors to the chr1 so sum adds perfectly to 1. 
assert(sum(list_chrom_probabilities)==1)

#perfrom the non-mutant site draw 
number_nonmuts = int(len(sites)*1.2)
chrom_draw = choice(list_chroms, number_nonmuts,p=list_chrom_probabilities)

print(tissue,"make the sites list with the chr# and site" )
for i in (range(1,23)): 
    chrom = "chr"+str(i)
    chrom_nchoose = list(chrom_draw).count("chr"+str(i))
    chrom_sites_chosen = random.sample(range(1, dict_lengths[chrom]), chrom_nchoose) #without duplucates 
    for j in chrom_sites_chosen: 
        sites.append([chrom,j,0])# the 0 if for the mutation status column. 0 = no 
        

# genral declarations before the big function ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

distance_max = max(list_of_surrounding_contexts)

fastas_dict = {}   # creating dictionary with fasta alignment, length of seq, 
print(tissue,"making the fastas dictionary")
for chrom in (list_chroms):
    filename_tmp = tmp_file_dir+"data/global/sequence/{c}.fa.gz".format(c=chrom)
    fastas_dict[chrom] = []
    with gzip.open(filename_tmp, "rt") as handle:
        fastas_dict[chrom].append(AlignIO.read(handle,"fasta"))
        alignment_tmp = fastas_dict[chrom][0]
        fastas_dict[chrom].append(len(str(alignment_tmp[0].seq)))

                 
#generate header ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        

header = "Chromosome"+ "\t"+"site"+"\t" +"triplet"+"\t"+"mutation_status"                        # creating the begining of the header 
for trackname in tracksColFile_dict.keys():   # the rest of the header is a function of tracks 
    if trackname in ["annotation","mappability","dist_rep_org_main","dist_rep_org_all","CpGisland"]: 
        header = header + "\t"+str(trackname)
    else: 
        for distance in list_of_surrounding_contexts:                                                # and distance (need a col for every track and for every distance value within ) 
            header = header + "\t"+str(trackname)+"-"+str(distance)
for distance in list_of_surrounding_contexts:                                                    # creating the end of the header assoc with no track (the seqeunce at different 
    header = header + "\t"+"Apercent-"+str(distance)+ "\t"+"Gpercent-"+str(distance)+ "\t"+"Cpercent-"+str(distance)+ "\t"+"Tpercent-"+str(distance)   # distace values) 
header = header +"\n"                                                                            # obviously needs to end with a \n 


filename = tmp_file_dir+'data/{a}/dataframes/{m}/predictorDf.txt'.format(a=tissue,m=model_desc)
with open(filename,"w") as f: 
    f.write(header)       
    


# big function ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def predictor_rowString(site): 
# site = sites[52429]
    row = []
    if site[1] <= distance_max or site[1]+distance_max >= fastas_dict[site[0]][1]:               # only use sites that will have values for site +- max distance (buffer). second element in fastas dict is the length 
        row.extend([str(site),"out of buffer range"])


    else: 
        row.extend([site[0], site[1]])
        alignment = fastas_dict[site[0]][0]                                                       #create the alingment from the list of fastas 

        #makes the triplet and tests alignment 
        if site[2] ==1:                                     #muts_by-chrom_dict is literally a dictionary containing list of sites that are mutations in that chrom. 
            mutation_row = mutations_df[(mutations_df.chromosome == site[0]) & (mutations_df.start == site[1])]  #get the row containing mut info out of the df 
            old_bp = mutation_row.ref.values[0]
            old_triplet = (str(alignment[0,site[1]-1])+str(old_bp)+str(alignment[0,site[1]+1])).upper()
            row.extend([old_triplet, 1])#usin the old/ref triplet in the df instead (the mutation happened to the old triplet! )

            seq_triplet = str(alignment[0,site[1]-1:site[1]+2].seq)
            if old_triplet.upper() != seq_triplet.upper():                                                          #testing that 
                row.append("discordant. triplet using daata = "+old_triplet+", seqeunce triplet = "+seq_triplet)
        else: 
            triplet= str(alignment[0,site[1]-1:site[1]+2].seq).upper()
            row.extend([triplet,0])

        for trackname,track_val in tracksColFile_dict.items():                 
            data_col = track_val[0] 
            global_or_tissue_specific = track_val[1]
            Na_is_0_or_NA = track_val[2]
            filename = tmp_file_dir+track_val[3]

            if trackname == "annotation": 
                if not [record for record in pysam.Tabixfile(filename).fetch(site[0], site[1], site[1]+1)]:                #if no value at that site 
                     row.append("not_transcribed")
                else: 
                    track_output = [record for record in pysam.Tabixfile(filename).fetch(site[0], site[1], site[1]+1)]
                    old_labels = [element.split()[2] for element in track_output]
                    converted_list, final_label,alien_labels = [],str(),[] 
                    for label in old_labels: 
                        if label in annotation_handling.annotation_conversion.keys():  #even though i though i controlled for it, occasiaonlyl there would eb anew label, so sontrol for this 
                            converted_list.append(annotation_handling.annotation_conversion[label])
                        else: 
                            alien_labels.append(label)
                    final_label = annotation_handling.annotation_priorityLabel(converted_list)
                    if len(alien_labels) != 0:    # if there is an "alien" annotation label, then just add that into the position and i can handle on ind basis later 
                        for label in alien_labels: 
                            final_label+="_"+label
                    row.append(final_label)
            elif trackname == "mappability": 
                if [record for record in pysam.Tabixfile(filename).fetch(site[0], site[1]-distance, site[1]+distance+1)]: row.append("mappable")
                else: row.append("not")
            elif trackname in ["dist_rep_org_main","dist_rep_org_all"]: 
                row.append(closest_val.shortest_distance(site,filename))
            elif trackname == "CpGisland": 
                if [record for record in pysam.Tabixfile(filename).fetch(site[0], site[1], site[1]+1)]: row.append("island")
                elif closest_val.shortest_distance(site,filename) <= 2000: row.append("shore")
                else: row.append("not")
            else: 
                for distance in list_of_surrounding_contexts: 
                    if not [record for record in pysam.Tabixfile(filename).fetch(site[0], site[1]-distance, site[1]+distance+1)]:                #if no value at that site 
                        if Na_is_0_or_NA == "Na=0": 
                            row.append(0)
                        else: 
                            row.extend(["NA"])                                                                                      
                    else:                                                   
                        track_output = [record for record in pysam.Tabixfile(filename).fetch(site[0], site[1]-distance, site[1]+distance+1)]
                        multiple_values = []

                        if tracksColFile_dict[trackname][0] == 4:
                            for element in track_output: 
                                multiple_values.append(float(element.split()[4]))
                            average_value = sum(multiple_values)/len(multiple_values) 
                            row.append(average_value)
                        elif tracksColFile_dict[trackname][0] == 3: 
                            for element in track_output: 
                                multiple_values.append(float(element.split()[3]))
                            average_value = sum(multiple_values)/len(multiple_values) 
                            row.append(average_value)
                        elif tracksColFile_dict[trackname][0] == 'binary': 
                            row.append(len(track_output))
                        else: 
                            error_log += ((str(site)+" ERROR: track coloumns not 4 or 5 or binary: "+trackname+"\n"))

        #sequence stuff                
        for distance in list_of_surrounding_contexts: 
            seq_around = str(alignment[0,site[1]-1:site[1]+2].seq)
            if seq_around != '': 
                seq_around = str(alignment[0,site[1]-distance:site[1]+distance+1].seq)
                Acount = seq_around.count('a')+seq_around.count("A")
                Gcount = seq_around.count('g')+seq_around.count("G")
                Ccount = seq_around.count('c')+seq_around.count("C")
                Tcount = seq_around.count('t')+seq_around.count("T")
                Apercent = Acount/len(seq_around)
                Gpercent = Gcount/len(seq_around)
                Cpercent = Ccount/len(seq_around)
                Tpercent = Tcount/len(seq_around)
                row.extend([Apercent, Gpercent, Cpercent, Tpercent])
            else: 
                row.extend(['NA','NA','NA','NA'])
                list_no_seq_at_site.append(site)

    row_string = str()
    for i in range(0,len(row)): 
        row_string = row_string+str(row[i])+"\t"
    row_string = row_string.rstrip("\t") # dont need to add the "\n" here as it is added below int he f.write 
#     row_string = row_string+"\n" 
    return row_string

#WRITE THE MODEL #
start_time = time.time()
print(tissue,"starting big loop")
def rowString_handler():
    p = multiprocessing.Pool(10)
    with open(filename, 'a') as f:
        for result in p.imap(predictor_rowString, sites[0:100]):
            f.write('%s\n' % result)

if __name__=='__main__':
    rowString_handler()
    
error_log += (("creating the df loop took "+str(time.time()-start_time)[0:4]+" seconds\n"))

#writing error log to file 
with open(tmp_file_dir+"data/{a}/dataframes/{m}/predictorDf_errorlog.txt".format(a=tissue,m=model_desc),"w") as f: 
      f.write(error_log)
     