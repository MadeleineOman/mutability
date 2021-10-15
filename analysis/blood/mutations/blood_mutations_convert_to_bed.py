## formatting hte combined mut file into a nice bed format (preparing for liftiver) 


mut_lines = open("data/blood/mutations/all_blood_mutations.txt").readlines()

with open("data/blood/mutations/all_blood_mutations_rearrangedFroLiftover.bed","w") as f:
    for line in mut_lines[1:]: 
        chrom = line.split("\t")[1]
        start = int(line.split("\t")[2])-1
        stop = start+1                            #add a "stop" line 
        file_id = line.split("\t")[0]
        row_string = chrom+"\t"+str(start)+"\t"+str(stop)+"\t"+file_id #rest of the cols are automated 
        for i in line.split("\t")[3:]: 
            i = i.replace(" ","_")   #get ride of all spaces 
            row_string+="|"+str(i) #thought to add | as the sep for the later fields as done in this conversion vcf to bed format https://genome.soe.ucsc.narkive.com/g5oNkLn5/liftover-hg19-to-hg18
        f.write(row_string) 