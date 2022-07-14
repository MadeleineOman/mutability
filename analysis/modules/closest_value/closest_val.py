import pysam

def shortest_distance(site,filename):
    track_output = [record for record in pysam.Tabixfile(filename).fetch(site[0])]
    all_distances = []
    for instance in track_output: 
        chrom,start,end = instance.split("\t")[0:3]
        if int(start)<site[1]<int(end): 
            all_distances.append(0)
        else: 
            all_distances.extend([abs(int(start)-site[1]),abs(int(end)-site[1])])
    return(min(all_distances))