import gzip
import numpy as np

def digest_GFF(gff):
    with gzip.open(gff, 'rb') as f:
        tss_data = {}
        bc=0
        for line in f:
            feature_type = line.strip().split()[2].decode('UTF-8')
            if feature_type != "transcript": continue
            chromosome = line.strip().split()[0].decode('UTF-8')
            start = int(line.strip().split()[3])
            end = int(line.strip().split()[4])
            strand = line.strip().split()[6].decode('UTF-8')
            if strand == "-": tss = end
            else: tss = start
            if chromosome not in tss_data:
                tss_data[chromosome] = {}
            tss_data[chromosome][tss] = {"strand":strand, "start":start, "end":end}
    return tss_data


def dist2nearest_tss(c, p,tss_data):
    tss_list = list(tss_data[c].keys())
    closest_idx, closest_tss = min(enumerate(tss_list), key=lambda x: abs(x[1]-p))
    strand = tss_data[c][closest_tss]["strand"]
    # after the last gene
    if closest_idx == len(tss_list) and closest_tss < p:
        if strand == "-": dist_to_tss = closest_tss-p
        elif strand == "+": dist_to_tss = p - closest_tss
        return  dist_to_tss
    #before the first gene
    elif closest_idx == 0 and closest_tss > p : 
        if   strand == "-": dist_to_tss = p-closest_tss
        elif strand == "+": dist_to_tss = closest_tss-p
        return  dist_to_tss
    # right on a TSS
    elif closest_tss - p == 0 : 
        #tss is the mutation
        dist_to_tss = 0
        return dist_to_tss
    #otherwise its in the middle of two genes so lets define the left and right
    # before continuing
    elif closest_tss - p < 0: # if closest is left of p
        right_idx, right_tss =  closest_idx+1, tss_list[closest_idx+1]
        left_idx,  left_tss = closest_idx, closest_tss
    elif closest_tss - p > 0:    #if closest is right of p
        right_idx, right_tss =closest_idx, closest_tss
        left_idx,  left_tss = closest_idx-1, tss_list[closest_idx-1]
    else:
        print("wtf? what else is possible??")
        return None

    left_strand  = tss_data[c][left_tss]["strand"]
    right_strand = tss_data[c][right_tss]["strand"]
    #####################################################################
    ### There are 4 scenarios of what genes p can lie between:
    #    1  -/+
    #    2  +/+
    #    3  -/-
    #    4  +/-

    # (1) -  < === p === > +
    # p must be intergenic and the upstream of both genes
    # it is therefore a negative value upstream of the nearer gene
    if left_strand == "-" and right_strand == "+":
        dist_to_tss = -1 * min([abs(p-left_tss),abs(p-right_tss)])
        return dist_to_tss
    # (2) +  > === p === > +
    # p could be in left gene, downstream of left or upstream of right 
    elif left_strand == "+" and right_strand == "+":
        left_end = tss_data[c][left_tss]["end"]
        if left_tss <= p <= left_end: # in the gene
            dist_to_tss = p-left_tss #positive
        else:
            if  abs(left_end-p) < abs(right_tss-p): #closer to left gene
                dist_to_tss = p-left_tss #positive
            elif abs(left_end-p) >= abs(right_tss-p): #closer to right tss
                dist_to_tss = p-right_tss #negative
            else:
                print("wtf")
        return dist_to_tss
    # (3) -  < === p === < -
    # p could be in right gene, downstream of right or upstream of left   
    elif left_strand == "-" and right_strand == "-":
        right_end = tss_data[c][left_tss]["start"]
        if right_end <= p <= right_tss: # in the gene
            dist_to_tss = right_tss-p #positive
        else:
            if abs(right_end-p) < abs(left_tss-p): #closer to right gene
                dist_to_tss = right_tss-p #positive
            elif abs(right_end-p) >=  abs(left_tss-p): #closer to left tss
                dist_to_tss = left_tss-p #negative
            else:print("WTF")
        return dist_to_tss
    # (4) +  > === p === < -
    # p could be in left gene, in the right gene or in between
    elif left_strand == "+" and right_strand == "-":
        left_end = tss_data[c][left_tss]["end"]
        right_end = tss_data[c][left_tss]["start"]
        if left_tss <= p <= left_end: # in the left gene
            dist_to_tss = p-left_tss #postive
        elif right_end <= p <= right_tss: # in the right gene
            dist_to_tss = right_tss - p #postive
        else: #in between
            dist_to_tss = min([abs(p-left_end),abs(p-right_end)])
        return dist_to_tss
    else:
        print("hmmm?")
        return None







def main():
    gff = "/research/projects/hsapiens/mutability/data/global/track_data/annotation/global_annotation_sorted.gff.gz"# Digest gff into a dictionary
    tss_data = digest_GFF(gff)
    #### run it on a few sites
    ## @Madeleine - you need to figure out how you want to feed mutations to this program. 
    # It is slow to call the script once per mutation because it will digest the GFF everytime
    #   you should digest the GFF once and then run it on many mutations.
    # You just need to replace c and p below with the chromosome and position, 
    #   where the position should be an integer and chromosomes are like the gff (chr1, chr4, chr12 etc)
    c = 'chr12'
    for p in np.random.randint(1, 100000000, 100):
        dist_to_tss =  dist2nearest_tss(c, p,tss_data)
        print(dist_to_tss)

if __name__ == "__main__":
    main()

