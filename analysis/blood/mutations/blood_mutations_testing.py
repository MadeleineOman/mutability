#ensuring that the "ref" col of the mutations dataset matches the seq pos, aslo that the produced daatset is 0-based
import gzip 
from Bio import SeqIO

chr3_muts_lines = open("data/blood/mutations/test_chr3_muts_hg18.bed").readlines()
for line in chr3_muts_lines: 
    assert line.split("\t")[0] == "chr3"
misaligned_lines = []
count_n = 0
with gzip.open("data/global/sequence/chr3.fa.gz", "rt") as handle: #as seen in https://stackoverflow.com/questions/42757283/seqio-parse-on-a-fasta-gz
    for record in SeqIO.parse(handle,"fasta"): 
        for line in chr3_muts_lines: 
            mut_pos = int(line.split("\t")[2])-1
            mut_data_ref = line.split("\t")[3].split("|")[1].upper()
            hg18_seq_base = record.seq[mut_pos].upper()
            if hg18_seq_base != "N":
                assert mut_data_ref == hg18_seq_base