rule all: 
    input: 
        "data/blood/track_data/DNAse/DNAse_chr1.bed.gz", "data/blood/track_data/H3k27/H3k27_chr1.bed.gz", 
        "data/blood/track_data/H3k4me1/H3k4me1_chr1.bed.gz", "data/blood/track_data/H3k4me3/H3k4me3_chr1.bed.gz"
    
        
rule DNAse_downloadWrangle:
    output:
        "data/blood/track_data/DNAse/DNAse_chr1.bed.gz"
    message: 
        "use the DNAse download and wrangle script"
    conda: 
        "conda_snakeSomMut_env.yml"
    shell: 
        "bash analysis/blood/track_data/DNAse/blood_DNAse_downloadWrangle.sh"
        
        
rule H3k27_downloadWrangle:
    output:
        "data/blood/track_data/H3k27/H3k27_chr1.bed.gz"
    message: 
        "use the H3k27 download and wrangle script"
    conda: 
        "conda_snakeSomMut_env.yml"
    shell: 
        "bash analysis/blood/track_data/H3k27/blood_H3k27_downloadWrangle.sh"


rule H3k4me1_downloadWrangle:
    output:
        "data/blood/track_data/H3k4me1/H3k4me1_chr1.bed.gz"
    message: 
        "use the H3k4me1 download and wrangle script"
    conda: 
        "conda_snakeSomMut_env.yml"
    shell: 
        "bash analysis/blood/track_data/H3k4me1/blood_H3k4me1_downloadWrangle.sh"
        
rule H3k4me3_downloadWrangle:
    output:
        "data/blood/track_data/H3k4me3/H3k4me3_chr1.bed.gz"
    message: 
        "use the H3k4me3 download and wrangle script"
    conda: 
        "conda_snakeSomMut_env.yml"
    shell: 
        "bash analysis/blood/track_data/H3k4me3/blood_H3k4me3_downloadWrangle.sh"