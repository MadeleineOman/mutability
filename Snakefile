rule all: 
    input: 
        "data/blood/track_data/DNAse/DNAse_chr1.bed.gz"
    
        
rule hello:
    output:
        "data/blood/track_data/DNAse/DNAse_chr1.bed.gz"
    message: 
        "use the DNAse download and wrangle script"
    conda: 
        "conda_mutability_env.yml"
    shell: 
        "bash analysis/blood/track_data/DNAse/blood_DNAse_downloadWrangle.sh"

