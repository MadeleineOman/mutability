rule all: 
    input: 
        "data/blood/track_data/DNAse/DNAse_chr1.bed.gz", "data/blood/track_data/H3k27/H3k27_chr1.bed.gz", 
        "data/blood/track_data/H3k4me1/H3k4me1_chr1.bed.gz", "data/blood/track_data/H3k4me3/H3k4me3_chr1.bed.gz", 
        "data/global/track_data/laminB1/laminB1_chr1.bed.gz", "data/blood/track_data/transcription/transcription.bed.gz", 
        "data/global/track_data/replication/replication.bed.gz", "data/global/track_data/recombination/recombination.bed.gz", 
        "data/global/track_data/phastcons/phastcons_chr1.bed.gz"
    

rule phastcons_downloadWrangle:
    output:
        "data/global/track_data/phastcons/phastcons_chr1.bed.gz"
    message: 
        "use the phastcons download and wrangle script"
    conda: 
        "conda_snakeSomMut_env.yml"
    shell: 
        "bash analysis/global/track_data/phastcons/global_phastcons_downloadWrangle.sh"

rule recombination_downloadWrangle:
    output:
        "data/global/track_data/recombination/recombination.bed.gz"
    message: 
        "use the recombination download and wrangle script"
    conda: 
        "conda_snakeSomMut_env.yml"
    shell: 
        "bash analysis/global/track_data/recombination/global_recombination_downloadWrangle.sh"
        
rule replication_downloadWrangle:
    output:
        "data/global/track_data/replication/replication.bed.gz"
    message: 
        "use the replication download and wrangle script"
    conda: 
        "conda_snakeSomMut_env.yml"
    shell: 
        "bash analysis/global/track_data/replication/global_replication_downloadWrangle.sh"


rule txn_downloadWrangle:
    output:
        "data/blood/track_data/transcription/transcription.bed.gz"
    message: 
        "use the txn download and wrangle script"
    conda: 
        "conda_snakeSomMut_env.yml"
    shell: 
        "bash analysis/blood/track_data/transcription/blood_transcription_downloadWrangle.sh"

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
        
        
rule laminB1_downloadWrangle:
    output:
        "data/global/track_data/laminB1/laminB1_chr1.bed.gz"
    message: 
        "use the laminB1 download and wrangle script"
    conda: 
        "conda_snakeSomMut_env.yml"
    shell: 
        "bash analysis/global/track_data/laminB1/global_laminB1_downloadWrangle.sh"