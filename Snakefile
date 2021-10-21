rule all: 
    input: 
        "data/blood/track_data/DNAse/DNAse_chr1.bed.gz", "data/blood/track_data/H3k27/H3k27_chr1.bed.gz", 
        "data/blood/track_data/H3k4me1/H3k4me1_chr1.bed.gz", "data/blood/track_data/H3k4me3/H3k4me3_chr1.bed.gz", 
         "data/blood/track_data/transcription/transcription.bed.gz", "data/blood/mutations/blood_mutations_hg18_sorted_tabdelim.bed",
          
         "data/global/sequence/chr1.fa.gz",  
         "data/global/track_data/laminB1/laminB1_chr1.bed.gz", "data/global/track_data/replication/replication.bed.gz",
         "data/global/track_data/recombination/recombination.bed.gz","data/global/track_data/phastcons/phastcons_chr1.bed.gz",
         "data/global/track_data/repeats/repeats.bed.gz",
         
         "data/germline/track_data/transcription/transcription_male_hg18_sorted.bed.gz",
         "data/germline/track_data/transcription/transcription_female_hg18_sorted.bed.gz"
         
         
    

rule germline_txn_download:
    output:
        "data/germline/track_data/transcription/transcription_male_hg18_sorted.bed.gz",
        "data/germline/track_data/transcription/transcription_female_hg18_sorted.bed.gz"
    message: 
        "use the ownlaod wrangle for germline"
    shell: 
        "bash analysis/germline/track_data/transcription/germline_transcription_downloadWrangle.sh" 


rule blood_mut_download:
    input: 
        "analysis/blood/mutations/DSMNC_list_all_files_blood.csv"
    output:
        "data/blood/mutations/all_blood_mutations.txt"
    message: 
        "use the ownlaod wrangle and test the blood mutations"
    shell: 
        "bash analysis/blood/mutations/blood_mutations_download.sh" 

rule blood_mut_prepareLiftover:
    input: 
         "data/blood/mutations/all_blood_mutations.txt"
    output:
        "data/blood/mutations/all_blood_mutations_rearrangedFroLiftover.bed"
    message: 
        "use the python script to prepare the blood mutations file for liftover"
    shell: 
        "python analysis/blood/mutations/blood_mutations_convert_to_bed.py"
        
rule blood_mut_liftoverSort:
    input: 
        "data/blood/mutations/all_blood_mutations_rearrangedFroLiftover.bed"
    output:
        "data/blood/mutations/blood_mutations_hg18_sorted_tabdelim.bed"
    conda: 
        "conda_snakeSomMut_env.yml"
    message: 
        "use script to liftover and sort. then oython script to tets mutations"
    shell: 
        "bash analysis/blood/mutations/blood_mutations_liftoverSort.sh; "
        "grep 'chr3' data/blood/mutations/blood_mutations_hg18_sorted.bed > data/blood/mutations/test_chr3_muts_hg18.bed; "
        "python analysis/blood/mutations/blood_mutations_testing.py"


rule seq_download:
    output:
        "data/global/sequence/chr1.fa.gz"
    message: 
        "use the hg18 seq download script"
    conda: 
        "conda_snakeSomMut_env.yml"
    shell: 
        "bash analysis/global/sequence/download_sequence.sh"

rule blood_methylation_downloadWrangle:
    output:
        "data/blood/track_data/methylation/methylation_CHG.bed.gz"
    message: 
        "use the blood methylation download and wrangle script"
    conda: 
        "conda_snakeSomMut_env.yml"
    shell: 
        "bash analysis/blood/track_data/methylation/blood_methylation_downloadWrangle.sh"

rule repeats_downloadWrangle:
    output:
        "data/global/track_data/repeats/repeats.bed.gz"
    message: 
        "use the repeats download and wrangle script"
    conda: 
        "conda_snakeSomMut_env.yml"
    shell: 
        "bash analysis/global/track_data/repeats/global_repeats_downloadWrangle.sh"

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


rule blood_txn_downloadWrangle:
    output:
        "data/blood/track_data/transcription/transcription.bed.gz"
    message: 
        "use the txn download and wrangle script"
    conda: 
        "conda_snakeSomMut_env.yml"
    shell: 
        "bash analysis/blood/track_data/transcription/blood_transcription_downloadWrangle.sh"

rule blood_DNAse_downloadWrangle:
    output:
        "data/blood/track_data/DNAse/DNAse_chr1.bed.gz"
    message: 
        "use the DNAse download and wrangle script"
    conda: 
        "conda_snakeSomMut_env.yml"
    shell: 
        "bash analysis/blood/track_data/DNAse/blood_DNAse_downloadWrangle.sh"
        
        
rule blood_H3k27_downloadWrangle:
    output:
        "data/blood/track_data/H3k27/H3k27_chr1.bed.gz"
    message: 
        "use the H3k27 download and wrangle script"
    conda: 
        "conda_snakeSomMut_env.yml"
    shell: 
        "bash analysis/blood/track_data/H3k27/blood_H3k27_downloadWrangle.sh"


rule blood_H3k4me1_downloadWrangle:
    output:
        "data/blood/track_data/H3k4me1/H3k4me1_chr1.bed.gz"
    message: 
        "use the H3k4me1 download and wrangle script"
    conda: 
        "conda_snakeSomMut_env.yml"
    shell: 
        "bash analysis/blood/track_data/H3k4me1/blood_H3k4me1_downloadWrangle.sh"
        
rule blood_H3k4me3_downloadWrangle:
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