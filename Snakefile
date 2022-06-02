tissues = ["skin","blood","liver","germline"]

rule all: 
    input: 
        #LIVER FILES
        "data/liver/dataframes/model2/predictorDf.txt",
        "analysis/liver/plots/model2/scatter_liver_on_liver.pdf",
        "data/liver/objects/model2/liver_model.RData",
        
        "data/liver/track_data/H3k4me1/H3k4me1.bed.gz",
        "data/liver/track_data/H3k4me3/H3k4me3.bed.gz",
        "data/liver/track_data/H3k27ac/H3k27ac.bed.gz",
        "data/liver/track_data/H3k36me3/H3k36me3.bed.gz",
        "data/liver/track_data/DNAse/DNAse.bed.gz",
        "data/liver/track_data/transcription/transcription.bed.gz",        
    
        #SKIN FILES 
        #"data/skin/dataframes/model2/predictorDf.txt",
        #"analysis/skin/plots/model2/scatter_skin_on_skin.pdf",
        #"data/skin/dataframes/model2/skin_on_skin_ProbabilityDf.csv",
        
        #GERMLINE FILES 
        "data/germline/dataframes/model2/predictorDf.txt",
        #"data/germline/dataframes/model2/germline_on_germline_ProbabilityDf.csv",
        #"analysis/germline/plots/model2/scatter_germline_on_germline.pdf",
        
        "data/germline/track_data/transcription/transcription_male_hg18_sorted.bed.gz",
        "data/germline/track_data/DNAse/DNAse_male_hg18_sorted.bed.gz",
        "data/germline/track_data/H3k27/H3k27ac_male_hg18_sorted.bed.gz",
            
        #BLOOD FILES 
        "data/blood/dataframes/model2/predictorDf.txt",
        "data/blood/dataframes/model2/blood_on_blood_ProbabilityDf.csv",
        #"analysis/blood/plots/model2/scatter_blood_on_blood.pdf",
        
        #"data/blood/mutations/blood_mutations_hg18_sorted_tabdelim.bed",
        "data/blood/track_data/DNAse/DNAse.bed.gz", 
        "data/blood/track_data/H3k27ac/H3k27ac.bed.gz", 
        "data/blood/track_data/H3k4me1/H3k4me1.bed.gz", 
        "data/blood/track_data/H3k4me3/H3k4me3.bed.gz", 
        "data/blood/track_data/transcription/transcription.bed.gz", 
        "data/blood/track_data/H3k27me3/H3k27me3.bed.gz",
        #"data/blood/dataframes/model2/predictorDf_2022_05_12.txt",
        
        #GLOBAL FILES
        "data/global/sequence/chr1.fa.gz",  
        "data/global/track_data/laminB1/laminB1_chr1.bed.gz", 
        "data/global/track_data/replication/replication.bed.gz",
        "data/global/track_data/recombination/recombination.bed.gz",
        "data/global/track_data/phastcons/phastcons_chr1.bed.gz",
        "data/global/track_data/repeats/repeats.bed.gz",
        
        #COLON FILES 
        #"data/colon/track_data/transcription/transcription.bed.gz",
        #"data/colon/track_data/H3k27ac/H3k27ac.bed.gz",
        #"data/colon/track_data/H3k27me3/H3k27me3.bed.gz",
        #"data/colon/track_data/H3k4me1/H3k4me1.bed.gz",
        #"data/colon/track_data/H3k4me3/H3k4me3.bed.gz",
        #"data/colon/track_data/H3k36me3/H3k36me3.bed.gz",
        #"data/colon/track_data/DNAse/DNAse.bed.gz"

#ANALYSIS RULES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rule createDF: 
    output: "data/{tissue}/dataframes/model2/predictorDf.txt"
    threads: 10
    conda: "conda_createDF.yml"
    shell: "python analysis/modules/createDF/createDF.py {wildcards.tissue} 'model2' '[1,100,10000]';"
        "grep 'buffer' data/{wildcards.tissue}/dataframes/model2/predictorDf.txt > data/{wildcards.tissue}/dataframes/model2/predictorDf_errorlog.txt;"
        "grep 'discord' data/{wildcards.tissue}/dataframes/model2/predictorDf.txt >> data/liver/dataframes/model2/predictorDf_errorlog.txt;"
        "grep -v 'discord' data/{wildcards.tissue}/dataframes/model2/predictorDf.txt > data/{wildcards.tissue}/dataframes/model2/predictorDf_noDiscord.txt;"
        "grep -v 'buffer' data/{wildcards.tissue}/dataframes/model2/predictorDf_noDiscord.txt >  data/{wildcards.tissue}/dataframes/model2/predictorDf.txt"

rule createModel: 
    input: "data/{tissue}/dataframes/{model}/predictorDf.txt"
    output:"data/{tissue}/objects/{model}/{tissue}_model.RData"
    conda: "conda_RcreateDfModel_env.yml"
    shell: "Rscript --vanilla analysis/modules/create_model/create_model.R {wildcards.tissue} {wildcards.model}"

rule predict: 
    input: "data/{tissue}/objects/{model}/{tissue}_model.RData"
    output: "data/{tissue}/dataframes/{model}/{tissue}_on_{tissue_predOn}_ProbabilityDf.csv"
    conda: "conda_RcreateDfModel_env.yml"
    shell: "Rscript --vanilla analysis/modules/create_model/predict.R {wildcards.tissue} {wildcards.tissue_predOn} {wildcards.model}"

rule plotting: 
    input: "data/{tissue}/dataframes/{model}/{tissue}_on_{tissue}_ProbabilityDf.csv"
    output: "analysis/{tissue}/plots/{model}/scatter_{tissue}_on_{tissue_predOn}.pdf"
    conda: "conda_Rplotting.yml"
    shell: "Rscript --vanilla analysis/modules/plotting_scatter/plotting_scatter.R {wildcards.model} 400 {wildcards.tissue} {wildcards.tissue_predOn}"    


#LIVER TRACK RULES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rule liver_H3k27ac_downloadWrangle:
    output: "data/liver/track_data/H3k27ac/H3k27ac.bed.gz"
    conda: "conda_snakeSomMut_env.yml"
    shell: "bash analysis/modules/download_encode/test.sh liver H3k27ac ENCFF764VSN.bigWig https://www.encodeproject.org/files/ENCFF764VSN/@@download/ENCFF764VSN.bigWig"

rule liver_H3k4me3_downloadWrangle:
    output: "data/liver/track_data/H3k4me3/H3k4me3.bed.gz"
    conda: "conda_snakeSomMut_env.yml"
    shell: "bash analysis/modules/download_encode/test.sh liver H3k4me3 ENCFF917LFF.bigWig https://www.encodeproject.org/files/ENCFF917LFF/@@download/ENCFF917LFF.bigWig"

rule liver_H3k4me1_downloadWrangle:
    output: "data/liver/track_data/H3k4me1/H3k4me1.bed.gz"
    conda: "conda_snakeSomMut_env.yml"
    shell: "bash analysis/modules/download_encode/test.sh liver H3k4me1 ENCFF042DZJ.bigWig https://www.encodeproject.org/files/ENCFF042DZJ/@@download/ENCFF042DZJ.bigWig"
        
rule liver_H3k36me3_downloadWrangle:
    output: "data/liver/track_data/H3k36me3/H3k36me3.bed.gz"
    conda: "conda_snakeSomMut_env.yml"
    shell: "bash analysis/modules/download_encode/test.sh liver H3k36me3 ENCFF376RWI.bigWig https://www.encodeproject.org/files/ENCFF376RWI/@@download/ENCFF376RWI.bigWig" 

rule liver_DNAse_downloadWrangle:
    output: "data/liver/track_data/DNAse/DNAse.bed.gz"
    conda: "conda_snakeSomMut_env.yml"
    shell: "bash analysis/modules/download_encode/test.sh liver DNAse ENCFF606YDZ.bigWig https://www.encodeproject.org/files/ENCFF606YDZ/@@download/ENCFF606YDZ.bigWig"

rule liver_transcription_downloadWrangle:
    output: "data/liver/track_data/transcription/transcription.bed.gz"
    conda: "conda_snakeSomMut_env.yml"
    shell: "bash analysis/modules/download_encode/test_twoFiles.sh liver transcription ENCFF186YCU ENCFF112IXW https://www.encodeproject.org/files/ENCFF186YCU/@@download/ENCFF186YCU.bigWig https://www.encodeproject.org/files/ENCFF112IXW/@@download/ENCFF112IXW.bigWig"
        
        
#SKIN TRACK RULES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





#GERMLINE TRACK RULES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        



rule germline_H3k27ac_download:
    output: "data/germline/track_data/H3k27/H3k27ac_male_hg18_sorted.bed.gz",
    shell: "bash analysis/germline/track_data/H3k27/germline_H3k27ac_downloadWrangle.sh" 

rule germline_dnase_download:
    output: "data/germline/track_data/DNAse/DNAse_male_hg18_sorted.bed.gz",
    shell: "bash analysis/germline/track_data/DNAse/germline_dnase_downloadWrangle.sh" 

rule germline_txn_download:
    output: "data/germline/track_data/transcription/transcription_male_hg18_sorted.bed.gz",
    shell: "bash analysis/germline/track_data/transcription/germline_transcription_downloadWrangle.sh" 


#BLOOD TRACK RULES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rule blood_createDF: #doesnt work yet 
    output: "data/blood/dataframes/model2/predictorDf_2022_05_12.txt"
    conda: "conda_snakeSomMut_env.yml"
    shell: "bash analysis/modules/createDF/createDF.py 'blood' 'model2' '[1,100,10000]' "
    
rule blood_mut_download:
    input: "analysis/blood/mutations/DSMNC_list_all_files_blood.csv"
    output: "data/blood/mutations/all_blood_mutations.txt"
    shell: "bash analysis/blood/mutations/blood_mutations_download.sh" 

rule blood_mut_prepareLiftover:
    input: "data/blood/mutations/all_blood_mutations.txt"
    output: "data/blood/mutations/all_blood_mutations_rearrangedFroLiftover.bed"
    shell: "python analysis/blood/mutations/blood_mutations_convert_to_bed.py"
        
rule blood_mut_liftoverSort:
    input: "data/blood/mutations/all_blood_mutations_rearrangedFroLiftover.bed"
    output: "data/blood/mutations/blood_mutations_hg18_sorted_tabdelim.bed"
    conda: "conda_snakeSomMut_env.yml"
    shell: 
        "bash analysis/blood/mutations/blood_mutations_liftoverSort.sh; "
        "grep 'chr3' data/blood/mutations/blood_mutations_hg18_sorted.bed > data/blood/mutations/test_chr3_muts_hg18.bed; "
        "python analysis/blood/mutations/blood_mutations_testing.py"

rule blood_methylation_downloadWrangle:
    output: "data/blood/track_data/methylation/methylation_CHG.bed.gz"
    conda: "conda_snakeSomMut_env.yml"
    shell: "bash analysis/blood/track_data/methylation/blood_methylation_downloadWrangle.sh"

rule blood_txn_downloadWrangle:
    output: "data/blood/track_data/transcription/transcription.bed.gz"
    message: "use the encode_download script for blood txn. may take a while"
    conda: "conda_snakeSomMut_env.yml"
    shell: "bash analysis/modules/download_encode/test.sh blood transcription ENCFF333QAU.bigWig https://www.encodeproject.org/files/ENCFF333QAU/@@download/ENCFF333QAU.bigWig"

rule blood_DNAse_downloadWrangle:
    output: "data/blood/track_data/DNAse/DNAse.bed.gz"
    message: "use the encode_download script for blood txn. may take a while"
    conda: "conda_snakeSomMut_env.yml"
    shell: "bash analysis/modules/download_encode/test.sh blood DNAse ENCFF850RIV.bigWig https://www.encodeproject.org/files/ENCFF850RIV/@@download/ENCFF850RIV.bigWig"       
        
rule blood_H3k27ac_downloadWrangle:
    output: "data/blood/track_data/H3k27ac/H3k27ac.bed.gz"
    message: "use the encode_download script for blood H3k27ac.may take a while (~4h?)"
    conda: "conda_snakeSomMut_env.yml"
    shell: "bash analysis/modules/download_encode/test.sh blood H3k27ac ENCFF311TAY.bigWig https://www.encodeproject.org/files/ENCFF311TAY/@@download/ENCFF311TAY.bigWig"
        
rule blood_H3k27me3_downloadWrangle:
    output: "data/blood/track_data/H3k27me3/H3k27me3.bed.gz"
    message: "use the encode_download script for blood H3k27me3.may take a while (~4h?)"
    conda: "conda_snakeSomMut_env.yml"
    shell: "bash analysis/modules/download_encode/test.sh blood H3k27me3 ENCFF810MPM.bigWig https://www.encodeproject.org/files/ENCFF810MPM/@@download/ENCFF810MPM.bigWig"

rule blood_H3k4me1_downloadWrangle:
    output: "data/blood/track_data/H3k4me1/H3k4me1.bed.gz"
    conda: "conda_snakeSomMut_env.yml"
    shell: "bash analysis/modules/download_encode/test.sh blood H3k4me1 ENCFF625RWJ.bigWig https://www.encodeproject.org/files/ENCFF625RWJ/@@download/ENCFF625RWJ.bigWig"
        
rule blood_H3k4me3_downloadWrangle:
    output: "data/blood/track_data/H3k4me3/H3k4me3.bed.gz"
    conda: "conda_snakeSomMut_env.yml"
    shell: "bash analysis/modules/download_encode/test.sh blood H3k4me3 ENCFF685DZI.bigWig https://www.encodeproject.org/files/ENCFF685DZI/@@download/ENCFF685DZI.bigWig"


#GLOBAL TRACK RULES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        
        
rule laminB1_downloadWrangle:
    output: "data/global/track_data/laminB1/laminB1_chr1.bed.gz"
    conda: "conda_snakeSomMut_env.yml"
    shell: "bash analysis/global/track_data/laminB1/global_laminB1_downloadWrangle.sh"

rule seq_download:
    output: "data/global/sequence/chr1.fa.gz"
    conda: "conda_snakeSomMut_env.yml"
    shell: "bash analysis/global/sequence/download_sequence.sh"

rule repeats_downloadWrangle:
    output: "data/global/track_data/repeats/repeats.bed.gz"
    conda: "conda_snakeSomMut_env.yml"
    shell: "bash analysis/global/track_data/repeats/global_repeats_downloadWrangle.sh"

rule phastcons_downloadWrangle:
    output: "data/global/track_data/phastcons/phastcons_chr1.bed.gz"
    conda: "conda_snakeSomMut_env.yml"
    shell: "bash analysis/global/track_data/phastcons/global_phastcons_downloadWrangle.sh"

rule recombination_downloadWrangle:
    output: "data/global/track_data/recombination/recombination.bed.gz"
    conda: "conda_snakeSomMut_env.yml"
    shell: "bash analysis/global/track_data/recombination/global_recombination_downloadWrangle.sh"
        
rule replication_downloadWrangle:
    output: "data/global/track_data/replication/replication.bed.gz"
    conda: "conda_snakeSomMut_env.yml"
    shell: "bash analysis/global/track_data/replication/global_replication_downloadWrangle.sh"
              
        
#COLON TRACK RULES~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~         

rule colon_histone_downloadWrangle:
    output:
         "data/colon/track_data/H3k27ac/H3k27ac.bed.gz",
         "data/colon/track_data/H3k27me3/H3k27me3.bed.gz",
         "data/colon/track_data/H3k4me1/H3k4me1.bed.gz",
         "data/colon/track_data/H3k4me3/H3k4me3.bed.gz",
         "data/colon/track_data/H3k36me3/H3k36me3.bed.gz",
    message: 
        "use the encode_download script for colon histones"
    conda: 
        "conda_snakeSomMut_env.yml"
    shell: 
        "bash analysis/modules/download_encode/test.sh colon H3k27ac ENCFF111DLN.bigWig https://www.encodeproject.org/files/ENCFF111DLN/@@download/ENCFF111DLN.bigWig"
        "bash analysis/modules/download_encode/test.sh colon H3k27me3"
        "bash analysis/modules/download_encode/test.sh colon H3k4me1 ENCFF623BIK.bigWig https://www.encodeproject.org/files/ENCFF623BIK/@@download/ENCFF623BIK.bigWig"
        "bash analysis/modules/download_encode/test.sh colon H3k4me3 ENCFF237VMY.bigWig https://www.encodeproject.org/files/ENCFF237VMY/@@download/ENCFF237VMY.bigWig"
        "bash analysis/modules/download_encode/test.sh colon H3k36me3 ENCFF165ZSZ.bigWig https://www.encodeproject.org/files/ENCFF165ZSZ/@@download/ENCFF165ZSZ.bigWig" 
        "bash analysis/modules/download_encode/test.sh colon DNAse ENCFF299OOV.bigWig https://www.encodeproject.org/files/ENCFF299OOV/@@download/ENCFF299OOV.bigWig"
        
rule colon_transcription_downloadWrangle:
    output: 
        "data/colon/track_data/transcription/transcription.bed.gz"
    conda: 
        "conda_snakeSomMut_env.yml"
    shell: 
        "bash analysis/modules/download_encode/test_twoFiles.sh colon transcription ENCFF304GQR.bigWig ENCFF036UAI.bigWig https://www.encodeproject.org/files/ENCFF304GQR/@@download/ENCFF304GQR.bigWig https://www.encodeproject.org/files/ENCFF036UAI/@@download/ENCFF036UAI.bigWig"


 