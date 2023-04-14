
rule all: 
    input:              
        #coef plot 
        "analysis/global/plots/model8/coefScatter_blood_on_germline_onlySignCoefs_equiv_toLowest_comb.pdf",
        
        #remake the model for the coef plotting --> equiv
        #"data/blood/dataframes/model8/blood_coefDF_equiv_toLowest.csv",
        #"data/germline/dataframes/model8/germline_coefDF_equiv_toLowest.csv",
        #"data/liver/dataframes/model8/liver_coefDF_equiv_toLowest.csv",
        #"data/skin/dataframes/model8/skin_coefDF.csv",
        
        #need the model for the prediction --> not equiv to lowest and need the r data 
        "data/blood/objects/model8/blood_model_fullModel.RData",
        #"data/germline/objects/model8/germline_model_fullModel.RData",
        #"data/liver/objects/model8/liver_model_fullModel.RData",
        #"data/skin/objects/model8/skin_model_fullModel.RData",
        
        #scatter plot 
        "analysis/blood/plots/model8/scatter_blood_on_blood_fullModel.pdf",
        #"analysis/germline/plots/model8/scatter_germline_on_germline_fullModel.pdf",
        #"analysis/liver/plots/model8/scatter_liver_on_liver_fullModel.pdf",
        #"analysis/skin/plots/model8/scatter_skin_on_skin_fullModel.pdf",
        

        #modelprep
        #"data/blood/dataframes/model8/blood_all_data_readyForPrediction.csv",
        #"data/germline/dataframes/model8/germline_all_data_readyForPrediction.csv",
        #"data/liver/dataframes/model8/liver_all_data_readyForPrediction.csv",
        #"data/skin/dataframes/model8/skin_all_data_readyForPrediction.csv",

        #new df 
        #"data/blood/dataframes/model7/predictorDf_allTissueSpecTracks.txt",
        #"data/germline/dataframes/model7/predictorDf_allTissueSpecTracks.txt",
        #"data/liver/dataframes/model7/predictorDf_allTissueSpecTracks.txt",
        #"data/skin/dataframes/model7/predictorDf_allTissueSpecTracks.txt",

        #SCATTER PLOT PREDICTION ON SELF        
        #"analysis/blood/plots/model6/scatter_blood_on_blood.pdf",
        #"analysis/liver/plots/model6/scatter_liver_on_liver.pdf",
        #"analysis/germline/plots/model6/scatter_germline_on_germline.pdf",
        #"analysis/skin/plots/model6/scatter_skin_on_skin.pdf",
        
        #SCATTER PLOT PREDICTION ON OTHER 
        #"analysis/skin/plots/model6/scatter_skin_on_blood.pdf",
        #"analysis/skin/plots/model6/scatter_skin_on_germline.pdf",
        #"analysis/skin/plots/model6/scatter_skin_on_liver.pdf",
        #"analysis/blood/plots/model6/scatter_blood_on_skin.pdf",
        #"analysis/blood/plots/model6/scatter_blood_on_germline.pdf",
        #"analysis/blood/plots/model6/scatter_blood_on_liver.pdf",
        #"analysis/germline/plots/model6/scatter_germline_on_skin.pdf",
        #"analysis/germline/plots/model6/scatter_germline_on_blood.pdf",
        #"analysis/germline/plots/model6/scatter_germline_on_liver.pdf",
        #"analysis/germline/plots/model6/scatter_germline_on_liver.pdf",
        #"analysis/liver/plots/model6/scatter_liver_on_blood.pdf",
        #"analysis/liver/plots/model6/scatter_liver_on_skin.pdf",
        #"analysis/liver/plots/model6/scatter_liver_on_germline.pdf",
        
        #COEF COMPARISON SCATTER PLOT 
        #"analysis/global/plots/model6/coefScatter_blood_on_skin.pdf",
        #"analysis/global/plots/model6/coefScatter_blood_on_germline.pdf",
        #"analysis/global/plots/model6/coefScatter_blood_on_liver.pdf",
        #"analysis/global/plots/model6/coefScatter_germline_on_skin.pdf",
        #"analysis/global/plots/model6/coefScatter_germline_on_skin.pdf",
        #"analysis/global/plots/model6/coefScatter_germline_on_liver.pdf",
        #"analysis/global/plots/model6/coefScatter_skin_on_liver.pdf",
        
        #PCA 
        #"analysis/global/plots/model4/pca_allData_mutsAndNonMuts.pdf"

#ANALYSIS RULES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rule dataPrep:
    input: "data/{tissue}/dataframes/{model}/predictorDf.txt"
    output:"data/{tissue}/dataframes/{model}/{tissue}_all_data_readyForPrediction{model_desc}.csv"
    conda: "conda_RcreateDfModel_env.yml"
    shell: "Rscript --vanilla analysis/modules/create_model/dataWrangle_modelPrep.R {wildcards.tissue} {wildcards.model} {wildcards.model_desc}" # 

rule createModel: 
    input: "data/{tissue}/dataframes/{model}/{tissue}_all_data_readyForPrediction{model_desc}.csv"
    output:"data/{tissue}/dataframes/{model}/{tissue}_coefDF{model_desc}.csv",
           "data/{tissue}/objects/{model}/{tissue}_model{model_desc}.RData"
    conda: "conda_RcreateDfModel_env.yml"
    shell: "Rscript --vanilla analysis/modules/create_model/create_linearModel.R {wildcards.tissue} {wildcards.model} 0 {wildcards.model_desc}" #tissue, model, n samples

rule predict: 
    input: "data/{tissue}/objects/{model}/{tissue}_model{model_desc}.RData"
    output: "data/{tissue}/dataframes/{model}/{tissue}_on_{tissue_predOn}_ProbabilityDf{model_desc}.csv"
    conda: "conda_RcreateDfModel_env.yml"
    shell: "Rscript --vanilla analysis/modules/create_model/predict.R {wildcards.tissue} {wildcards.tissue_predOn} {wildcards.model} {wildcards.model_desc}" 
 
rule plotting_scatter: 
    input: "data/{tissue}/dataframes/{model}/{tissue}_on_{tissue_predOn}_ProbabilityDf{model_desc}.csv"
    output: "analysis/{tissue}/plots/{model}/scatter_{tissue}_on_{tissue_predOn,[A-Za-z]+}{model_desc}.pdf" #the ,[A-Za-z]+ enforces only match letters
    conda: "conda_Rplotting.yml"
    shell: "Rscript --vanilla analysis/modules/plotting_scatter/plotting_scatter.R {wildcards.model} 400 {wildcards.tissue} {wildcards.tissue_predOn} {wildcards.model_desc}"  

rule plotting_coef: 
    input: "data/{tissue}/dataframes/{model}/{tissue}_coefDF{model_desc}.csv",
           "data/{tissue_predOn}/dataframes/{model}/{tissue_predOn}_coefDF{model_desc}.csv"
    output: "analysis/global/plots/{model}/coefScatter_{tissue}_on_{tissue_predOn}_onlySignCoefs{model_desc}_comb.pdf"
    conda: "conda_Rplotting.yml"
    shell: "Rscript --vanilla analysis/modules/plotting_coef/plotting_coef_pairwise.R {wildcards.tissue} {wildcards.tissue_predOn} {wildcards.model} {wildcards.model_desc}"   

rule pca: 
    input: 
        "data/blood/dataframes/{model}/blood_all_data_readyForPrediction.csv", 
        "data/liver/dataframes/{model}/liver_all_data_readyForPrediction.csv",
        "data/germline/dataframes/{model}/germline_all_data_readyForPrediction.csv",
        "data/skin/dataframes/{model}/skin_all_data_readyForPrediction.csv"
    output: "analysis/global/plots/{model}/pca_allData_mutsAndNonMuts.pdf"
    conda: "conda_Rplotting.yml"
    shell: "Rscript --vanilla analysis/modules/pca/pca.R {wildcards.model}" 



#CREATE DF RULES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rule createDF_blood: 
    input: 
        #"data/global/sequence/chr1.fa.gz",  
        "data/global/track_data/laminB1/laminB1_chr1.bed.gz", 
        "data/global/track_data/replication/replication.bed.gz",
        "data/global/track_data/recombination/recombination.bed.gz",
        #"data/global/track_data/phastcons/phastcons_chr1.bed.gz",
        "data/global/track_data/repeats/repeats.bed.gz",
        "data/global/track_data/mappability/mappability.bed.gz",
        "data/blood/track_data/DNAse/DNAse.bed.gz", 
        "data/blood/track_data/H3k27ac/H3k27ac.bed.gz", 
        "data/blood/track_data/H3k4me1/H3k4me1.bed.gz", 
        "data/blood/track_data/H3k4me3/H3k4me3.bed.gz", 
        "data/blood/track_data/transcription/transcription.bed.gz", 
        "data/blood/track_data/H3k27me3/H3k27me3.bed.gz",
        "data/blood/track_data/methylation/methylation.bed.gz",
    output: "data/blood/dataframes/{model}/predictorDf{model_desc}.txt"
    threads: 10
    conda: "conda_createDF.yml"
    shell: "python analysis/modules/createDF/createDF.py blood {wildcards.model} '[0,100,10000]' {wildcards.model_desc};"
        "grep 'buffer' data/blood/dataframes/{wildcards.model}/predictorDf{wildcards.model_desc}.txt >> data/blood/dataframes/{wildcards.model}/predictorDf{wildcards.model_desc}_errorlog.txt;"
        "grep 'discord' data/blood/dataframes/{wildcards.model}/predictorDf{wildcards.model_desc}.txt >> data/blood/dataframes/{wildcards.model}/predictorDf{wildcards.model_desc}_errorlog.txt;"
        "grep -v 'discord' data/blood/dataframes/{wildcards.model}/predictorDf{wildcards.model_desc}.txt > data/blood/dataframes/{wildcards.model}/predictorDf{wildcards.model_desc}_noDiscord.txt;"
        "grep -v 'buffer' data/blood/dataframes/{wildcards.model}/predictorDf{wildcards.model_desc}_noDiscord.txt >  data/blood/dataframes/{wildcards.model}/predictorDf{wildcards.model_desc}.txt"
        
rule createDF_liver: 
    input:  
        #"data/global/sequence/chr1.fa.gz",  
        "data/global/track_data/laminB1/laminB1_chr1.bed.gz", 
        "data/global/track_data/replication/replication.bed.gz",
        "data/global/track_data/recombination/recombination.bed.gz",
        #"data/global/track_data/phastcons/phastcons_chr1.bed.gz",
        "data/global/track_data/repeats/repeats.bed.gz",
        "data/global/track_data/mappability/mappability.bed.gz",
        "data/liver/track_data/H3k4me1/H3k4me1.bed.gz",
        "data/liver/track_data/H3k4me3/H3k4me3.bed.gz",
        "data/liver/track_data/H3k27ac/H3k27ac.bed.gz",
        "data/liver/track_data/H3k36me3/H3k36me3.bed.gz",
        "data/liver/track_data/DNAse/DNAse.bed.gz",
        "data/liver/track_data/transcription/transcription.bed.gz", 
        "data/liver/track_data/methylation/methylation.bed.gz",
    output: "data/liver/dataframes/{model}/predictorDf{model_desc}.txt"
    threads: 10
    conda: "conda_createDF.yml"
    shell: "python analysis/modules/createDF/createDF.py liver {wildcards.model} '[0,100,10000]' {wildcards.model_desc};"
        "grep 'buffer' data/liver/dataframes/{wildcards.model}/predictorDf{wildcards.model_desc}.txt >> data/liver/dataframes/{wildcards.model}/predictorDf{wildcards.model_desc}_errorlog.txt;"
        "grep 'discord' data/liver/dataframes/{wildcards.model}/predictorDf{wildcards.model_desc}.txt >> data/liver/dataframes/{wildcards.model}/predictorDf{wildcards.model_desc}_errorlog.txt;"
        "grep -v 'discord' data/liver/dataframes/{wildcards.model}/predictorDf{wildcards.model_desc}.txt > data/liver/dataframes/{wildcards.model}/predictorDf{wildcards.model_desc}_noDiscord.txt;"
        "grep -v 'buffer' data/liver/dataframes/{wildcards.model}/predictorDf{wildcards.model_desc}_noDiscord.txt >  data/liver/dataframes/{wildcards.model}/predictorDf{wildcards.model_desc}.txt"

        
rule createDF_germline: 
    input:  
        #"data/global/sequence/chr1.fa.gz",  
        "data/global/track_data/laminB1/laminB1_chr1.bed.gz", 
        "data/global/track_data/replication/replication.bed.gz",
        "data/global/track_data/recombination/recombination.bed.gz",
        #"data/global/track_data/phastcons/phastcons_chr1.bed.gz",
        "data/global/track_data/repeats/repeats.bed.gz",
        "data/global/track_data/mappability/mappability.bed.gz",
        "data/germline/track_data/transcription/transcription_male_hg18_sorted.bed.gz",
        "data/germline/track_data/DNAse/DNAse_male_hg18_sorted.bed.gz",
        "data/germline/track_data/H3k27/H3k27ac_male_hg18_sorted.bed.gz",
        "data/germline/track_data/methylation/methylation.bed.gz",
    output: "data/germline/dataframes/{model}/predictorDf{model_desc}.txt"
    threads: 10
    conda: "conda_createDF.yml"
    shell: "python analysis/modules/createDF/createDF.py germline {wildcards.model} '[0,100,10000]' {wildcards.model_desc};"
        "grep 'buffer' data/germline/dataframes/{wildcards.model}/predictorDf{wildcards.model_desc}.txt >> data/germline/dataframes/{wildcards.model}/predictorDf{wildcards.model_desc}_errorlog.txt;"
        "grep 'discord' data/germline/dataframes/{wildcards.model}/predictorDf{wildcards.model_desc}.txt >> data/germline/dataframes/{wildcards.model}/predictorDf{wildcards.model_desc}_errorlog.txt;"
        "grep -v 'discord' data/germline/dataframes/{wildcards.model}/predictorDf{wildcards.model_desc}.txt > data/germline/dataframes/{wildcards.model}/predictorDf{wildcards.model_desc}_noDiscord.txt;"
        "grep -v 'buffer' data/germline/dataframes/{wildcards.model}/predictorDf{wildcards.model_desc}_noDiscord.txt >  data/germline/dataframes/{wildcards.model}/predictorDf{wildcards.model_desc}.txt"

        
rule createDF_skin: 
    input:   
        #"data/global/sequence/chr1.fa.gz",  
        "data/global/track_data/laminB1/laminB1_chr1.bed.gz", 
        "data/global/track_data/replication/replication.bed.gz",
        "data/global/track_data/recombination/recombination.bed.gz",
        #"data/global/track_data/phastcons/phastcons_chr1.bed.gz",
        "data/global/track_data/repeats/repeats.bed.gz",
        "data/global/track_data/mappability/mappability.bed.gz",
        "data/skin/track_data/methylation/methylation.bed.gz",
    output: "data/skin/dataframes/{model}/predictorDf{model_desc}.txt"
    threads: 10
    conda: "conda_createDF.yml"
    shell: "python analysis/modules/createDF/createDF.py skin {wildcards.model} '[0,100,10000]' {wildcards.model_desc};"
        "grep 'buffer' data/skin/dataframes/{wildcards.model}/predictorDf{wildcards.model_desc}.txt >> data/skin/dataframes/{wildcards.model}/predictorDf{wildcards.model_desc}_errorlog.txt;"
        "grep 'discord' data/skin/dataframes/{wildcards.model}/predictorDf{wildcards.model_desc}.txt >> data/skin/dataframes/{wildcards.model}/predictorDf{wildcards.model_desc}_errorlog.txt;"
        "grep -v 'discord' data/skin/dataframes/{wildcards.model}/predictorDf{wildcards.model_desc}.txt > data/skin/dataframes/{wildcards.model}/predictorDf{wildcards.model_desc}_noDiscord.txt;"
        "grep -v 'buffer' data/skin/dataframes/{wildcards.model}/predictorDf{wildcards.model_desc}_noDiscord.txt >  data/skin/dataframes/{wildcards.model}/predictorDf{wildcards.model_desc}.txt"




#LIVER TRACK RULES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rule liver_ctcf_downloadWrangle:
    output: "data/liver/track_data/ctcf/ctcf.bed.gz"
    conda: "conda_snakeSomMut_env.yml"
    shell: "bash analysis/modules/download_encode/test.sh liver ctcf ENCFF005YBS.bigWig https://www.encodeproject.org/files/ENCFF005YBS/@@download/ENCFF005YBS.bigWig"

rule liver_methylation_downloadWrangle:
    output: "data/liver/track_data/methylation/methylation.bed.gz"
    conda: "conda_snakeSomMut_env.yml"
    shell: "bash analysis/modules/download_encode/methylation_threeFilesCombine.sh liver ENCFF863DRW ENCFF244CRP ENCFF505YJQ"

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

rule skin_ctcf_downloadWrangle:
    output: "data/skin/track_data/ctcf/ctcf.bed.gz"
    conda: "conda_snakeSomMut_env.yml"
    shell: "bash analysis/modules/download_encode/test.sh skin ctcf ENCFF848HOS.bigWig https://www.encodeproject.org/files/ENCFF848HOS/@@download/ENCFF848HOS.bigWig"

rule skin_methylation_downloadWrangle:
    output: "data/skin/track_data/methylation/methylation.bed.gz"
    conda: "conda_snakeSomMut_env.yml"
    shell: "bash analysis/modules/download_encode/methylation_threeFilesCombine.sh skin ENCFF019FIO ENCFF221CNC ENCFF706NQD"


#GERMLINE TRACK RULES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        

rule germline_ctcf_downloadWrangle:
    output: "data/germline/track_data/ctcf/ctcf.bed.gz"
    conda: "conda_snakeSomMut_env.yml"
    shell: "bash analysis/modules/download_encode/test.sh germline ctcf ENCFF883KMQ.bigWig https://www.encodeproject.org/files/ENCFF883KMQ/@@download/ENCFF883KMQ.bigWig"

rule germline_methylation_downloadWrangle:
    output: "data/germline/track_data/methylation/methylation.bed.gz"
    conda: "conda_snakeSomMut_env.yml"
    shell: "bash analysis/modules/download_encode/methylation_threeFilesCombine.sh germline ENCFF602MJL ENCFF082IED  ENCFF887HIO"

rule germline_H3k27ac_download:
    output: "data/germline/track_data/H3k27/H3k27ac_male_hg18_sorted.bed.gz"
    shell: "bash analysis/germline/track_data/H3k27/germline_H3k27ac_downloadWrangle.sh" 

rule germline_dnase_download:
    output: "data/germline/track_data/DNAse/DNAse_male_hg18_sorted.bed.gz"
    shell: "bash analysis/germline/track_data/DNAse/germline_dnase_downloadWrangle.sh" 

rule germline_txn_download:
    output: "data/germline/track_data/transcription/transcription_male_hg18_sorted.bed.gz",
    shell: "bash analysis/germline/track_data/transcription/germline_transcription_downloadWrangle.sh" 


#BLOOD TRACK RULES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
rule blood_ctcf_downloadWrangle:
    output: "data/blood/track_data/ctcf/ctcf.bed.gz"
    conda: "conda_snakeSomMut_env.yml"
    shell: "bash analysis/modules/download_encode/test.sh blood ctcf !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" 

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
    output: "data/blood/track_data/methylation/methylation.bed.gz"
    conda: "conda_snakeSomMut_env.yml"
    shell: "bash analysis/modules/download_encode/methylation_threeFilesCombine.sh blood ENCFF187PPX ENCFF348PBG ENCFF233CYB"

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

rule annotation_downloadWrangle:
    output: "data/global/track_data/annotation/annotation.bed.gz"
    conda: "conda_snakeSomMut_env.yml"
    shell: "bash analysis/global/track_data/annotation/global_annotation_downloadWrangle.sh"

rule mappability_downloadWrangle:
    output: "data/global/track_data/mappability/mappability.bed.gz"
    conda: "conda_snakeSomMut_env.yml"
    shell: "bash analysis/global/track_data/mappability/global_mappability_downloadWrangle.sh"

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


#SCRATCH 
rule createModel_rmCpG: 
    input: "data/{tissue}/dataframes/{model}/predictorDf.txt"
    output:"data/{tissue}/dataframes/{model}/{tissue}_coefDF_boot10k_noCpG.csv",
           "data/{tissue}/dataframes/{model}/{tissue}_all_data_readyForPrediction_boot10k_noCpG.csv"
    conda: "conda_RcreateDfModel_env.yml"
    shell: "Rscript --vanilla analysis/modules/create_model/create_linearModel.R  {wildcards.tissue} {wildcards.model} 10000 FALSE FALSE TRUE FALSE "#equiv, remove_CpG, remove triplets, remove TCXCCX
    
rule createModel_rmTrip: 
    input: "data/{tissue}/dataframes/{model}/predictorDf.txt"
    output:"data/{tissue}/dataframes/{model}/{tissue}_coefDF_boot10k_noTriplets.csv",
           "data/{tissue}/dataframes/{model}/{tissue}_all_data_readyForPrediction_boot10k_noTriplets.csv"
    conda: "conda_RcreateDfModel_env.yml"
    shell: "Rscript --vanilla analysis/modules/create_model/create_linearModel.R  {wildcards.tissue} {wildcards.model} 10000 FALSE TRUE FALSE FALSE "#equiv, remove_CpG, remove triplets, remove TCXCCX

rule createModel_bloodEquiv: 
    input: "data/{tissue}/dataframes/{model}/predictorDf.txt"
    output:"data/{tissue}/dataframes/{model}/{tissue}_coefDF_bloodEquiv_boot10k.csv",
           "data/{tissue}/dataframes/{model}/{tissue}_all_data_readyForPrediction_bloodEquiv_boot10k.csv"
    conda: "conda_RcreateDfModel_env.yml"
    shell: "Rscript --vanilla analysis/modules/create_model/create_linearModel.R  {wildcards.tissue} {wildcards.model} 10000 TRUE FALSE FALSE FALSE"#equiv, remove_CpG, remove triplets, remove TCXCCX

rule createModel_bloodEquiv_rmCpG: 
    input: "data/{tissue}/dataframes/{model}/predictorDf.txt"
    output:"data/{tissue}/dataframes/{model}/{tissue}_coefDF_bloodEquiv_boot10k_noCpG.csv",
           "data/{tissue}/dataframes/{model}/{tissue}_all_data_readyForPrediction_bloodEquiv_boot10k_noCpG.csv"
    conda: "conda_RcreateDfModel_env.yml"
    shell: "Rscript --vanilla analysis/modules/create_model/create_linearModel.R {wildcards.tissue} {wildcards.model} 10000 TRUE TRUE FALSE FALSE"#equiv, remove_CpG, remove triplets, remove TCXCCX

rule createModel_bloodEquiv_rmTCXCCX: 
    input: "data/{tissue}/dataframes/{model}/predictorDf.txt"
    output:"data/{tissue}/dataframes/{model}/{tissue}_coefDF_bloodEquiv_boot10k_noTCXCCX.csv",
           "data/{tissue}/dataframes/{model}/{tissue}_all_data_readyForPrediction_bloodEquiv_boot10k_noTCXCCX.csv"
    conda: "conda_RcreateDfModel_env.yml"
    shell: "Rscript --vanilla analysis/modules/create_model/create_linearModel.R {wildcards.tissue} {wildcards.model} 10000 TRUE FALSE FALSE TRUE"#equiv, remove_CpG, remove triplets, remove TCXCCX
    
rule createModel_bloodEquiv_rmCpG_rmTCXCCX: 
    input: "data/{tissue}/dataframes/{model}/predictorDf.txt"
    output:"data/{tissue}/dataframes/{model}/{tissue}_coefDF_bloodEquiv_boot10k_noCpG_noTCXCCX.csv",
           "data/{tissue}/dataframes/{model}/{tissue}_all_data_readyForPrediction_bloodEquiv_boot10k_noCpG_noTCXCCX.csv"
    conda: "conda_RcreateDfModel_env.yml"
    shell: "Rscript --vanilla analysis/modules/create_model/create_linearModel.R {wildcards.tissue} {wildcards.model} 10000 TRUE TRUE FALSE TRUE"#equiv, remove_CpG, remove triplets, remove TCXCCX
    
rule createModel_bloodEquiv_rmTrip: 
    input: "data/{tissue}/dataframes/{model}/predictorDf.txt"
    output:"data/{tissue}/dataframes/{model}/{tissue}_coefDF_bloodEquiv_boot10k_noTriplets.csv",
           "data/{tissue}/dataframes/{model}/{tissue}_all_data_readyForPrediction_bloodEquiv_boot10k_noTriplets.csv"
    conda: "conda_RcreateDfModel_env.yml"
    shell: "Rscript --vanilla analysis/modules/create_model/create_linearModel.R {wildcards.tissue} {wildcards.model} 10000 TRUE FALSE TRUE FALSE"#equiv, remove_CpG, remove triplets, remove TCXCCX