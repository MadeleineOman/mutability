library(dplyr)


args = commandArgs(trailingOnly=TRUE)
tissue = args[1]
tissue_predOn = args[2]
model_name = args[3]
if ("_equiv_toLowest" %in% args){ equiv_toLowest = TRUE
}else{equiv_toLowest = FALSE }
if ("_noTriplets" %in% args){exclude_triplet = TRUE
}else{ exclude_triplet = FALSE }
if ("_noCpG" %in% args){exclude_CpG = TRUE
}else{exclude_CpG = FALSE }
if ("_noTCX_CCX" %in% args){ exclude_TCX_CCX = TRUE
}else{ exclude_TCX_CCX = FALSE }


# tissue = "blood"
# tissue_predOn = "skin"
# model_name = "model6"
# equiv_toLowest=TRUE
# exclude_triplet = FALSE
# exclude_CpG= TRUE

tmp_file_path = ""

model_desc_modify = ""

if ("_fullModel" %in% args){model_desc_modify = paste(model_desc_modify,"_fullModel",sep="")}

if (equiv_toLowest==TRUE){
    model_desc_modify = paste(model_desc_modify,"_equiv_toLowest",sep="")}
if (exclude_CpG==TRUE){
    model_desc_modify = paste(model_desc_modify,"_noCpG",sep="")}
if (exclude_triplet==TRUE){
    model_desc_modify = paste(model_desc_modify,"_noTriplets",sep="")}


if ((tissue_predOn == "liver")&&(tissue != "liver")){
    load(paste(tmp_file_path,"data/",tissue,"/objects/",model_name,"/",tissue,"_forLiver_model",model_desc_modify,".RData",sep=""))#model
}else{
    load(paste(tmp_file_path,"data/",tissue,"/objects/",model_name,"/",tissue,"_model",model_desc_modify,".RData",sep=""))#model
}
if ((tissue_predOn != "liver")&&(tissue == "liver")){
    load(paste(tmp_file_path,"data/",tissue_predOn,"/objects/",model_name,"/",tissue_predOn,"_forLiver_samples_sites_test",model_desc_modify,".RData",sep=""))#sample_sites_test
    all_data <- read.csv(paste(tmp_file_path,"data/",tissue_predOn,"/dataframes/",model_name,"/",tissue_predOn,"_forLiver_all_data_readyForPrediction",model_desc_modify,".csv",sep=""),header=TRUE)
}else{
    load(paste(tmp_file_path,"data/",tissue_predOn,"/objects/",model_name,"/",tissue_predOn,"_samples_sites_test",model_desc_modify,".RData",sep=""))#sample_sites_test
    all_data <- read.csv(paste(tmp_file_path,"data/",tissue_predOn,"/dataframes/",model_name,"/",tissue_predOn,"_all_data_readyForPrediction",model_desc_modify,".csv",sep=""),header=TRUE)
}



#rerplacing with "not_transcribed" is fine as they are the last level before non_transcribed in the annotation module     
probs <- predict.glm(model, all_data[sample_sites_test,], type="response")
probs_df <- data.frame(x = probs)
#binding the predicted proabilities to the OG data (but the test sites only)
probs_df <- probs_df %>%
    bind_cols(all_data[sample_sites_test,])
#renaming the glm_probs coloumn
colnames(probs_df ) <- replace(colnames(probs_df ), 1, "glm_probs")
#writing to file 
filename = paste(tmp_file_path,"data/",tissue,"/dataframes/",model_name,"/",tissue,"_on_",tissue_predOn,"_ProbabilityDf",model_desc_modify,".csv",sep="")#this sep is for the filename string
write.csv(probs_df ,filename, row.names = FALSE)
