library(dplyr)

args = commandArgs(trailingOnly=TRUE)
tissue = args[1]
# tissue = "liver"
tissue_predOn = args[2]
# tissue_predOn = "liver"
model_name = args[3]
# model_name = "model4"
tmp_file_path = ""

if ((tissue_predOn == "liver")&&(tissue != "liver")){
    load(paste(tmp_file_path,"data/",tissue,"/objects/",model_name,"/",tissue,"_forLiver_model.RData",sep=""))#model
}else{
    load(paste(tmp_file_path,"data/",tissue,"/objects/",model_name,"/",tissue,"_model.RData",sep=""))#model
}
if ((tissue_predOn != "liver")&&(tissue == "liver")){
    load(paste(tmp_file_path,"data/",tissue_predOn,"/objects/",model_name,"/",tissue_predOn,"_forLiver_samples_sites_test.RData",sep=""))#sample_sites_test
    all_data <- read.csv(paste(tmp_file_path,"data/",tissue_predOn,"/dataframes/",model_name,"/",tissue_predOn,"_forLiver_all_data_readyForPrediction.csv",sep=""),header=TRUE)
}else{
    load(paste(tmp_file_path,"data/",tissue_predOn,"/objects/",model_name,"/",tissue_predOn,"_samples_sites_test.RData",sep=""))#sample_sites_test
    all_data <- read.csv(paste(tmp_file_path,"data/",tissue_predOn,"/dataframes/",model_name,"/",tissue_predOn,"_all_data_readyForPrediction.csv",sep=""),header=TRUE)
}
all_data$annotation <- gsub("protein_binding","not_transcribed",all_data$annotation)#there are so few protein binding sites that we may as well omit 
#rerplacing with "not_transcribed" is fine as they are the last level before non_transcribed in the annotation module     
probs <- predict.glm(model, all_data[sample_sites_test,], type="response")
probs_df <- data.frame(x = probs)
#binding the predicted proabilities to the OG data (but the test sites only)
probs_df <- probs_df %>%
    bind_cols(all_data[sample_sites_test,])
#renaming the glm_probs coloumn
colnames(probs_df ) <- replace(colnames(probs_df ), 1, "glm_probs")
#writing to file 
filename = paste(tmp_file_path,"data/",tissue,"/dataframes/",model_name,"/",tissue,"_on_",tissue_predOn,"_ProbabilityDf.csv",sep="")#this sep is for the filename string
write.csv(probs_df ,filename, row.names = FALSE)


