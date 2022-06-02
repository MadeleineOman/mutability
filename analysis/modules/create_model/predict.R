library(glmnet)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
tissue = args[1]
tissue_predOn = args[2]
model_name = args[3]

load(paste("data/",tissue,"/objects/",model_name,"/",tissue,"_model.RData",sep=""))#model

load(paste("data/",tissue_predOn,"/objects/",model_name,"/",tissue_predOn,"_matrix_test.RData",sep=""))#matrix_test
load(paste("data/",tissue_predOn,"/objects/",model_name,"/",tissue_predOn,"_samples_sites_test.RData",sep=""))#sample_sites_test
all_data <- read.table(paste("data/",tissue_predOn,"/dataframes/",model_name,"/predictorDf.txt",sep=""), sep="\t",header=TRUE)

#doing the prediction 
probs = predict(model$glmnet.fit, newx = matrix_test, s = model$lambda.min, type = "response")
probs_df <- data.frame(x = probs)

#binding the predicted proabilities to the OG data (but the test sites only)
probs_df <- probs_df %>%
    bind_cols(all_data[sample_sites_test,])

#renaming the glm_probs coloumn
colnames(probs_df ) <- replace(colnames(probs_df ), 1, "glm_probs")

#writing to file 
filename = paste("data/",tissue,"/dataframes/",model_name,"/",tissue,"_on_",tissue_predOn,"_ProbabilityDf.csv",sep="")#this sep is for the filename string
write.csv(probs_df ,filename, row.names = FALSE)


