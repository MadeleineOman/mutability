library(dplyr)
library(glmnet)
library(stringr)
library(stringi)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
# tissue = args[1]
# model_name = args[2]
tissue = "liver"
# tissue_predOn = "liver"
model_name = "model10"
tmp_file_path = "../../../"

all_data <- read.csv("../../../data/liver/dataframes/model10/liver_all_data_readyForPrediction_equiv_toLowest.csv")

# n_muts = nrow(filter(all_data,all_data$mutation_status==1))
# only_muts = all_data %>% filter(mutation_status == 1)
# only_nonMuts = all_data %>% filter(mutation_status == 0)
# test = only_nonMuts[sample(nrow(all_data), nrow(all_data)),]

#create the predictor and reponse for the model input. any NA OMit? can we use the same indexies? 
predictor_matrix = model.matrix(mutation_status~., all_data)[,-1]   #got this from the book, idk . removes NA coloumn and (takes out) the reponse coloum cool cool cool 
response = all_data[['mutation_status']]
stopifnot(nrow(all_data)==nrow(predictor_matrix))


#subsetting into training and testing 
#training 
sample_sites_train <- sample (1:nrow(predictor_matrix), nrow(predictor_matrix)/2) 
matrix_train = predictor_matrix[sample_sites_train,]
response_train = response[sample_sites_train]
#testing
sample_sites_test = c(1:nrow(predictor_matrix))[-sample_sites_train]
matrix_test = predictor_matrix[sample_sites_test,]
response_test = response[sample_sites_test]

#making the model 
model <- cv.glmnet(matrix_train, response_train,alpha =0, family="binomial", nfolds = 10, type.measure = "mse") 

#plotting the cv 
# filename = paste(tmp_file_path,"analysis/",tissue,"/plots/",model_name,"/",tissue,"cvPlot.jpeg", sep="")
# jpeg(filename)
# plot(model,main=tissue)
# dev.off() 

#printing the lambda values
# string_to_print = paste("minimum lambda is ",model$lambda.min,sep=" ")
# cat(string_to_print,file=error_output_file,sep="\n",append=TRUE)
# string_to_print = paste("1se lambda is ",model$lambda.1se,sep=" ")
# cat(string_to_print,file=error_output_file,sep="\n",append=TRUE)

jpeg((paste(tmp_file_path,"analysis/global/plots/model10/liver_regularization_CV.jpeg", sep="")))
plot(model,main=tissue)
dev.off()


#saving the model variables 
#model
filename = paste(tmp_file_path,"data/",tissue,"/objects/",model_name,"/",tissue,"_model_lasso.RData", sep="")
save(model, file=filename)
#testing sites 
filename =paste(tmp_file_path,"data/",tissue,"/objects/",model_name,"/",tissue,"_samples_sites_test_lasso.RData", sep="")
save(sample_sites_test , file=filename)
#testing matrix 
filename = paste(tmp_file_path,"data/",tissue,"/objects/",model_name,"/",tissue,"_matrix_test_lasso.RData", sep="")
save(matrix_test , file=filename)

#analyze the coeficients~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
coefs <- coef(model,s="lambda.min")

coef_df <- data.frame(name = coefs@Dimnames[[1]][coefs@i + 1], coefficient = coefs@x)
#getting the coef output (weird matrix format) into a daatframe : https://stackoverflow.com/questions/27801130/extracting-coefficient-variable-names-from-glmnet-into-a-data-frame

#this si the model without regularization 
normalModel_coef <- read.csv("../../../data/liver/dataframes/model10/liver_coefDF_equiv_toLowest.csv")

comb_data <- merge(normalModel_coef,coef_df)
comb_data <- comb_data %>% filter(name!="(Intercept)")


ggplot(comb_data,aes(value,coefficient))+
    geom_point()+
    theme_light()+
    geom_smooth(method='lm', formula= y~x)+
    xlab("Beta coeficient from model without regularization")+
    ylab("Beta coeficient from model with regularization")+
    theme(
    axis.text = element_text(size = 20, family = 'Helvetica', color = 'black'),
    axis.title = element_text(size = 20, family = 'Helvetica')
    ) 
ggsave(paste(tmp_file_path,"analysis/global/plots/model10/liver_with_without_regularization.jpeg",sep=""))

#now get the R value 
summary(lm(comb_data$coefficient~ comb_data$value))