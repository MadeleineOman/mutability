library(dplyr)
library(glmnet)
library(stringr)
library(stringi)

args = commandArgs(trailingOnly=TRUE)
tissue = args[1]
model_name = args[2]
# tissue = "liver"
# tissue_predOn = "liver"
# model_name = "model2"

error_output_file = paste("data/",tissue,"/objects/",model_name,"/",tissue,"_create_model_text_output.txt",sep="")
cat("output for the create_model notebook",file=error_output_file,sep="\n")

#import the data 
input_filePath = paste("data/",tissue,"/dataframes/",model_name,"/predictorDf.txt",sep="")
all_data <- read.table(input_filePath, header = TRUE,sep="\t")

#this if conditional accounts for the accidental mix up in the smallest scale: 
#sometimes the sclae is 0-bases around (at the site), while sometimes its 1-base around (triplet) 
if (".0" %in% unique(str_extract(colnames(all_data),"[.][0-9]*"))) {
        all_data <- all_data %>%
            mutate(GC_content.0 = Gpercent.0+Cpercent.0) %>%
            mutate(GC_content.100 = Gpercent.100+Cpercent.100) %>% 
            mutate(GC_content.10000 = Gpercent.10000+Cpercent.10000)
        all_data <- all_data[,!(names(all_data) %in% c("site",'Apercent.0','Gpercent.0','Cpercent.0','Tpercent.0',
                                               'Apercent.100','Gpercent.100','Cpercent.100','Tpercent.100',
                                               'Apercent.10000','Gpercent.10000','Cpercent.10000','Tpercent.10000'))]
}else if (".1" %in% unique(str_extract(colnames(all_data),"[.][0-9]*"))){
    all_data <- all_data %>%
        mutate(GC_content.1 = Gpercent.1+Cpercent.1) %>%
        mutate(GC_content.100 = Gpercent.100+Cpercent.100) %>% 
        mutate(GC_content.10000 = Gpercent.10000+Cpercent.10000)
    all_data <- all_data[,!(names(all_data) %in% c("site",'Apercent.1','Gpercent.1','Cpercent.1','Tpercent.1',
                                               'Apercent.100','Gpercent.100','Cpercent.100','Tpercent.100',
                                               'Apercent.10000','Gpercent.10000','Cpercent.10000','Tpercent.10000'))]
}

#editing the triplets 
all_data$triplet <- toupper(all_data$triplet)
#filter for rows that dont have NNN as the triplet --> write the details to file 
cat(paste(nrow(all_data[all_data$triplet == "NNN",]), "rows removed due to N in triplet, ",sep=" "),file=error_output_file,sep="\n",append=TRUE)
all_data <- all_data[all_data$triplet != "NNN",]
cat(paste(nrow(all_data),"rows left",sep=" "),file=error_output_file,sep="\n",append=TRUE)

#how many mutations and sites 
#printing the mutation ratio to file 
string_to_print = paste("mut/normal ratio: ",nrow(all_data[all_data$mutation_status == "1",])/nrow(all_data),sep=" ")
cat(string_to_print,file=error_output_file,sep="\n",append=TRUE)
#printing the total mutations included to file 
string_to_print = paste("n muts = ", nrow(all_data[all_data$mutation_status == "1",]),sep=" ")
cat(string_to_print,file=error_output_file,sep="\n",append=TRUE)
#printing the total rows included to file 
string_to_print = paste("total nrow", nrow(all_data),sep=" ")
cat(string_to_print,file=error_output_file,sep="\n",append=TRUE)

#converting 64-->32 triplets 
rc_removeAG <- function(dna){
    middle_base = substr(dna, 2, 2)
    if(middle_base %in% c("A","G")){
        dna <- stri_reverse(chartr("acgtACGT", "tgcaTGCA", dna))}
    return(dna)
}#substring slicing https://www.johnmyleswhite.com/notebook/2009/02/25/text-processing-in-r/
all_data$triplet <- unlist(lapply(as.character(all_data$triplet),rc_removeAG)) #need unlistto turn the list into a vector 


#check for problems with levels / column tye etc. 
all_data$triplet <- as.character(all_data$triplet) #need to convert to char so the as.factor properly reducs to 64 levels after the removal of 'triplets" longer than 3
string_to_print = paste(nrow(all_data[nchar(as.character(all_data$triplet))!=3,]),"rows removed due to triplet larger than 3 in length ",sep=" ")
cat(string_to_print,file=error_output_file,sep="\n",append=TRUE)
all_data <- all_data[nchar(all_data$triplet)==3,] #make sure only including triplets 
string_to_print=paste(nrow(all_data)," rows left", sept = " ")
cat(string_to_print,file=error_output_file,sep="\n",append=TRUE)

#factorizing the columns 
all_data$mutation_status <- as.factor(all_data$mutation_status)
all_data$triplet <- as.factor(all_data$triplet)
all_data$Chromosome <- as.factor(all_data$Chromosome)

#checking that the factor variables are correct (mutation status, triplets, chroms) --> will raise error if not true 
stopifnot(length(levels(all_data$mutation_status)) == 2)
stopifnot(length(levels(all_data$triplet))==32)
stopifnot(length(levels(all_data$Chromosome))==22)

#if germline, then edit out the female and make the male the default ( i used to include both female and male) 
if (tissue == "germline"){
    all_data <- select(all_data,-matches("female")) 
    colnames(all_data)<- str_replace_all(colnames(all_data), "_male", "")
}

#na omit and checking how many rows were removed 
string_to_print = paste(nrow(all_data)- nrow(na.omit(all_data)), "lost to NA values, ", nrow(na.omit(all_data))," rows remain after",sep=" ")
cat(string_to_print,file=error_output_file,sep="\n",append=TRUE)
all_data <- na.omit(all_data)

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
filename = paste("analysis/",tissue,"/plots/",model_name,"/",tissue,"cvPlot.jpeg", sep="")
jpeg(filename)
plot(model,main=tissue)
dev.off() 

#printing the lambda values
string_to_print = paste("minimum lambda is ",model$lambda.min,sep=" ")
cat(string_to_print,file=error_output_file,sep="\n",append=TRUE)
string_to_print = paste("1se lambda is ",model$lambda.1se,sep=" ")
cat(string_to_print,file=error_output_file,sep="\n",append=TRUE)

#saving the model variables 
#model
filename = paste("data/",tissue,"/objects/",model_name,"/",tissue,"_model.RData", sep="")
save(model, file=filename)
#testing sites 
filename =paste("data/",tissue,"/objects/",model_name,"/",tissue,"_samples_sites_test.RData", sep="")
save(sample_sites_test , file=filename)
#testing matrix 
filename = paste("data/",tissue,"/objects/",model_name,"/",tissue,"_matrix_test.RData", sep="")
save(matrix_test , file=filename)


#analyze the coeficients~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
coefs <- coef(model,s="lambda.min")

coef_df <- data.frame(name = coefs@Dimnames[[1]][coefs@i + 1], coefficient = coefs@x)
#getting the coef output (weird matrix format) into a daatframe : https://stackoverflow.com/questions/27801130/extracting-coefficient-variable-names-from-glmnet-into-a-data-frame

meanColValues <- apply(matrix_train,2,mean)
stderrColValues <- apply(matrix_train,2,sd)/sqrt(nrow(matrix_train))
meanStderr_colValues <-data.frame(meanColValues,stderrColValues)
meanStderr_colValues<- tibble::rownames_to_column(meanStderr_colValues, "name") # https://stackoverflow.com/questions/29511215/convert-row-names-into-first-column

coef_values_df <-merge(coef_df,meanStderr_colValues) #merging the coeff and mean/stderr dataframes 

coef_values_df$bx <-coef_values_df$meanColValues* coef_values_df$coefficient #creating the col that combines the coef (b) witht he values (x)
coef_values_df$bx_stderr <-coef_values_df$stderrColValues* coef_values_df$coefficient

coef_df_ordered <- coef_values_df[order(-coef_values_df$coefficient),]
#https://www.statmethods.net/management/sorting.html

filename = paste("data/",tissue,"/dataframes/",model_name,"/",tissue,"_coefDF.csv",sep="")#this sep is for the filename string
write.csv(coef_df_ordered,filename,row.names=FALSE)