library(dplyr)
library(stringr)
library(stringi)
library(corrplot)
library(rlang)

args = commandArgs(trailingOnly=TRUE)
tissue = args[1]
model_name = args[2]
# tissue = "germline"
# model_name = "model3"
tmp_file_path = ""

error_output_file = paste(tmp_file_path,"data/",tissue,"/objects/",model_name,"/",tissue,"_create_model_text_output.txt",sep="")
cat("output for the create_model notebook",file=error_output_file,sep="\n")

#import the data 
input_filePath = paste(tmp_file_path,"data/",tissue,"/dataframes/",model_name,"/predictorDf.txt",sep="")
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

#editing the annotation colun --> ignored into non_transcribed 
all_data$annotation <- gsub("ignored","not_transcribed",all_data$annotation)
all_data$annotation <- gsub("protein_binding","not_transcribed",all_data$annotation)#there are so few protein binding sites that we may as well omit 
#rerplacing with "not_transcribed" is fine as they are the last level before non_transcribed in the annotation module 

#factorizing the columns 
all_data$mutation_status <- as.factor(all_data$mutation_status)
all_data$triplet <- as.factor(all_data$triplet)
all_data$Chromosome <- as.factor(all_data$Chromosome)
all_data$annotation <- as.factor(all_data$annotation)

#checking that the factor variables are correct (mutation status, triplets, chroms) --> will raise error if not true 
# stopifnot(length(levels(all_data$mutation_status)) == 2)
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

#printing the categorical mut ratios
cat("chromosome mutation distributions",file=error_output_file,sep="\n",append=TRUE)
test_df<-as.data.frame(table(all_data$mutation_status, all_data$Chromosome))
colnames(test_df) <- c("mut_stat","Chromosome","Freq")
for (cur_annot in levels(all_data$Chromosome)) {
    string_to_print = (paste(cur_annot,round(filter(test_df, (Chromosome==cur_annot)&(mut_stat ==1))$Freq / sum(filter(test_df, (Chromosome==cur_annot))$Freq),digits=3),sep="  \t"))
    cat(string_to_print,file=error_output_file,sep="\n",append=TRUE)
} 
cat("annotation mutation distributions",file=error_output_file,sep="\n",append=TRUE)
test_df<-as.data.frame(table(all_data$mutation_status, all_data$annotation))
colnames(test_df) <- c("mut_stat","annotation","Freq")
for (cur_annot in levels(all_data$annotation)) {
    string_to_print = (paste(cur_annot,round(filter(test_df, (annotation==cur_annot)&(mut_stat ==1))$Freq / sum(filter(test_df, (annotation==cur_annot))$Freq),digits=3),sep="  \t"))
    cat(string_to_print,file=error_output_file,sep="\n",append=TRUE)
} 
cat("triplet mutation distributions",file=error_output_file,sep="\n",append=TRUE)
test_df<-as.data.frame(table(all_data$mutation_status, all_data$triplet))
colnames(test_df) <- c("mut_stat","triplet","Freq")
for (cur_annot in levels(all_data$triplet)) {
    string_to_print = (paste(cur_annot,round(filter(test_df, (triplet==cur_annot)&(mut_stat ==1))$Freq / sum(filter(test_df, (triplet==cur_annot))$Freq),digits=3),sep="  \t"))
    cat(string_to_print,file=error_output_file,sep="\n",append=TRUE)
} 


#distributions ~~~~~~~~~~~~~~~~~~~~~
#getting only numeric 
num_cols <- colnames(all_data[,unlist(lapply(all_data, is.numeric), use.names = FALSE)])#https://stackoverflow.com/questions/5863097/selecting-only-numeric-columns-from-a-data-frame

pdf(paste(paste(tmp_file_path,"analysis/",tissue,"/plots/",model_name,"/",tissue,"_predictor_distributions.pdf", sep="")))
for (cur_pred in num_cols){
    cur_data = all_data[,cur_pred]
    title = paste(cur_pred,"\nsd:",round(sd(cur_data),digits=2), " mean:", round(mean(cur_data),digits=2), " median:",round(median(cur_data),digits=2)," min,max:",round(min(cur_data),digits=5),",",max(round(cur_data,digits=2)))
    hist(cur_data, main=title,breaks=20)
}
dev.off()
pdf(paste(paste(tmp_file_path,"analysis/",tissue,"/plots/",model_name,"/",tissue,"_predictor_distributions_remove0.pdf", sep="")))
for (cur_pred in num_cols){
    cur_data = filter(all_data,!!sym(cur_pred)>0)[,cur_pred]# get r to handle strign as value https://stackoverflow.com/questions/48219732/pass-a-string-as-variable-name-in-dplyrfilter
    title = paste(cur_pred,"\nsd:",round(sd(cur_data),digits=2), " mean:", round(mean(cur_data),digits=2), " median:",round(median(cur_data),digits=2)," min,max:",round(min(cur_data),digits=5),",",max(round(cur_data,digits=2)))
    hist(cur_data, main=title,breaks=20)
}
dev.off()

#STANDARDIZING~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (tissue=="liver") {
all_data[4:33] <- as.data.frame(scale(all_data[4:33], center=TRUE, scale=TRUE))
all_data[35:37] <- as.data.frame(scale(all_data[35:37], center=TRUE, scale=TRUE))
} else {
all_data[4:36] <- as.data.frame(scale(all_data[4:36], center=TRUE, scale=TRUE))
all_data[38:40] <- as.data.frame(scale(all_data[38:40], center=TRUE, scale=TRUE))
}

sample_sites_train <- sample (1:nrow(all_data), nrow(all_data)/2) 
sample_sites_test = c(1:nrow(all_data))[-sample_sites_train]

#creating and saving the corplot 
num_data <- na.omit(subset(all_data, select=-c(Chromosome,triplet,annotation)))
num_data$mutation_status <- as.numeric(num_data$mutation_status)
num_cor <- cor(num_data)
filename = paste(tmp_file_path,"analysis/",tissue,"/plots/",model_name,"/",tissue,"_corPlot.pdf", sep="")
pdf(filename)
corrplot(num_cor,method = "color",type="upper",number.cex=0.5,tl.cex = 0.5,main=tissue)
dev.off() 

#CORRELATION ANALYSIS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#create the correlation df and 
cor_df <- as.data.frame(which(num_cor >=0.8,arr.ind=T))
cor_df$col_name <- colnames(num_cor)[cor_df$col]  #insert tyhe col names 
testFunc <- function(a, b) round(num_cor[a,b],digits=3)#https://stackoverflow.com/questions/15059076/call-apply-like-function-on-each-row-of-dataframe-with-multiple-arguments-from-e
cor_df$cor_val <- apply(cor_df[,c('row','col')], 1, function(y) testFunc(y['row'],y['col'])) #insert the value of the correlatnio 
#exclude the self corelations and convert row names into a col 
cor_df<- filter(as.data.frame(cor_df),row!=col)[,c("col_name","cor_val")]
cor_df<- tibble::rownames_to_column(cor_df, "row_name")
#reorder the var names into alphabetical order 
testFunc <- function(a,b) paste(str_sort(c(a,b)),sep="-")[1]
cor_df$first_name <- apply(cor_df,1,function(y) testFunc(y['row_name'],y['col_name']))
testFunc <- function(a,b) paste(str_sort(c(a,b)),sep="-")[2]
cor_df$last_name <- apply(cor_df,1,function(y) testFunc(y['row_name'],y['col_name']))
#subset cols
cor_df <- cor_df[,c("first_name","last_name","cor_val")]
#doing soem strign replace to maker the distinct below applicable 
cor_df <- data.frame(lapply(cor_df, function(x) {gsub("1[.]1", "1", x)}))
cor_df <- data.frame(lapply(cor_df, function(x) {gsub("[.]2", "", x)}))
cor_df <- data.frame(lapply(cor_df, function(x) {gsub("0[.]1", "0", x)}))
#remove the duplicates and reorder for legibility 
cor_df <- distinct(cor_df)
cor_df <- cor_df[order(cor_df$first_name),]
#save to file 
filename = paste(tmp_file_path,"data/",tissue,"/dataframes/",model_name,"/",tissue,"_fullModel_corTable_above0.8.csv", sep="")
write.csv(cor_df,filename, row.names = FALSE)
                          
#REMOVING CORRELATED VARIABLES 
if (tissue =="liver"){ #liver need its own b/c no H3k36me3/H3k27me3
    all_data <- all_data[,c("Chromosome","triplet","mutation_status",
                            'H3k27.10000','H3k4me1.10000','H3k4me3.10000',"Transcription.10000","Replication.10000","recombination.10000","DNAse.100",
                            'laminB1.1','laminB1.100','laminB1.10000', "Repeats.1","Repeats.100","Repeats.10000","GC_content.1","GC_content.100","GC_content.10000","annotation")]
}else{
        all_data <- all_data[,c("Chromosome","triplet","mutation_status",
                            'H3k27.10000','H3k4me1.10000','H3k4me3.10000','H3k27me3.10000','H3k36me3.10000',"Transcription.10000","Replication.10000","recombination.10000","DNAse.100",
                            'laminB1.1','laminB1.100','laminB1.10000', "Repeats.1","Repeats.100","Repeats.10000","GC_content.1","GC_content.100","GC_content.10000","annotation")]
}

#MAKING A MODEL ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
model <- glm(mutation_status~., data=all_data[sample_sites_train,],family="binomial")

#saving the model variables 
filename = paste(tmp_file_path,"data/",tissue,"/objects/",model_name,"/",tissue,"_model.RData", sep="") #model 
save(model, file=filename)
filename =paste(tmp_file_path,"data/",tissue,"/objects/",model_name,"/",tissue,"_samples_sites_test.RData", sep="")#sample sites test 
save(sample_sites_test , file=filename)
filename =paste(tmp_file_path,"data/",tissue,"/dataframes/",model_name,"/",tissue,"_all_data_readyForPrediction.csv", sep="")
write.csv(all_data,filename)
                          
#checing the assumptions.. though idk how to read 
pdf(paste(tmp_file_path,"analysis/",tissue,"/plots/",model_name,"/",tissue,"_model_diagnostics.pdf", sep=""))
plot(model)
dev.off() 
#checking for linearity one at a time... 
pdf(paste(tmp_file_path,"analysis/",tissue,"/plots/",model_name,"/",tissue,"_model_diagnostics_checkingLinearity.pdf", sep=""))
for (pred_name in colnames(all_data)){
    plot(model$linear.predictors ~ all_data[sample_sites_train,][,pred_name],main=pred_name)
}
dev.off()

#analyze the coeficients~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#first the whole model 
coefs <- coef(summary(model))
coef_df <- as.data.frame(coefs)
colnames(coef_df) <- c("value","std_err","z_val","p_val")
coef_df<- tibble::rownames_to_column(coef_df, "name") # https://stackoverflow.com/questions/29511215/convert-row-names-into-first-column
coef_df_ordered <- coef_df[order(-coef_df$value),]#https://www.statmethods.net/management/sorting.html
filename = paste(tmp_file_path,"data/",tissue,"/dataframes/",model_name,"/",tissue,"_coefDF.csv",sep="")#this sep is for the filename string
write.csv(coef_df_ordered,filename,row.names=FALSE)

#BOOTSTRAP                         
#my bootstrapping help https://stackoverflow.com/questions/54749641/bootstrapping-with-glm-model
coef_df<- as.data.frame("names" <- names(coef(model)))
colnames(coef_df) <- c("name")

print("bootstrap started")
plm <- proc.time()
for (i in 1:100){
    sample_sites_train_cur <- sample(sample_sites_train,size=length(sample_sites_train),replace=TRUE)
    cur_data <- all_data[sample_sites_train_cur,]
    model<- glm(mutation_status~., data=cur_data,family="binomial")

    cur_coef_df <- as.data.frame(coef(model))
    colnames(cur_coef_df) <- c(paste("value_",i,sep=""))
    cur_coef_df<- tibble::rownames_to_column(cur_coef_df, "name") # https://stackoverflow.com/questions/29511215/convert-row-names-into-first-column
    coef_df <- merge(coef_df,cur_coef_df)
}
paste("bootstrap took",proc.time()-plm,"seconds",sep=" ")

#getting the mean and quantiles into a df for plotting 
coef_df$mean_est <- rowMeans(coef_df[,2:101],na.rm=TRUE)
quantile_df <- data.frame(coef_df$name,t(apply(coef_df[,2:101],1,quantile,c(0.025,0.975),na.rm=TRUE)))#https://stackoverflow.com/questions/54749641/bootstrapping-with-glm-model
colnames(quantile_df) <- c("name","quant025","quant975")
bootstrap_df<- merge(coef_df,quantile_df)
filename <- paste(tmp_file_path,"data/",tissue,"/dataframes/",model_name,"/",tissue,"_coefDF_bootstrap.csv",sep="")
write.csv(bootstrap_df,filename,row.names=FALSE)

#distribtuion analysis 
sign_preds = replicate(nrow(bootstrap_df),"not")
sign_preds[bootstrap_df$quant025 * bootstrap_df$quant975 >0 ] = "yes"
bootstrap_df$significant <- sign_preds
num_cols <- colnames(all_data[,unlist(lapply(all_data, is.numeric), use.names = FALSE)])#https://stackoverflow.com/questions/5863097/selecting-only-numeric-columns-from-a-data-frame
#making the pllots 
pdf(paste(paste(tmp_file_path,"analysis/",tissue,"/plots/",model_name,"/",tissue,"_predictor_distributions_unCor_includedModel.pdf", sep="")))
for (cur_pred in num_cols){
    cur_data = all_data[,cur_pred]
    title = paste(cur_pred,"\nsd:",round(sd(cur_data),digits=2), " mean:", round(mean(cur_data),digits=2), " median:",round(median(cur_data),digits=2),
                  " min,max:",round(min(cur_data),digits=5),",",max(round(cur_data,digits=2)),
                  "\nboostrap mean est:",round(filter(bootstrap_df, name==cur_pred)$mean_est,digits=5), "singificant?:", filter(bootstrap_df, name==cur_pred)$significant)
    hist(cur_data, main=title,breaks=20)
}
dev.off()
                          
                          
#CREATING MODELS FOR THE LIVER TISSUE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (tissue != "liver"){
    all_data <- all_data[,c("Chromosome","triplet","mutation_status",
                            'H3k27.10000','H3k4me1.10000','H3k4me3.10000','H3k27me3.10000','H3k36me3.10000',"Transcription.10000","Replication.10000","recombination.10000","DNAse.100",
                            'laminB1.1','laminB1.100','laminB1.10000', "Repeats.1","Repeats.100","Repeats.10000","GC_content.1","GC_content.100","GC_content.10000","annotation")]

    #one whole model for 
    model <- glm(mutation_status~., data=all_data[sample_sites_train,],family="binomial")

    #saving the model variables 
    filename = paste(tmp_file_path,"data/",tissue,"/objects/",model_name,"/",tissue,"_forLiver_model.RData", sep="") #model 
    save(model, file=filename)
    filename =paste(tmp_file_path,"data/",tissue,"/objects/",model_name,"/",tissue,"_forLiver_samples_sites_test.RData", sep="")#sample sites test 
    save(sample_sites_test , file=filename)
    filename =paste(tmp_file_path,"data/",tissue,"/dataframes/",model_name,"/",tissue,"_forLiver_all_data_readyForPrediction.csv", sep="")
    write.csv(all_data,filename)

    #analyze the coeficients~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    coefs <- coef(summary(model))
    coef_df <- as.data.frame(coefs)
    colnames(coef_df) <- c("value","std_err","z_val","p_val")
    coef_df<- tibble::rownames_to_column(coef_df, "name") # https://stackoverflow.com/questions/29511215/convert-row-names-into-first-column
    coef_df_ordered <- coef_df[order(-coef_df$value),]#https://www.statmethods.net/management/sorting.html
    filename = paste(tmp_file_path,"data/",tissue,"/dataframes/",model_name,"/",tissue,"_forLiver_coefDF.csv",sep="")#this sep is for the filename string
    write.csv(coef_df_ordered,filename,row.names=FALSE)

    #MODEL DIAGNOSTICS 
    #checing the assumptions.. though idk how to read 
    pdf(paste(tmp_file_path,"analysis/",tissue,"/plots/",model_name,"/",tissue,"_forLiver_model_diagnostics.pdf", sep=""))
    plot(model)
    dev.off() 
    #checking for linearity one at a time... 
    pdf(paste(tmp_file_path,"analysis/",tissue,"/plots/",model_name,"/",tissue,"_forLiver_model_diagnostics_checkingLinearity.pdf", sep=""))
    for (pred_name in colnames(all_data)){
        plot(model$linear.predictors ~ all_data[sample_sites_train,][,pred_name],main=pred_name)
    }
    dev.off()

    #BOOTSTRAP 
    #my bootstrapping help https://stackoverflow.com/questions/54749641/bootstrapping-with-glm-model
    coef_df<- as.data.frame("names" <- names(coef(model)))
    colnames(coef_df) <- c("name")
    print("started forLiver bootstrap")
    plm <- proc.time()
    for (i in 1:100){
        sample_sites_train_cur <- sample(sample_sites_train,size=length(sample_sites_train),replace=TRUE)
        cur_data <- all_data[sample_sites_train_cur,]
        model<- glm(mutation_status~., data=cur_data,family="binomial")

        cur_coef_df <- as.data.frame(coef(model))
        colnames(cur_coef_df) <- c(paste("value_",i,sep=""))
        cur_coef_df<- tibble::rownames_to_column(cur_coef_df, "name") # https://stackoverflow.com/questions/29511215/convert-row-names-into-first-column
        coef_df <- merge(coef_df,cur_coef_df)
    }
    paste("forLiver bootstrap took",proc.time()-plm,"seconds",sep=" ")

    #getting the mean and quantiles into a df for plotting 
    coef_df$mean_est <- rowMeans(coef_df[,2:101],na.rm=TRUE)
    quantile_df <- data.frame(coef_df$name,t(apply(coef_df[,2:101],1,quantile,c(0.025,0.975),na.rm=TRUE)))#https://stackoverflow.com/questions/54749641/bootstrapping-with-glm-model
    colnames(quantile_df) <- c("name","quant025","quant975")
    bootstrap_df<- merge(coef_df,quantile_df)
    filename <- paste(tmp_file_path,"data/",tissue,"/dataframes/",model_name,"/",tissue,"_forLiver_coefDF_bootstrap.csv",sep="")
    write.csv(bootstrap_df,filename,row.names=FALSE)

    #creating a significant column 
    sign_preds = replicate(nrow(bootstrap_df),"not")
    sign_preds[bootstrap_df$quant025 * bootstrap_df$quant975 >0 ] = "yes"
    bootstrap_df$significant <- sign_preds

    num_cols <- colnames(all_data[,unlist(lapply(all_data, is.numeric), use.names = FALSE)])#https://stackoverflow.com/questions/5863097/selecting-only-numeric-columns-from-a-data-frame
    #making the pllots 
    pdf(paste(paste(tmp_file_path,"analysis/",tissue,"/plots/",model_name,"/",tissue,"_forLiver_predictor_distributions_unCor_includedModel.pdf", sep="")))
    for (cur_pred in num_cols){
        cur_data = all_data[,cur_pred]
        title = paste(cur_pred,"\nsd:",round(sd(cur_data),digits=2), " mean:", round(mean(cur_data),digits=2), " median:",round(median(cur_data),digits=2),
                      " min,max:",round(min(cur_data),digits=5),",",max(round(cur_data,digits=2)),
                      "\nboostrap mean est:",round(filter(bootstrap_df, name==cur_pred)$mean_est,digits=5), "singificant?:", filter(bootstrap_df, name==cur_pred)$significant)
        hist(cur_data, main=title,breaks=20)
    }
    dev.off()
    #removing the 0s doesnt make any sense for standardized predictors 
}
