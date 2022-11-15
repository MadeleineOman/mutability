library(dplyr)
library(stringr)
library(stringi)
library(corrplot)
library(rlang)

args = commandArgs(trailingOnly=TRUE)
tissue = args[1]
model_name = args[2]
# tissue = "blood"
# model_name = "model6"
tmp_file_path = ""


error_output_file = paste(tmp_file_path,"data/",tissue,"/objects/",model_name,"/",tissue,"_create_model_text_output.txt",sep="")
cat("output for the create_model notebook",file=error_output_file,sep="\n")

#import the data 
input_filePath = paste(tmp_file_path,"data/",tissue,"/dataframes/",model_name,"/predictorDf.txt",sep="")
all_data <- read.table(input_filePath, header = TRUE,sep="\t")

#chnaging the 1-surround to 0 surround (--> due to a mess up, only blood )
# if (tissue=="blood"){
#     all_data <- all_data %>% 
#   rename(Apercent.0=Apercent.1,Gpercent.0=Gpercent.1,Cpercent.0=Cpercent.1,Tpercent.0=Tpercent.1,
#          methylation_coverage.0=methylation_coverage.1,methylation_precent.0=methylation_precent.1,
#          H3k27.0=H3k27.1,H3k27me3.0=H3k27me3.1,H3k4me1.0=H3k4me1.1,H3k4me3.0=H3k4me3.1,H3k36me3.0=H3k36me3.1,
#          Transcription.0=Transcription.1,recombination.0=recombination.1,Repeats.0=Repeats.1,DNAse.0=DNAse.1,laminB1.0=laminB1.1)
# }

#
all_data$annotation <- gsub("protein_binding","not_transcribed",all_data$annotation)#there are so few protein binding sites that we may as well omit 

#removing unamppable points 
mappable_mut_summary <- all_data %>% 
        group_by(mappability,mutation_status)%>%
    summarise(count=n())
unmapped_muts = (filter(mappable_mut_summary, (mappability == "not")&(mutation_status == 1))$count)
unmapped_nonMuts = (filter(mappable_mut_summary, (mappability == "not")&(mutation_status == 0))$count)
string_to_print = paste(tissue,":",unmapped_muts," mutations and",unmapped_nonMuts,"nonMut sites removed by mappability filter out of ",nrow(all_data),"sites",sep=" ")
cat(string_to_print,file=error_output_file,sep="\n",append=TRUE)
all_data <- filter(all_data,mappability=="mappable")
all_data <- all_data[,!(names(all_data) %in% c("mappability"))]

#gc content 
all_data <- all_data %>%
    mutate(GC_content.0 = Gpercent.0+Cpercent.0) %>%
    mutate(GC_content.100 = Gpercent.100+Cpercent.100) %>% 
    mutate(GC_content.10000 = Gpercent.10000+Cpercent.10000)
all_data <- all_data[,!(names(all_data) %in% c('Apercent.0','Gpercent.0','Cpercent.0','Tpercent.0',
                                       'Apercent.100','Gpercent.100','Cpercent.100','Tpercent.100',
                                       'Apercent.10000','Gpercent.10000','Cpercent.10000','Tpercent.10000'))]


#create df to test coverage in the methylation download notebook 
bases = c("A","T","C","G")
c_trips = c()
for (b1 in bases){
    for (b3 in bases){
        c_trips <- append(c_trips,paste(b1,"C",b3,sep=""))
    }
}
c_trips_df <- filter(all_data, triplet%in%c_trips)[,c("Chromosome","site","triplet",'methylation_precent.0','methylation_coverage.0')]
write.csv(c_trips_df,paste(tmp_file_path,"data/",tissue,"/track_data/methylation/predDf_ctrips.csv",sep=""))

# # plotting METHYLATION coverage  in dataset
# # create the distributions for methylation 
# all_data$test_meth <- as.integer(all_data$methylation_precent.0)
# all_data$test_meth_cov <- as.integer(all_data$methylation_coverage.0)
# to_plot = ((filter(all_data,GC_content.0==1&test_meth_cov <100))$test_meth_cov)
# to_plot_percent = ((filter(all_data,GC_content.0==1&test_meth_cov <100))$test_meth)
# hist(as.integer(to_plot),main=paste(tissue,"coverage",sep=" "),breaks=40)
# hist(as.integer(to_plot_percent),main=paste(tissue,"percent methylation",sep=" "),breaks=40)


#METHYLATION                                 
#create the distributions for methylation 
# all_data$test_meth <- as.integer(all_data$methylation_precent.0)
# all_data$test_meth_cov <- as.integer(all_data$methylation_coverage.0)
# to_plot = ((filter(all_data,GC_content.0==1&test_meth_cov <100))$test_meth_cov)
# hist(as.integer(to_plot),main=paste(tissue,"coverage",sep=" "),breaks=40)
#create the methylabel column 


all_data <- all_data %>% 
      mutate(methylable = ifelse((methylation_precent.0=="no_percent_data" & GC_content.0 == 0), "no",#not a C/G site = no 
               ifelse((methylation_precent.0=="no_percent_data" & GC_content.0 == 1), NA, # C/G site with no data = NA 
               ifelse(methylation_precent.0=="0.0","no","yes")))) #C/G site with 0 = no, >0 = yes 


#https://stackoverflow.com/questions/24459752/can-dplyr-package-be-used-for-conditional-mutating
if(nrow(all_data %>% filter(methylable=="yes" & GC_content.0 == 0)) != 0 ){ #making sure only gc sites are methylated 
        methyl_messed<- filter(all_data,methylable=="yes" & GC_content.0 == 0)[,c("Chromosome","site","triplet","methylation_precent.0","methylable")]
        paste(nrow(methyl_messed),"rows have non-zero methylation values for non-c sites",sep=" ")
        filename = paste(tmp_file_path,"data/",tissue,"/dataframes/",model_name,"/",tissue,"_methyl_messed.csv",sep="")
        write.csv(methyl_messed,filename, row.names = FALSE)
        all_data <- all_data %>% 
          mutate(methylable = ifelse((methylation_precent.0=="infsuff_coverage" & GC_content.0 == 0), "no", methylable))
}
#get methylation summary stats to print to output file 
methyl_mut_summary <- all_data %>% 
        group_by(methylable,mutation_status)%>%
    summarise(count=n())
cantmeth_muts <- (filter(methyl_mut_summary, is.na(methylable)&(mutation_status == 1)))$count
cantmeth_NonMuts <- (filter(methyl_mut_summary, is.na(methylable)&(mutation_status == 0)))$count
string_to_print = paste(tissue,":",cantmeth_muts," mutations and",cantmeth_NonMuts,"nonMut sites removed by no methylable coverage out of ",nrow(all_data),"sites",sep=" ")
cat(string_to_print,file=error_output_file,sep="\n",append=TRUE)
#removing methylation columns 
all_data <- all_data[,!(names(all_data) %in% c("site","methylation_precent.0","methylation_precent.100","methylation_precent.10000",
                                               "methylation_coverage.0","methylation_coverage.100","methylation_coverage.10000"))]
all_data <- all_data[!(is.na(all_data$methylable)),] #getting rid of the methylation nas 


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
all_data$CpGisland <- as.factor(all_data$CpGisland)
all_data$methylable <- as.factor(all_data$methylable)

#checking that the factor variables are correct (mutation status, triplets, chroms) --> will raise error if not true 
stopifnot(length(levels(all_data$mutation_status)) == 2)
stopifnot(length(levels(all_data$triplet))==32)
stopifnot(length(levels(all_data$Chromosome))==22)
stopifnot(length(levels(all_data$CpGisland))==3)
stopifnot(length(levels(all_data$methylable))==2)

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


#making categorical columns into seperate numerical columns ready for standardizing 
all_data$dummy <- 1#need a dumym column that the model.matrix can remove (as it has to remove a column aparently)
muts_col <- all_data$mutation_status#saving the column of mutations so i can add it back later 
non_num_preds = c("mutation_status")
preds_to_standardize<-!(names(all_data) %in% non_num_preds)
all_data <-(data.frame((model.matrix(dummy~., all_data[, preds_to_standardize])[,-1] )))
all_data$mutation_status <- muts_col



#STANDARDIZING~~~~~~~~~~~~~~~~~~~~~~~~~~~~
non_num_preds = c("mutation_status")
preds_to_standardize<-!(names(all_data) %in% non_num_preds)#yes i have to do this again because i chnaged the order of the mutant column 
all_data[, preds_to_standardize]<- (as.data.frame(scale(all_data[, preds_to_standardize],center=TRUE,scale=TRUE)))

sample_sites_train <- sample (1:nrow(all_data), nrow(all_data)/2) 
sample_sites_test = c(1:nrow(all_data))[-sample_sites_train]


#creating and saving the corplot 
# num_data <- na.omit(subset(all_data, select=-c(Chromosome,triplet,annotation,CpGisland))) --> dont need this anymore as correlating ther categoricals too 
num_data <-all_data
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
correlated_vars_toRm = c('H3k27.1','H3k27.0','H3k27.100',
                         'H3k4me1.0','H3k4me1.1','H3k4me1.100',
                         'H3k4me3.0','H3k4me3.1','H3k4me3.100',
                         'H3k27me3.0','H3k27me3.1','H3k27me3.100',
                         'H3k36me3.0','H3k36me3.1','H3k36me3.100',
                         "Transcription.0","Transcription.1","Transcription.100",
                         "recombination.0","recombination.1","recombination.100",
                         "DNAse.0","DNAse.1","DNAse.10000",
                         "GC_content.0",
                         "Repeats.100")
all_data <- all_data[,!(names(all_data) %in% correlated_vars_toRm)]
#redoing rhe corplot after emoving correlations 

num_data <-all_data
num_data$mutation_status <- as.numeric(num_data$mutation_status)
num_cor <- cor(num_data)
filename = paste(tmp_file_path,"analysis/",tissue,"/plots/",model_name,"/",tissue,"_corPlot_pruned.pdf", sep="")
pdf(filename)
corrplot(num_cor,method = "color",type="upper",number.cex=0.5,tl.cex = 0.5,main=tissue)
dev.off() 
                          
                          
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

paste(tissue,"bootstrap started",sep=" ")
plm <- Sys.time()
for (i in 1:100){
    sample_sites_train_cur <- sample(sample_sites_train,size=length(sample_sites_train),replace=TRUE)
    cur_data <- all_data[sample_sites_train_cur,]
    model<- glm(mutation_status~., data=cur_data,family="binomial")

    cur_coef_df <- as.data.frame(coef(model))
    colnames(cur_coef_df) <- c(paste("value_",i,sep=""))
    cur_coef_df<- tibble::rownames_to_column(cur_coef_df, "name") # https://stackoverflow.com/questions/29511215/convert-row-names-into-first-column
    coef_df <- merge(coef_df,cur_coef_df)
}
paste(tissue,"bootstrap took",Sys.time()-plm,"s/mins/hs, you guess",sep=" ")

#getting the mean and quantiles into a df for plotting 
coef_df$mean_est <- rowMeans(coef_df[,2:101],na.rm=TRUE)
quantile_df <- data.frame(coef_df$name,t(apply(coef_df[,2:101],1,quantile,c(0.025,0.975),na.rm=TRUE)))#https://stackoverflow.com/questions/54749641/bootstrapping-with-glm-model
colnames(quantile_df) <- c("name","quant025","quant975")
bootstrap_df<- merge(coef_df,quantile_df)
bootstrap_df <- bootstrap_df[order(-bootstrap_df$mean_est),]#
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
    all_data <- all_data[,!(names(all_data) %in% c("H3k27me3.10000"))]
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
    paste(tissue,"started forLiver bootstrap",sep=" ")
    plm <- Sys.time()
    for (i in 1:100){
        sample_sites_train_cur <- sample(sample_sites_train,size=length(sample_sites_train),replace=TRUE)
        cur_data <- all_data[sample_sites_train_cur,]
        model<- glm(mutation_status~., data=cur_data,family="binomial")

        cur_coef_df <- as.data.frame(coef(model))
        colnames(cur_coef_df) <- c(paste("value_",i,sep=""))
        cur_coef_df<- tibble::rownames_to_column(cur_coef_df, "name") # https://stackoverflow.com/questions/29511215/convert-row-names-into-first-column
        coef_df <- merge(coef_df,cur_coef_df)
    }
    paste(tissue,"forLiver bootstrap took",Sys.time()-plm,"s/mins/hs, you guess",sep=" ")

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
                      "\nboostrap mean est:",round(filter(bootstrap_df, name==cur_pred)$mean_est,digits=5), "significant?:", filter(bootstrap_df, name==cur_pred)$significant)
        hist(cur_data, main=title,breaks=20)
    }
    dev.off()
    #removing the 0s doesnt make any sense for standardized predictors 
}   
                          
