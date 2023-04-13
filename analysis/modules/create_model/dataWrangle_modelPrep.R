library(dplyr)
library(stringr)
library(stringi)
library(corrplot)
library(rlang)
library(ggplot2)


# tissue = "germline"
# model_name = "model8"
equiv_toLowest=FALSE
exclude_CpG=FALSE 
exclude_triplet=FALSE
exclude_TCX_CCX = FALSE
NA_omit=TRUE
mappability_col_keep = FALSE

args = commandArgs(trailingOnly=TRUE)
tissue = args[1]
model_name = args[2]
if ("_equiv_toLowest" %in% args){ equiv_toLowest = TRUE
}else{equiv_toLowest = FALSE }
if ("_noTriplets" %in% args){exclude_triplet = TRUE
}else{ exclude_triplet = FALSE }
if ("_noCpG" %in% args){allTissueSpecTracks = TRUE
}else{allTissueSpecTracks = FALSE }
if ("_noTCX_CCX" %in% args){ exclude_TCX_CCX = TRUE
}else{ exclude_TCX_CCX = FALSE }
if ("_allTissueSpecTracks" %in% args){allTissueSpecTracks = TRUE
}else{allTissueSpecTracks = FALSE }

tmp_file_path = ""
options(warn=1)

model_desc_modify = ""

#import the data 
input_filePath = paste(tmp_file_path,"data/",tissue,"/dataframes/",model_name,"/predictorDf",model_desc_modify,".txt",sep="")
all_data <- read.table(input_filePath, header = TRUE,sep="\t")
if (equiv_toLowest==TRUE){
    blood_nrow = nrow(read.table(paste(tmp_file_path,"data/blood/dataframes/",model_name,"/predictorDf",model_desc_modify,".txt",sep=""),header = TRUE,sep="\t"))
    all_data <-all_data[sample(nrow(all_data), blood_nrow), ]
    model_desc_modify = paste(model_desc_modify,"_equiv_toLowest",sep="")}



if (exclude_CpG==TRUE){
    all_data <-all_data[!str_detect(all_data$triplet,"CG"),]
    model_desc_modify = paste(model_desc_modify,"_noCpG",sep="")}
if (exclude_TCX_CCX==TRUE){
    model_desc_modify = paste(model_desc_modify,"_noTCX_CCX",sep="")}

# print(summary(as.factor(all_data$triplet)))
if (exclude_triplet==TRUE){
    all_data <- all_data[,!(names(all_data) %in% c('triplet'))]
    model_desc_modify = paste(model_desc_modify,"_noTriplets",sep="")}
# print(colnames(all_data))

#setting up the error output file 
error_output_file = paste(tmp_file_path,"data/",tissue,"/objects/",model_name,"/",tissue,"_create_model_text_output",model_desc_modify,".txt",sep="")
cat("output for the create_model notebook",file=error_output_file,sep="\n")

#chnaging the 1-surround to 0 surround (--> due to a mess up, only blood )
# if (tissue=="blood"){
#     all_data <- all_data %>% 
#   rename(Apercent.0=Apercent.1,Gpercent.0=Gpercent.1,Cpercent.0=Cpercent.1,Tpercent.0=Tpercent.1,
#          methylation_coverage.0=methylation_coverage.1,methylation_precent.0=methylation_precent.1,
#          H3k27.0=H3k27.1,H3k27me3.0=H3k27me3.1,H3k4me1.0=H3k4me1.1,H3k4me3.0=H3k4me3.1,H3k36me3.0=H3k36me3.1,
#          Transcription.0=Transcription.1,recombination.0=recombination.1,Repeats.0=Repeats.1,DNAse.0=DNAse.1,laminB1.0=laminB1.1)
# }

#remove the protein binding flag (if it exists in the version of the predictor data: not existant in later versions)
all_data$annotation <- gsub("protein_binding","not_transcribed",all_data$annotation)#there are so few protein binding sites that we may as well omit 

#removing unamppable points "
mappable_mut_summary <- all_data %>% 
    group_by(mappability,mutation_status)%>%
    summarise(count=n())
unmapped_muts = (filter(mappable_mut_summary, (mappability == "not")&(mutation_status == 1))$count)
unmapped_nonMuts = (filter(mappable_mut_summary, (mappability == "not")&(mutation_status == 0))$count)
string_to_print = paste(tissue,":",unmapped_muts," mutations and",unmapped_nonMuts,"nonMut sites removed by mappability filter out of ",nrow(all_data),"sites",sep=" ")
cat(string_to_print,file=error_output_file,sep="\n",append=TRUE)
if (mappability_col_keep==FALSE){
    all_data <- filter(all_data,mappability=="mappable")
    all_data <- all_data[,!(names(all_data) %in% c("mappability"))]
}else{model_desc_modify=paste(model_desc_modify,"_withMappableCol",sep="")}

#handling gc content"
all_data <- all_data %>%
    mutate(GC_content.0 = Gpercent.0+Cpercent.0) %>%
    mutate(GC_content.100 = Gpercent.100+Cpercent.100) %>% 
    mutate(GC_content.10000 = Gpercent.10000+Cpercent.10000)
all_data <- all_data[,!(names(all_data) %in% c('Apercent.0','Gpercent.0','Cpercent.0','Tpercent.0',
                                       'Apercent.100','Gpercent.100','Cpercent.100','Tpercent.100',
                                       'Apercent.10000','Gpercent.10000','Cpercent.10000','Tpercent.10000'))]

#editing the triplets 
if (exclude_triplet==FALSE){
    all_data$triplet <- toupper(all_data$triplet)
    #filter for rows that dont have NNN as the triplet --> write the details to file 
    cat(paste(nrow(all_data[all_data$triplet == "NNN",]), "rows removed due to N in triplet, ",sep=" "),file=error_output_file,sep="\n",append=TRUE)
    all_data <- all_data[all_data$triplet != "NNN",]
    cat(paste(nrow(all_data),"rows left",sep=" "),file=error_output_file,sep="\n",append=TRUE)
}


#create df to test coverage in the methylation download notebook 
bases = c("A","T","C","G")
c_trips = c()
for (b1 in bases){
    for (b3 in bases){
        c_trips <- append(c_trips,paste(b1,"C",b3,sep=""))
    }
}


all_data <- all_data[,!str_detect(colnames(all_data),"methyl")]
# model_desc_modify = paste(model_desc_modify,"_noMethyl",sep="")





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
if (exclude_triplet==FALSE){
    rc_removeAG <- function(dna){
        middle_base = substr(dna, 2, 2)
        if(middle_base %in% c("A","G")){
            dna <- stri_reverse(chartr("acgtACGT", "tgcaTGCA", dna))}
        return(dna)
    }#substring slicing https://www.johnmyleswhite.com/notebook/2009/02/25/text-processing-in-r/
    all_data$triplet <- unlist(lapply(as.character(all_data$triplet),rc_removeAG)) #need unlistto turn the list into a vector 
}

if (exclude_TCX_CCX==TRUE){
    all_data <-all_data[!str_detect(all_data$triplet,"[TC]C[ATCG]"),]}

#check for problems with levels / column tye etc. 
if (exclude_triplet==FALSE){
    all_data$triplet <- as.character(all_data$triplet) #need to convert to char so the as.factor properly reducs to 64 levels after the removal of 'triplets" longer than 3
    string_to_print = paste(nrow(all_data[nchar(as.character(all_data$triplet))!=3,]),"rows removed due to triplet larger than 3 in length ",sep=" ")
    cat(string_to_print,file=error_output_file,sep="\n",append=TRUE)
    all_data <- all_data[nchar(all_data$triplet)==3,] #make sure only including triplets 
    string_to_print=paste(nrow(all_data)," rows left", sept = " ")
    cat(string_to_print,file=error_output_file,sep="\n",append=TRUE)
}

#editing the annotation colun --> ignored into non_transcribed 
all_data$annotation <- gsub("ignored","not_transcribed",all_data$annotation)
all_data$annotation <- gsub("protein_binding","not_transcribed",all_data$annotation)#there are so few protein binding sites that we may as well omit 
#rerplacing with "not_transcribed" is fine as they are the last level before non_transcribed in the annotation module 

#factorizing the columns 
all_data$mutation_status <- as.factor(all_data$mutation_status)
if (exclude_triplet==FALSE){all_data$triplet <- as.factor(as.character(all_data$triplet))}
all_data$Chromosome <- as.factor(all_data$Chromosome)
all_data$annotation <- as.factor(all_data$annotation)
all_data$CpGisland <- as.factor(all_data$CpGisland)



#checking that the factor variables are correct (mutation status, triplets, chroms) --> will raise error if not true 
stopifnot(length(levels(all_data$mutation_status)) == 2)
if (exclude_triplet==FALSE){
    if(exclude_CpG==TRUE&exclude_TCX_CCX==FALSE){
        stopifnot(length(levels(all_data$triplet))==28)}
    if(exclude_CpG==FALSE&exclude_TCX_CCX==TRUE){
        stopifnot(length(levels(all_data$triplet))==24)}
    if(exclude_CpG==TRUE&exclude_TCX_CCX==TRUE){
        stopifnot(length(levels(all_data$triplet))==22)}
    if(exclude_CpG==FALSE&exclude_TCX_CCX==FALSE){
        stopifnot(length(levels(all_data$triplet))==32)}}
stopifnot(length(levels(all_data$Chromosome))==22)
stopifnot(length(levels(all_data$CpGisland))==3)


#if germline, then edit out the female and make the male the default ( i used to include both female and male) 
if (tissue == "germline"){
    all_data <- select(all_data,-matches("female")) 
    colnames(all_data)<- str_replace_all(colnames(all_data), "_male", "")
}

#na omit and checking how many rows were removed 
string_to_print = paste(nrow(all_data)- nrow(na.omit(all_data)), "lost to NA values, ", nrow(na.omit(all_data))," rows remain after",sep=" ")
cat(string_to_print,file=error_output_file,sep="\n",append=TRUE)
if (NA_omit ==TRUE){
    all_data <- na.omit(all_data)
    }

#writing frequencies and counts for the categories of factor predictors 
#get he factor names 
factor_names <- names(Filter(is.factor,all_data))
factor_names <-factor_names[ !factor_names == 'mutation_status']
#for loop to print for each factor 
predictor_freqs_count_filePath <- paste(tmp_file_path,"data/",tissue,"/objects/",model_name,"/",tissue,"_create_model_predFreqs",model_desc_modify,".csv",sep="")
cat("predictor,category,frequency,n_sites",file=predictor_freqs_count_filePath,sep="\n")
for (predictor_name in factor_names){
#     string_to_print<-(paste(predictor_name," mutation distributions",sep=" "))
#     cat(string_to_print,file=predictor_freqs_count_filePath,sep="\n",append=TRUE)
    test_df<-as.data.frame(table(all_data[,"mutation_status"], all_data[,predictor_name]))#https://stackoverflow.com/questions/14596420/how-to-get-value-by-column-name-in-r
    colnames(test_df) <- c("mut_stat",predictor_name,"Freq")
    for (cur_annot in levels(all_data[,predictor_name])){
        n_sites = sum(test_df[test_df[predictor_name]==cur_annot,]$Freq)
        n_muts <- as.integer(test_df[test_df[predictor_name]==cur_annot&test_df["mut_stat"]=="1"][3])#the third colun is the one with the count 
        string_to_print <-((paste(predictor_name,cur_annot,n_muts/n_sites,n_sites,sep=",")))
        cat(string_to_print,file=predictor_freqs_count_filePath,sep="\n",append=TRUE)}}


#print("about to do the distributions")

#distributions ~~~~~~~~~~~~~~~~~~~~~
#getting only numeric 
num_cols <- colnames(all_data[,unlist(lapply(all_data, is.numeric), use.names = FALSE)])#https://stackoverflow.com/questions/5863097/selecting-only-numeric-columns-from-a-data-frame

pdf(paste(paste(tmp_file_path,"analysis/",tissue,"/plots/",model_name,"/",tissue,"_predictor_distributions",model_desc_modify,".pdf", sep="")))
for (cur_pred in num_cols){
    cur_data = all_data[,cur_pred]
    title = paste(cur_pred,"\nsd:",round(sd(cur_data),digits=2), " mean:", round(mean(cur_data),digits=2), " median:",round(median(cur_data),digits=2)," min,max:",round(min(cur_data),digits=5),",",max(round(cur_data,digits=2)))
    hist(cur_data, main=title,breaks=20)
}
dev.off()
pdf(paste(paste(tmp_file_path,"analysis/",tissue,"/plots/",model_name,"/",tissue,"_predictor_distributions_remove0",model_desc_modify,".pdf", sep="")))
for (cur_pred in num_cols){
    cur_data = filter(all_data,!!sym(cur_pred)>0)[,cur_pred]# get r to handle strign as value https://stackoverflow.com/questions/48219732/pass-a-string-as-variable-name-in-dplyrfilter
    title = paste(cur_pred,"\nsd:",round(sd(cur_data),digits=2), " mean:", round(mean(cur_data),digits=2), " median:",round(median(cur_data),digits=2)," min,max:",round(min(cur_data),digits=5),",",max(round(cur_data,digits=2)))
    hist(cur_data, main=title,breaks=20)
}
dev.off()

#making categorical columns into seperate numerical columns ready for standardizing 
all_data$dummy <- 1#need a dumym column that the model.matrix can remove (as it has to remove a column aparently)
muts_col <- all_data$mutation_status#saving the column of mutations so i can add it back later
if (tissue=="global"){
    tissue_col <- all_data$tissue
    non_num_preds = c("mutation_status","tissue")
}else{non_num_preds = c("mutation_status")}#saving the column of mutations (and tissue, if relevant) so i can add it back later

preds_to_standardize<-!(names(all_data) %in% non_num_preds)
all_data <-(data.frame((model.matrix(dummy~., all_data[, preds_to_standardize])[,-1] )))
all_data$mutation_status <- muts_col
if (tissue=="global"){
    all_data$tissue <- tissue_col}

#print("standardizing")
#STANDARDIZING~~~~~~~~~~~~~~~~~~~~~~~~~~~~
non_num_preds = c("mutation_status")
if (tissue=="global"){
    non_num_preds = c("mutation_status","tissue")
}
preds_to_standardize<-!(names(all_data) %in% non_num_preds)#yes i have to do this again because i chnaged the order of the mutant column 
all_data[, preds_to_standardize]<- (as.data.frame(scale(all_data[, preds_to_standardize],center=TRUE,scale=TRUE)))


#CORRELATION ANALYSIS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#print("creating and saving the corplot ")
# num_data <- na.omit(subset(all_data, select=-c(Chromosome,triplet,annotation,CpGisland))) --> dont need this anymore as correlating ther categoricals too 
num_data <-all_data[,!(names(all_data) %in% non_num_preds)]
num_data$mutation_status <- as.numeric(all_data$mutation_status)
num_cor <- cor(num_data)
filename = paste(tmp_file_path,"analysis/",tissue,"/plots/",model_name,"/",tissue,"_corPlot",model_desc_modify,".pdf", sep="")
pdf(filename)
corrplot(num_cor,method = "color",type="upper",number.cex=0.5,tl.cex = 0.5,main=tissue)
dev.off() 


#print("create the correlation df and pruning for new corplot")
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
filename = paste(tmp_file_path,"data/",tissue,"/dataframes/",model_name,"/",tissue,"_fullModel_corTable_above0.8",model_desc_modify,".csv", sep="")
write.csv(cor_df,filename, row.names = FALSE)
                          
#REMOVING CORRELATED VARIABLES 
correlated_vars = c('H3k27.1','H3k27.0','H3k27.100',
                         'H3k4me1.0','H3k4me1.1','H3k4me1.100',
                         'H3k4me3.0','H3k4me3.1','H3k4me3.100',
                         'H3k27me3.0','H3k27me3.1','H3k27me3.100',
                         'H3k36me3.0','H3k36me3.1','H3k36me3.100',
                         "Transcription.0","Transcription.1","Transcription.100",
                         "recombination.0","recombination.1","recombination.100",
                         "DNAse.0","DNAse.1","DNAse.10000",
                         "GC_content.0",
                         "Repeats.100")
all_data <- all_data[,!(names(all_data) %in% correlated_vars)]
                          

#VIF ANALYSIS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
x_data <- all_data[,!(names(all_data) %in% c("mutation_status","mutation_status1"))]
pred_names <-1:length((colnames(x_data)))
adj_rsqs <-1:length((colnames(x_data)))
vifs <-1:length((colnames(x_data)))
i=1
for (predictor_name in (colnames(x_data))){
    formula = paste(predictor_name,"~.",sep="")#write the formula as a string to feed into the lm function 
    adj_rsq <- summary(lm(formula,data=x_data))$adj.r.squared#extract r sqared 
    vif <- 1/(1-adj_rsq)                   #calculate vif 
    pred_names[i] <- predictor_name        #save the predictor name to the vector 
    adj_rsqs[i] <- adj_rsq                 #save the adjisted r squared to the vector 
    vifs[i] <- vif                         #save the vif to a vecotr 
    i=i+1                                  #iterate the index for saving the current r^2/vif 
    }
inter_cor_vif_df <- data.frame(pred_names,adj_rsqs,vifs)
inter_cor_vif_df<- inter_cor_vif_df[order(-vifs),]

#saving the model variables 
filename =paste(tmp_file_path,"data/",tissue,"/dataframes/",model_name,"/",tissue,"_vif",model_desc_modify,".csv", sep="")
write.csv(inter_cor_vif_df,filename) 

num_data <-all_data[,!(names(all_data) %in% non_num_preds)]
num_data$mutation_status <- as.numeric(all_data$mutation_status)
num_cor <- cor(num_data)
filename = paste(tmp_file_path,"analysis/",tissue,"/plots/",model_name,"/",tissue,"_corPlot_pruned",model_desc_modify,".pdf", sep="")
pdf(filename)
corrplot(num_cor,method = "color",type="upper",number.cex=0.5,tl.cex = 0.5,main=tissue)
dev.off() 
                          
                          
#saving the model variables 
filename =paste(tmp_file_path,"data/",tissue,"/dataframes/",model_name,"/",tissue,"_all_data_readyForPrediction",model_desc_modify,".csv", sep="")
write.csv(all_data,filename)    