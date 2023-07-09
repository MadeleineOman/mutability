library(dplyr)
library(stringr)
library(stringi)
library(corrplot)
library(rlang)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
tissue = args[1]
model_name = args[2]
n_bootstrap = as.integer(args [3])
equiv_toLowest=FALSE
exclude_CpG=FALSE
exclude_triplet=FALSE
exclude_TCX_CCX = FALSE
if ("_equiv_toLowest" %in% args){
    equiv_toLowest=TRUE
    print("equiv_toLowest=TRUE")}
if ("_exclude_CpG" %in% args){exclude_CpG=TRUE}
if ("_exclude_TCX_CCX" %in% args){exclude_TCX_CCX=TRUE}
if ("_exclude_triplet" %in% args){exclude_triplet=TRUE}

# tissue = "liver"
# model_name = "model9"
n_bootstrap = 0


tmp_file_path = ""

model_desc_modify = ""


if (equiv_toLowest==TRUE){
    model_desc_modify = paste(model_desc_modify,"_equiv_toLowest",sep="")}else{
    model_desc_modify = paste(model_desc_modify,"_fullModel",sep="")}
if (exclude_CpG==TRUE){
    all_data <-all_data[!str_detect(all_data$triplet,"CG"),]
    model_desc_modify = paste(model_desc_modify,"_noCpG",sep="")}
if (exclude_TCX_CCX==TRUE){
    model_desc_modify = paste(model_desc_modify,"_noTCX_CCX",sep="")}
if (exclude_triplet==TRUE){
    all_data <- all_data[,!(names(all_data) %in% c('triplet'))]
    model_desc_modify = paste(model_desc_modify,"_noTriplets",sep="")}


filename =paste(tmp_file_path,"data/",tissue,"/dataframes/",model_name,"/",tissue,"_all_data_readyForPrediction",model_desc_modify,".csv", sep="")
all_data <- read.csv(filename)
all_data <- all_data[,!(names(all_data) %in% c("X","X.1"))]

#MAKING A MODEL ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sample_sites_train <- sample (1:nrow(all_data), nrow(all_data)/2) 
sample_sites_test = c(1:nrow(all_data))[-sample_sites_train]
model <- glm(mutation_status~., data=all_data[sample_sites_train,],family="binomial")
               



                          
           
#saving the model variables 
filename = paste(tmp_file_path,"data/",tissue,"/objects/",model_name,"/",tissue,"_model",model_desc_modify,".RData", sep="") #model 
save(model, file=filename)
filename =paste(tmp_file_path,"data/",tissue,"/objects/",model_name,"/",tissue,"_samples_sites_test",model_desc_modify,".RData", sep="")#sample sites test 
save(sample_sites_test , file=filename)
filename =paste(tmp_file_path,"data/",tissue,"/dataframes/",model_name,"/",tissue,"_all_data_readyForPrediction",model_desc_modify,".csv", sep="")
write.csv(all_data,filename)
                          
#checing the assumptions.. though idk how to read 
pdf(paste(tmp_file_path,"analysis/",tissue,"/plots/",model_name,"/",tissue,"_model_diagnostics",model_desc_modify,".pdf", sep=""))
plot(model)
dev.off() 
#checking for linearity one at a time... 
pdf(paste(tmp_file_path,"analysis/",tissue,"/plots/",model_name,"/",tissue,"_model_diagnostics_checkingLinearity",model_desc_modify,".pdf", sep=""))
for (pred_name in colnames(all_data)){
    plot(model$linear.predictors ~ all_data[sample_sites_train,][,pred_name],main=pred_name)
}
dev.off()

#analyze the coeficients~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#first the whole model 
coefs <- coef(summary(model))
coef_df <- as.data.frame(coefs)
colnames(coef_df) <- c("value","std_err","z_val","p_val")
coef_df$t_stat <- coef_df$value/coef_df$std_err
coef_df<- tibble::rownames_to_column(coef_df, "name") # https://stackoverflow.com/questions/29511215/convert-row-names-into-first-column
coef_df_ordered <- coef_df[order(-coef_df$value),]#https://www.statmethods.net/management/sorting.html
filename = paste(tmp_file_path,"data/",tissue,"/dataframes/",model_name,"/",tissue,"_coefDF",model_desc_modify,".csv",sep="")#this sep is for the filename string
write.csv(coef_df_ordered,filename,row.names=FALSE)

coef_df_ordered <- coef_df_ordered %>%
    mutate(type = case_when(        
        str_detect(name, 'percent') ~ 'sequence',
        str_detect(name, 'Chromosome') ~ 'sequence',
        str_detect(name, 'Repeats') ~ 'sequence',
        str_detect(name, 'site') ~ 'sequence',
        str_detect(name, 'annot') ~ 'sequence',
        str_detect(name, 'content') ~ 'sequence',
        str_detect(name, 'CpG') ~ 'sequence',

        str_detect(name, 'triplet') ~ 'triplet',

        str_detect(name, 'DNAse') ~ 'tissue_specific',
        str_detect(name, 'Transcription') ~ 'tissue_specific',
        str_detect(name, 'H3k') ~ 'tissue_specific',
        str_detect(name, 'methyl') ~ 'tissue_specific',

        str_detect(name, 'recomb') ~ 'global',
        str_detect(name, 'lamin') ~ 'global',
        str_detect(name, 'Replication') ~ 'global',
        str_detect(name, 'dist_rep') ~ 'global',

        str_detect(name, 'Intercept') ~ 'intercept',
    ))
stopifnot(nrow(coef_df_ordered[is.na(coef_df_ordered$type),])==0)#making sure all column names have a type

#GGPLOT PLOTTING THE COEF VALUES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#setting the colors for each predictor type 
colors_ggplot <- c(sequence = "#7CAE00", triplet = "#C77CFF", global ='#F8766D', tissue_specific = '#00BFC4', intercept = "darkgrey")# https://stackoverflow.com/questions/17180115/manually-setting-group-colors-for-ggplot2


#bar plot of all the coefs 
ggplot(coef_df_ordered[coef_df_ordered$p_val<0.05,], aes(x = reorder(name, -value), y = value,fill=type))+
    geom_bar(stat="identity")+
    xlab(paste(tissue,"predictors",sep=" "))+
    ylab("Coeficient estimate of all significnat predictors ")+
    #https://www.tutorialspoint.com/how-to-display-negative-labels-below-bars-in-barplot-using-ggplot2-in-r
    geom_text(aes(y=value+0.03*sign(value),label=name),angle = 90,size=2)+#bars next to bar value 
    #geom_text(aes(y=-0.05*sign(mean_est),label=name),angle = 90,size=2)+#labels next to the bottom (0 line) of bars 
    theme(text = element_text(size=7),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#     scale_fill_manual(values=c(orange_ggplot,"darkgrey", green_ggplot, blue_ggplot,purple_ggplot))
    scale_fill_manual(values=colors_ggplot)
 ggsave(paste(tmp_file_path,"analysis/",tissue,"/plots/",model_name,"/",tissue,"_coef_barplot_all",model_desc_modify,".pdf",sep=""))
#bar plot of the top 20 predictors 
coef_df_ordered_top10 <- head(coef_df_ordered[order(-abs(coef_df_ordered$value)),],n=20)
ggplot(coef_df_ordered_top10, aes(x = reorder(name, -abs(value)), y = value,fill=type))+
    geom_bar(stat="identity")+
    xlab(paste(tissue,"predictors",sep=" "))+
    ylab("Coeficient value +- SE")+ 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size=14),axis.text.y = element_text(size=16))+
    theme(axis.title=element_text(size=14))+
    scale_fill_manual(values=colors_ggplot)+
    geom_errorbar(aes(ymin=value-std_err, ymax=value+std_err), width=.2,position=position_dodge(.9),color="darkslategrey")+ #http://www.sthda.com/english/wiki/ggplot2-error-bars-quick-start-guide-r-software-and-data-visualization
    theme(axis.text.x = element_text(face="bold", size=14, angle=45))
ggsave(paste(tmp_file_path,"analysis/",tissue,"/plots/",model_name,"/",tissue,"_coef_barplot_top20",model_desc_modify,".pdf",sep=""))

#NOW FOR THE T STATISTIC 
ggplot(coef_df_ordered[coef_df_ordered$p_val<0.05,], aes(x = reorder(name, -value), y = t_stat,fill=type))+
    geom_bar(stat="identity")+
    xlab(paste(tissue,"predictors ordered by coeficient value",sep=" "))+
    ylab("T statistic estimate of all significnat predictors ")+
    #https://www.tutorialspoint.com/how-to-display-negative-labels-below-bars-in-barplot-using-ggplot2-in-r
    geom_text(aes(y=t_stat+1*sign(t_stat),label=name),angle = 90,size=2)+#bars next to bar value 
    #geom_text(aes(y=-0.05*sign(mean_est),label=name),angle = 90,size=2)+#labels next to the bottom (0 line) of bars 
    theme(text = element_text(size=7),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#     scale_fill_manual(values=c(orange_ggplot,"darkgrey", green_ggplot, blue_ggplot,purple_ggplot))
    scale_fill_manual(values=colors_ggplot)
 ggsave(paste(tmp_file_path,"analysis/",tissue,"/plots/",model_name,"/",tissue,"_tstat_barplot_all",model_desc_modify,".pdf",sep=""))



#CREATING MODELS FOR THE LIVER TISSUE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (tissue != "liver"){
    all_data <- all_data[,!(names(all_data) %in% c("H3k27me3.10000"))]
    #one whole model for 
    model <- glm(mutation_status~., data=all_data[sample_sites_train,],family="binomial")

    #saving the model variables 
    filename = paste(tmp_file_path,"data/",tissue,"/objects/",model_name,"/",tissue,"_forLiver_model",model_desc_modify,".RData", sep="") #model  
    save(model, file=filename)
    filename =paste(tmp_file_path,"data/",tissue,"/objects/",model_name,"/",tissue,"_forLiver_samples_sites_test",model_desc_modify,".RData", sep="")#sample sites test 
    save(sample_sites_test , file=filename)
    filename =paste(tmp_file_path,"data/",tissue,"/dataframes/",model_name,"/",tissue,"_forLiver_all_data_readyForPrediction",model_desc_modify,".csv", sep="")
    write.csv(all_data,filename)

    #analyze the coeficients~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    coefs <- coef(summary(model))
    coef_df <- as.data.frame(coefs)
    colnames(coef_df) <- c("value","std_err","z_val","p_val")
    coef_df$t_stat <- coef_df$value/coef_df$std_err
    coef_df<- tibble::rownames_to_column(coef_df, "name") # https://stackoverflow.com/questions/29511215/convert-row-names-into-first-column
    coef_df_ordered <- coef_df[order(-coef_df$value),]#https://www.statmethods.net/management/sorting.html
    filename = paste(tmp_file_path,"data/",tissue,"/dataframes/",model_name,"/",tissue,"_forLiver_coefDF",model_desc_modify,".csv",sep="")#this sep is for the filename string
    write.csv(coef_df_ordered,filename,row.names=FALSE)

    #MODEL DIAGNOSTICS 
    #checing the assumptions.. though idk how to read 
    pdf(paste(tmp_file_path,"analysis/",tissue,"/plots/",model_name,"/",tissue,"_forLiver_model_diagnostics",model_desc_modify,".pdf", sep=""))
    plot(model)
    dev.off() 
    #checking for linearity one at a time... 
    pdf(paste(tmp_file_path,"analysis/",tissue,"/plots/",model_name,"/",tissue,"_forLiver_model_diagnostics_checkingLinearity",model_desc_modify,".pdf", sep=""))
    for (pred_name in colnames(all_data)){
        plot(model$linear.predictors ~ all_data[sample_sites_train,][,pred_name],main=pred_name)
    }
    dev.off()    
}   
   