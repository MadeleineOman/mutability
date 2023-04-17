
library(dplyr)
library(stringr)
library(ggplot2)
library(patchwork)
library(tidyr)

args = commandArgs(trailingOnly=TRUE)
model_name = args[1]
if ("_equiv_toLowest" %in% args){
    equiv_toLowest=TRUE
    print("equiv_toLowest=TRUE")}
if ("_exclude_CpG" %in% args){exclude_CpG=TRUE}
if ("_exclude_triplet" %in% args){exclude_triplet=TRUE}
if ("_fullModel" %in% args){fullModel=TRUE}

# model_name = "model8"
# equiv_toLowest=TRUE
# exclude_triplet = FALSE
# exclude_CpG= FALSE

tmp_file_path = ""


model_desc_modify = ""
if (equiv_toLowest==TRUE){model_desc_modify = paste(model_desc_modify,"_equiv_toLowest",sep="")}
if (fullModel==TRUE){model_desc_modify = paste(model_desc_modify,"_fullModel",sep="")}
if (exclude_CpG==TRUE){model_desc_modify = paste(model_desc_modify,"_noCpG",sep="")}
if (exclude_triplet==TRUE){model_desc_modify = paste(model_desc_modify,"_noTriplets",sep="")}


#load in the data 
blood_coefs <- read.csv(paste(tmp_file_path,"data/blood/dataframes/",model_name,"/blood_forLiver_coefDF",model_desc_modify,".csv",sep=""))
germline_coefs <- read.csv(paste(tmp_file_path,"data/germline/dataframes/",model_name,"/germline_forLiver_coefDF",model_desc_modify,".csv",sep=""))
skin_coefs <- read.csv(paste(tmp_file_path,"data/skin/dataframes/",model_name,"/skin_forLiver_coefDF",model_desc_modify,".csv",sep=""))
liver_coefs <- read.csv(paste(tmp_file_path,"data/liver/dataframes/",model_name,"/liver_coefDF",model_desc_modify,".csv",sep=""))

#rename columns tp include the tissue in the colanme (for later merging)
for(i in 2:ncol(blood_coefs)) names(blood_coefs)[names(blood_coefs) == colnames(blood_coefs)[i]] = paste("blood_",colnames(blood_coefs)[i],sep="")
for(i in 2:ncol(germline_coefs)) names(germline_coefs)[names(germline_coefs) == colnames(germline_coefs)[i]] = paste("germline_",colnames(germline_coefs)[i],sep="")
for(i in 2:ncol(liver_coefs)) names(liver_coefs)[names(liver_coefs) == colnames(liver_coefs)[i]] = paste("liver_",colnames(liver_coefs)[i],sep="")
for(i in 2:ncol(skin_coefs)) names(skin_coefs)[names(skin_coefs) == colnames(skin_coefs)[i]] = paste("skin_",colnames(skin_coefs)[i],sep="")

#merge data together 
all_coefs <- merge(blood_coefs,germline_coefs)#https://www.statmethods.net/management/merging.html
all_coefs <- merge(all_coefs,liver_coefs)
all_coefs <- merge(all_coefs,skin_coefs)


#create col that categorizes predictors into types 
all_coefs<-all_coefs %>%
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
#create col that categorizes color for plotting 
all_coefs<-all_coefs %>%
    mutate(color_name = case_when(        
        str_detect(type, 'sequence') ~ "#7CAE00",
        str_detect(type, 'triplet') ~ "#C77CFF",
        str_detect(type, 'tissue_specific') ~ '#00BFC4',
        str_detect(type, 'global') ~ '#F8766D',
        str_detect(type, 'intercept') ~ "darkgrey",
    ))

#create conversion vector for predictor type and ggploting colour                           
colors_ggplot <- c(sequence = "#7CAE00", triplet = "#C77CFF", global ='#F8766D', tissue_specific = '#00BFC4', intercept = "darkgrey")

#gwet only significant predictors 
all_coefs_sign =filter(all_coefs,(blood_p_val<0.01 | germline_p_val<0.01| liver_p_val<0.01| skin_p_val<0.01))

#selecting columns for simpler data processing 
all_coefs_violin = all_coefs[,(names(all_coefs) %in% c("name","blood_value","liver_value","germline_value","skin_value","type","color_name"))]
all_coefs_violin_sign = all_coefs_sign[,(names(all_coefs_sign) %in% c("name","blood_value","liver_value","germline_value","skin_value","type","color_name"))]
all_tStat_violin = all_coefs[,(names(all_coefs) %in% c("name","blood_t_stat","liver_t_stat","germline_t_stat","skin_t_stat","type","color_name"))]
all_tStat_violin_sign = all_coefs_sign[,(names(all_coefs_sign) %in% c("name","blood_t_stat","liver_t_stat","germline_t_stat","skin_t_stat","type","color_name"))]

#getting the average fro all tissues --> need this to sort the violin plot 
all_coefs_violin$av_coef_val <- rowMeans(subset(all_coefs_violin, select = c("blood_value", "germline_value","liver_value","skin_value")), na.rm = TRUE)
all_coefs_violin_sign$av_coef_val <- rowMeans(subset(all_coefs_violin_sign, select = c("blood_value", "germline_value","liver_value","skin_value")), na.rm = TRUE)
all_tStat_violin$av_tStat_val <- rowMeans(subset(all_tStat_violin, select = c("blood_t_stat", "germline_t_stat","liver_t_stat","skin_t_stat")), na.rm = TRUE)
all_tStat_violin_sign$av_tStat_val <- rowMeans(subset(all_tStat_violin_sign, select = c("blood_t_stat", "germline_t_stat","liver_t_stat","skin_t_stat")), na.rm = TRUE)

#wide to long data http://www.cookbook-r.com/Manipulating_data/Converting_data_between_wide_and_long_format/
all_coefs_violin <- gather(all_coefs_violin,source,coef_value,blood_value:skin_value)
all_coefs_violin_sign <- gather(all_coefs_violin_sign,source,coef_value,blood_value:skin_value)
all_tStat_violin <- gather(all_tStat_violin,source,tStat_value,blood_t_stat:skin_t_stat)
all_tStat_violin_sign <- gather(all_tStat_violin_sign,source,tStat_value,blood_t_stat:skin_t_stat)


#PLOTTING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#T STATISTIC PLOTTING

p<-ggplot(all_tStat_violin, aes(x =reorder(name, -av_tStat_val), y = tStat_value,color=type)) + #http://www.sthda.com/english/wiki/ggplot2-violin-plot-quick-start-guide-r-software-and-data-visualization
# p<-ggplot(all_coefs_violin, aes(x =name, y = coef_value,color=type)) + #http://www.sthda.com/english/wiki/ggplot2-violin-plot-quick-start-guide-r-software-and-data-visualization
    geom_abline(intercept=0, slope=0, col="grey",linetype="dashed")+
    geom_violin()+
    geom_point()+
    xlab(paste("predictors",sep=" "))+
    ylab("predictor value distribution across all tissues")+
    theme(text = element_text(size=14),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    scale_color_manual(values=colors_ggplot)
ggsave(plot = p, width = 15, height = 10, dpi = 300, filename = paste(tmp_file_path,"analysis/global/plots/",model_name,"/coefViolinPlot_all_tStatDev",model_desc_modify,"_all.pdf",sep=""))

p<-ggplot(all_tStat_violin_sign, aes(x =reorder(name, -av_tStat_val), y = tStat_value,color=type)) + #http://www.sthda.com/english/wiki/ggplot2-violin-plot-quick-start-guide-r-software-and-data-visualization
# p<-ggplot(all_coefs_violin, aes(x =name, y = coef_value,color=type)) + #http://www.sthda.com/english/wiki/ggplot2-violin-plot-quick-start-guide-r-software-and-data-visualization
    geom_abline(intercept=0, slope=0, col="grey",linetype="dashed")+
    geom_violin()+
    geom_point()+
    xlab(paste("Significant predictors",sep=" "))+
    ylab("T statistic value distribution across all tissues")+
    theme(text = element_text(size=14),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    scale_color_manual(values=colors_ggplot)
ggsave(plot = p, width = 15, height = 10, dpi = 300, filename = paste(tmp_file_path,"analysis/global/plots/",model_name,"/coefViolinPlot_all_tStatDev",model_desc_modify,"_onlySign.pdf",sep=""))


#COEFICIENT PLOTS  
p<-ggplot(all_coefs_violin, aes(x =reorder(name, -av_coef_val), y = coef_value,color=type)) + #http://www.sthda.com/english/wiki/ggplot2-violin-plot-quick-start-guide-r-software-and-data-visualization
# p<-ggplot(all_coefs_violin, aes(x =name, y = coef_value,color=type)) + #http://www.sthda.com/english/wiki/ggplot2-violin-plot-quick-start-guide-r-software-and-data-visualization
    geom_abline(intercept=0, slope=0, col="grey",linetype="dashed")+
    geom_violin()+
    geom_point()+
    xlab(paste("predictors",sep=" "))+
    ylab("predictor value distribution across all tissues")+
    theme(text = element_text(size=14),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    scale_color_manual(values=colors_ggplot)
ggsave(plot = p, width = 15, height = 10, dpi = 300, filename = paste(tmp_file_path,"analysis/global/plots/",model_name,"/coefViolinPlot_all_coefDev",model_desc_modify,"_all.pdf",sep=""))

p<-ggplot(all_coefs_violin_sign, aes(x =reorder(name, -av_coef_val), y = coef_value,color=type)) + #http://www.sthda.com/english/wiki/ggplot2-violin-plot-quick-start-guide-r-software-and-data-visualization
# p<-ggplot(all_coefs_violin, aes(x =name, y = coef_value,color=type)) + #http://www.sthda.com/english/wiki/ggplot2-violin-plot-quick-start-guide-r-software-and-data-visualization
    geom_abline(intercept=0, slope=0, col="grey",linetype="dashed")+
    geom_violin()+
    geom_point()+
    xlab(paste("predictors",sep=" "))+
    ylab("predictor value distribution across all tissues")+
    theme(text = element_text(size=14),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    scale_color_manual(values=colors_ggplot)
ggsave(plot = p, width = 15, height = 10, dpi = 300, filename = paste(tmp_file_path,"analysis/global/plots/",model_name,"/coefViolinPlot_all_coefDev",model_desc_modify,"_onlySign.pdf",sep=""))

