args = commandArgs(trailingOnly=TRUE)
# tissue = "blood"
# tissue_predOn = "skin"
# model_name = "model8"
equiv_toLowest=FALSE
exclude_triplet = FALSE
exclude_CpG= FALSE

tissue = args[1]
tissue_predOn = args[2]
model_name = args[3]
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


model_desc_modify = ""
if (equiv_toLowest==TRUE){
    model_desc_modify = paste(model_desc_modify,"_equiv_toLowest",sep="")}

if (exclude_CpG==TRUE){
    model_desc_modify = paste(model_desc_modify,"_noCpG",sep="")}
if (exclude_triplet==TRUE){
    model_desc_modify = paste(model_desc_modify,"_noTriplets",sep="")}

library(dplyr)
library(stringr)
library(ggplot2)
library(patchwork)

#reading in data
if (tissue_predOn =="liver"){ #oif conditional as different files to compare to liver 
    tissue_coefs <- read.csv(paste(tmp_file_path,"data/",tissue,"/dataframes/",model_name,"/",tissue,"_forLiver_coefDF",model_desc_modify,".csv",sep=""))
}else{tissue_coefs <- read.csv(paste(tmp_file_path,"data/",tissue,"/dataframes/",model_name,"/",tissue,"_coefDF",model_desc_modify,".csv",sep=""))}
if (tissue=="liver"){
   tissue_predOn_coefs <- read.csv(paste(tmp_file_path,"data/",tissue_predOn,"/dataframes/",model_name,"/",tissue_predOn,"_forLiver_coefDF",model_desc_modify,".csv",sep=""))
}else{tissue_predOn_coefs <- read.csv(paste(tmp_file_path,"data/",tissue_predOn,"/dataframes/",model_name,"/",tissue_predOn,"_coefDF",model_desc_modify,".csv",sep=""))}


colnames(tissue_coefs) <- c("name",'tissue_est','tissue_stdErr','tissue_zVal','tissue_pVal','tissue_tStat')
colnames(tissue_predOn_coefs) <-c("name",'tissue_predOn_est','tissue_predOn_stdErr','tissue_predOn_zVal','tissue_predOn_pVal','tissue_predOn_tStat')

all_coefs <- merge(tissue_coefs,tissue_predOn_coefs,by="name")#for some reason the dim b/t germ and blood are different so need a merge... 

#create col that categorizes 
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
        
#         str_detect(name, 'Intercept') ~ '(intercept)',
    ))



all_sign_coefs <- filter(all_coefs[order(-all_coefs$tissue_est),],(tissue_pVal<0.05)|(tissue_predOn_pVal<0.05))


colors_ggplot <- c(sequence = "#7CAE00", triplet = "#C77CFF", global ='#F8766D', tissue_specific = '#00BFC4', intercept = "darkgrey")

#plotting
#  
#writing some infor on how wel they match to file 
error_output_file = paste(tmp_file_path,"analysis/global/plots/",model_name,"/coefScatter_",tissue,"_on_",tissue_predOn,"_textOutput",model_desc_modify,".txt",sep="")
MAE = mean(abs(all_coefs$tissue_est - all_coefs$tissue_predOn_est),na.rm=TRUE)
string_to_print = paste("mean absolute error is",round(MAE,4),sep=" ")
cat(string_to_print,file=error_output_file,sep="\n",append=TRUE)
fit<-lm(tissue_est~tissue_predOn_est,data=all_coefs)
string_to_print = paste("r-squared is",summary(fit)$r.squared,sep=" ")
cat(string_to_print,file=error_output_file,sep="\n",append=TRUE)





#plotting

p1<-ggplot(all_sign_coefs, aes(y = tissue_est, x = tissue_predOn_est,label=name)) +
    geom_point(aes(fill=type,color=type) )+
    theme_light()+
#     scale_y_continuous(breaks=seq(-0.75,1.75,0.75)) +
#     scale_x_continuous(breaks=seq(-0.75,1.75,0.75)) +
    geom_abline(intercept=0, slope=1, col="grey",linetype="dashed")+ #http://www.sthda.com/english/wiki/ggplot2-line-types-how-to-change-line-types-of-a-graph-in-r-software
    geom_abline(intercept=0, slope=0, col="cornflowerblue")+
    geom_vline(xintercept=0,col="cornflowerblue") + #http://www.sthda.com/english/wiki/ggplot2-add-straight-lines-to-a-plot-horizontal-vertical-and-regression-lines
    labs(
        x = paste(tissue_predOn," coefficient values",sep=""),
        y = paste(tissue," coefficient values",sep=""),
        color = "Predictor class"#https://stackoverflow.com/questions/14622421/how-to-change-legend-title-in-ggplot
        ) +
    theme(
    axis.text = element_text(size = 20, family = 'Helvetica', color = 'black'),
    axis.title = element_text(size = 20, family = 'Helvetica'),
    legend.title = element_text(size=12, family = 'Helvetica'),
    legend.text = element_text(size = 12, family = 'Helvetica')
    )  +
    scale_color_manual(values = colors_ggplot)+
    geom_text(hjust=0, vjust=0,aes(color=type),size=2)#https://stackoverflow.com/questions/15624656/label-points-in-geom-point
# ggsave(p1, paste(tmp_file_path,"analysis/global/plots/",model_name,"/coefScatter_",tissue,"_on_",tissue_predOn,"_onlySignCoefs_label",model_desc_modify,".pdf",sep=""))
p2<-ggplot(all_sign_coefs, aes(y = tissue_est, x = tissue_predOn_est,label=name)) +
    geom_point(aes(fill=type,color=type) )+
    theme_light()+
#     scale_y_continuous(breaks=seq(-0.75,1.75,0.75)) +
#     scale_x_continuous(breaks=seq(-0.75,1.75,0.75)) +
    geom_abline(intercept=0, slope=1, col="grey",linetype="dashed")+ #http://www.sthda.com/english/wiki/ggplot2-line-types-how-to-change-line-types-of-a-graph-in-r-software
    geom_abline(intercept=0, slope=0, col="cornflowerblue")+
    geom_vline(xintercept=0,col="cornflowerblue") + #http://www.sthda.com/english/wiki/ggplot2-add-straight-lines-to-a-plot-horizontal-vertical-and-regression-lines
    labs(
        x = paste(tissue_predOn," coefficient values",sep=""),
        y = paste(tissue," coefficient values",sep=""),
        color = "Predictor class"#https://stackoverflow.com/questions/14622421/how-to-change-legend-title-in-ggplot
        ) +
    theme(
    axis.text = element_text(size = 20, family = 'Helvetica', color = 'black'),
    axis.title = element_text(size = 20, family = 'Helvetica'),
    legend.title = element_text(size=12, family = 'Helvetica'),
    legend.text = element_text(size = 12, family = 'Helvetica')
    )  +
    geom_pointrange(aes(ymin=tissue_est-tissue_stdErr, ymax=tissue_est+tissue_stdErr, color=type),alpha=0.4)+
    geom_pointrange(aes(xmin=tissue_predOn_est-tissue_predOn_stdErr, xmax=tissue_predOn_est+tissue_predOn_stdErr, color=type),alpha=0.4)+
    scale_color_manual(values = colors_ggplot)
# ggsave(p2,paste(tmp_file_path,"analysis/global/plots/",model_name,"/coefScatter_",tissue,"_on_",tissue_predOn,"_onlySignCoefs",model_desc_modify,".pdf",sep=""))
p<-p1+p2
ggsave(plot = p, width = 18, height = 10, dpi = 200, filename =paste(tmp_file_path,"analysis/global/plots/",model_name,"/coefScatter_",tissue,"_on_",tissue_predOn,"_onlySignCoefs",model_desc_modify,"_comb.pdf",sep=""))


MAE = mean(abs(all_sign_coefs$tissue_est - all_sign_coefs$tissue_predOn_est),na.rm=TRUE)
string_to_print = paste("mean absolute error is for sign coefs only ",round(MAE,4),sep=" ")
cat(string_to_print,file=error_output_file,sep="\n",append=TRUE)
fit<-lm(tissue_est~tissue_predOn_est,data=all_sign_coefs)
string_to_print = paste("r-squared is for sign coefs only",summary(fit)$r.squared,sep=" ")
cat(string_to_print,file=error_output_file,sep="\n",append=TRUE)