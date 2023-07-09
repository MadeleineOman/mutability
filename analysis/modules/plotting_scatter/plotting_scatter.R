args = commandArgs(trailingOnly=TRUE)
model_name = args[1]
bin_size_desiredMin = as.integer(args[2])
tissue = args[3]
tissue_predOn = args[4]
if ("_fullModel" %in% args){ fullModel = TRUE
}else{fullModel = FALSE }
if ("_equiv_toLowest" %in% args){ equiv_toLowest = TRUE
}else{equiv_toLowest = FALSE }
if ("_noTriplets" %in% args){exclude_triplet = TRUE
}else{ exclude_triplet = FALSE }
if ("_noCpG" %in% args){exclude_CpG = TRUE
}else{exclude_CpG = FALSE }
if ("_noTCX_CCX" %in% args){ exclude_TCX_CCX = TRUE
}else{ exclude_TCX_CCX = FALSE }



# model_name = "model8"
# bin_size_desiredMin = 400
# tissue = "germline"
# tissue_predOn = "germline"
# fullModel=TRUE
# equiv_toLowest=FALSE
# exclude_triplet = FALSE
# exclude_CpG= FALSE

tmp_pathToFiles = ""


model_desc_modify = ""

if (fullModel ==TRUE){
    model_desc_modify = paste(model_desc_modify,"_fullModel",sep="")}

if (equiv_toLowest==TRUE){
    model_desc_modify = paste(model_desc_modify,"_equiv_toLowest",sep="")}

if (exclude_CpG==TRUE){
    model_desc_modify = paste(model_desc_modify,"_noCpG",sep="")}
if (exclude_triplet==TRUE){
    model_desc_modify = paste(model_desc_modify,"_noTriplets",sep="")} #an arbitrary variable i use to convert b/t the directory of the creator notebook (play mode) and the home directory (snakemake executible)


# 
library(ggplot2)
library(stringr)

#getting the tissue values
# tissue = str_extract(prob_df_file_name,"/(liver|skin|blood|germline)")
# tissue= str_remove(tissue,"/")
# tissue_predOn = str_extract(prob_df_file_name,"on_(liver|skin|blood|germline)")
# tissue_predOn = str_remove(tissue_predOn, "on_")

#improting the files 
prob_df_file_path = paste(tmp_pathToFiles,"data/",tissue,"/dataframes/",model_name,"/",tissue,"_on_",tissue_predOn,"_ProbabilityDf",model_desc_modify,".csv",sep="")
probs_df <- read.csv(prob_df_file_path)



error_output_file = paste(tmp_pathToFiles,"analysis/",tissue,"/plots/",model_name,"/",tissue,"_on_",tissue_predOn,"_createPlot_textOutput",model_desc_modify,".txt",sep="")

mse = round(sum((probs_df$glm_probs - probs_df$mutation_status)*(probs_df$glm_probs - probs_df$mutation_status))/nrow(probs_df),3)
string_to_print = paste("model mse predicting directly on mutations is",mse)
cat(string_to_print,file=error_output_file,sep="\n")
mae = round(sum(abs(probs_df$glm_probs - probs_df$mutation_status))/nrow(probs_df),3)
string_to_print = paste("model mae predicting directly on mutations is",mae)
cat(string_to_print,file=error_output_file,sep="\n",append=TRUE)



#sorting df by proabbility 
probs_df_sorted <- probs_df[order(probs_df$glm_probs, decreasing = FALSE),]

#printing some info about the df 
string_to_print = paste("skin on skin: ",sum(probs_df$mutation_status),"muts,",round(sum(probs_df$mutation_status)/nrow(probs_df),2),"ratio",sep=" ")
cat(string_to_print,file=error_output_file,sep="\n",append=TRUE)

#making the indexing vectors for plotting ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#creating the bins 
nbin = floor(nrow(probs_df_sorted)/bin_size_desiredMin) #=2
extra_leftover = nrow(probs_df_sorted)-nbin*bin_size_desiredMin #=3
extra_toAddToBins = floor(extra_leftover/nbin)#=1
bin_size = bin_size_desiredMin+extra_toAddToBins
extra_to_add_first_bin = nrow(probs_df_sorted)-(bin_size_desiredMin+extra_toAddToBins)*nbin
stopifnot(extra_to_add_first_bin<=nbin)
string_to_print = paste("to comply with exactly equal bin sizes, sites",bin_size,"to",bin_size+extra_to_add_first_bin,"are ommitted out of",nrow(probs_df),"total sites, with",nbin,"bins sized",bin_size,sep=" ")
cat(string_to_print,file=error_output_file,sep="\n",append=TRUE)

#making the vector for plotting 
vector_indexing = c(1, bin_size+extra_to_add_first_bin)
for (i in 1:(nbin-1)){#adjust the finish value so the vector finishes at the nrow of the df 
    vector_indexing <- append(vector_indexing, vector_indexing[i+1]+bin_size)
}
stopifnot(nrow(probs_df_sorted)==tail(vector_indexing,n=1)) #check that the vector finishes at the end of the dataframe

#making the plotting df ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plotting_df <- NULL
for (x in vector_indexing){
    x_end = x+bin_size
    average_predicted_probs = sum(probs_df_sorted$glm_probs[x:x_end])/bin_size
    predicted_probs_sd = sd(probs_df_sorted$glm_probs[x:x_end])
    predicted_probs_sterr = predicted_probs_sd/sqrt(bin_size)
    stderr_lower = average_predicted_probs-predicted_probs_sterr
    stderr_upper = average_predicted_probs+predicted_probs_sterr
    average_proportion_mutations = sum(probs_df_sorted$mutation_status[x:x_end])/bin_size
    row = c(average_predicted_probs,stderr_lower,stderr_upper,average_proportion_mutations)
    plotting_df <-rbind(plotting_df,row)
}
colnames(plotting_df) <- c("average_predicted_probs", "stderr_lower", "stderr_upper","average_proportion_mutations")
rownames(plotting_df) <- NULL
plotting_df <- as.data.frame(plotting_df)

#plotting ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ggplot(plotting_df, aes(x=average_predicted_probs, y=average_proportion_mutations))+
    geom_point()+

    scale_y_continuous(breaks=seq(0,1,0.1)) +
    scale_x_continuous(breaks=seq(0,1,0.1)) +
    theme_light() +
    labs(
        x = 'Predicted mutability bin',
        y = 'Number of mutations / total sites'
        ) +
    theme(
        axis.text = element_text(size = 20, family = 'Helvetica', color = 'black'),
        axis.title = element_text(size = 20, family = 'Helvetica')
        )   +
    geom_abline(intercept=0, slope=1, col="blue")+
    geom_errorbarh(plotting_df, mapping=aes(xmin=stderr_lower, xmax =stderr_upper, y = average_proportion_mutations))+
    ggtitle(paste(tissue,"model on",tissue_predOn,sep=" ")) +
    theme(plot.title = element_text(size = 30, face = "bold",hjust = 0.5))
filename = paste(tmp_pathToFiles,"analysis/",tissue,"/plots/",model_name,"/scatter_",tissue,"_on_",tissue_predOn, model_desc_modify,".pdf",sep="")
ggsave(filename)

#printing some information to file 
string_to_print = paste(round(sum(probs_df$mutation_status)/nrow(probs_df),2),"mutation Ratio",sep=" ")
cat(string_to_print,file=error_output_file,sep="\n",append=TRUE)
string_to_print = paste(nrow(probs_df)," totalSites_",sep=" ")
cat(string_to_print,file=error_output_file,sep="\n",append=TRUE)
string_to_print = paste(bin_size,"_sizedBinsNoOverlap",sep=" ")
cat(string_to_print,file=error_output_file,sep="\n",append=TRUE)
string_to_print = "errorBars = standard error"
cat(string_to_print,file=error_output_file,sep="\n",append=TRUE)


string_to_print = paste("accuracy of model on bins is",1-mean(abs(plotting_df$average_predicted_probs - plotting_df$average_proportion_mutations),na.rm=TRUE),sep=" ")
cat(string_to_print,file=error_output_file,sep="\n",append=TRUE)
string_to_print = paste("MAE of the model on the bins is",mean(abs(plotting_df$average_predicted_probs - plotting_df$average_proportion_mutations),na.rm=TRUE),sep=" ")
cat(string_to_print,file=error_output_file,sep="\n",append=TRUE)

#making the scatter linear model and printing detaisl to file 
fit<-lm(average_proportion_mutations~average_predicted_probs,data=plotting_df)
string_to_print = paste("r-squared is",summary(fit)$r.squared,sep=" ")
cat(string_to_print,file=error_output_file,sep="\n",append=TRUE)
string_to_print = paste("adjusted r-squared is",summary(fit)$adj.r.squared,sep=" ")
cat(string_to_print,file=error_output_file,sep="\n",append=TRUE)