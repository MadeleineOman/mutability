library(dplyr)
library(stringi)
library(stats)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
model_name = args[1]
amount_to_red_by = 100
tmp_file_path = ""

#load data 
blood_df <- read.csv(paste(tmp_file_path,"data/blood/dataframes/",model_name,"/blood_forLiver_all_data_readyForPrediction.csv",sep=""))
germline_df <- read.csv(paste(tmp_file_path,"data/germline/dataframes/",model_name,"/germline_forLiver_all_data_readyForPrediction.csv",sep=""))
skin_df <- read.csv(paste(tmp_file_path,"data/skin/dataframes/",model_name,"/skin_forLiver_all_data_readyForPrediction.csv",sep=""))
liver_df <- read.csv(paste(tmp_file_path,"data/liver/dataframes/",model_name,"/liver_all_data_readyForPrediction.csv",sep=""))

#creationg the column that will differntiate the tissues 
blood_df$tissue <- "blood"
germline_df$tissue <- "germline"
skin_df$tissue <- "skin"
liver_df$tissue <- "liver"
#combining tissues
all_data <- rbind(blood_df,liver_df)
all_data <- rbind(all_data,germline_df)
all_data <- rbind(all_data,skin_df)
#factorizing columns 
all_data$tissue <- as.factor(all_data$tissue)
all_data$annotation<- as.factor(all_data$annotation)
all_data$CpGisland <- as.factor(all_data$CpGisland)
all_data$mutation_status <- as.factor(all_data$mutation_status)
all_data$Chromosome <- as.factor(all_data$Chromosome)
all_data$triplet <- as.factor(all_data$triplet)
#removing the rowname column 
all_data <- all_data[,!(names(all_data) %in% c("X"))]

#MAKING THE EQUIV DF 
#finding the min nrow to cut down by 
tissue_nrows = c(nrow(blood_df),nrow(germline_df),nrow(skin_df),nrow(liver_df))
min_nrow <- min(tissue_nrows)
rm(skin_df)
rm(germline_df)
rm(blood_df)
rm(liver_df)
#cutting down the sections to make sure i have all the same # rows 
equiv_blood = (filter(all_data,tissue == "blood")[sample(nrow(filter(all_data,tissue == "blood")), min_nrow), ])
equiv_germline = (filter(all_data,tissue == "germline")[sample(nrow(filter(all_data,tissue == "germline")), min_nrow), ])
equiv_skin = (filter(all_data,tissue == "skin")[sample(nrow(filter(all_data,tissue == "skin")), min_nrow), ])
equiv_liver = (filter(all_data,tissue == "liver")[sample(nrow(filter(all_data,tissue == "liver")), min_nrow), ])
#combine the equivalent dfs 
equiv_data <- rbind(equiv_blood,equiv_germline)
equiv_data <- rbind(equiv_data,equiv_skin)
equiv_data <- rbind(equiv_data,equiv_liver)

model_descs <- c("allData_mutsAndNonMuts","allData_muts","allData_NonMuts","equivData_mutsAndNonMuts","equivData_muts","equivData_NonMuts")
for (model_desc in model_descs){
    if(grepl("allData",model_desc)){
        data_to_use <- all_data
    }else{
        data_to_use <- equiv_data
    }
    if(grepl("mutsAndNonMuts",model_desc)){
        data_to_use <- data_to_use #nothing to do! 
    }else if(grepl("NonMuts",model_desc)){
        data_to_use <- data_to_use %>% filter(mutation_status == "0")
        data_to_use<- data_to_use[,!(names(data_to_use) %in% c("mutation_status"))]
    }else {
        data_to_use <- data_to_use %>% filter(mutation_status == "1")
        data_to_use<- data_to_use[,!(names(data_to_use) %in% c("mutation_status"))]
    }
        
    #create the predictor and reponse for the model input. any NA OMit? can we use the same indexies? 
    predictor_matrix = model.matrix(Chromosome~., data_to_use)[,-1]  #got this from the book, idk . removes NA coloumn and (takes out) the reponse coloum cool cool cool 
    stopifnot(nrow(data_to_use)==nrow(predictor_matrix))

    #perfomr the pca 
    pca <- prcomp(predictor_matrix, scale=TRUE)

    #summarise the components for the scrern plot 
    pca_var <- pca$sdev^2
    pca_var_percent <- round(pca_var/sum(pca_var)*100,1)
    jpeg(paste(tmp_file_path,"analysis/global/plots/",model_name,"/pca_scree_plot_",model_desc,".jpeg",sep=""))
    barplot(pca_var_percent, main="Scree plot",ylab="percent variation explained",xlab="principal component")
    dev.off() 

    #create the pca plotting data frame --> using the first two components 
    pca_data <- data.frame(Tissue=data_to_use$tissue, X=pca$x[,1], Y=pca$x[,2])
    pca_data_reduced =  pca_data[sample(nrow(pca_data), nrow(pca_data)/amount_to_red_by), ]
    #create the coodinates for the mean of each tissue 
    germline_x_mean <- mean(filter(pca_data,Tissue == "germline")$X)
    germline_y_mean <- mean(filter(pca_data,Tissue == "germline")$Y)
    germline_x_stderr <- sd(filter(pca_data,Tissue == "germline")$X)/sqrt(min_nrow)
    germline_y_stderr <- sd(filter(pca_data,Tissue == "germline")$Y)/sqrt(min_nrow)
    skin_x_mean <- mean(filter(pca_data,Tissue == "skin")$X)
    skin_y_mean <- mean(filter(pca_data,Tissue == "skin")$Y)
    skin_x_stderr <- sd(filter(pca_data,Tissue == "skin")$X)/sqrt(min_nrow)
    skin_y_stderr <- sd(filter(pca_data,Tissue == "skin")$Y)/sqrt(min_nrow)
    blood_x_mean <- mean(filter(pca_data,Tissue == "blood")$X)
    blood_y_mean <- mean(filter(pca_data,Tissue == "blood")$Y)
    blood_x_stderr <- sd(filter(pca_data,Tissue == "blood")$X)/sqrt(min_nrow)
    blood_y_stderr <- sd(filter(pca_data,Tissue == "blood")$Y)/sqrt(min_nrow)
    liver_x_mean <- mean(filter(pca_data,Tissue == "liver")$X)
    liver_y_mean <- mean(filter(pca_data,Tissue == "liver")$Y)
    liver_x_stderr <- sd(filter(pca_data,Tissue == "liver")$X)/sqrt(min_nrow)
    liver_y_stderr <- sd(filter(pca_data,Tissue == "liver")$Y)/sqrt(min_nrow)

    #plotting 
    ggplot(data=pca_data_reduced, aes(x=X, y=Y))+ 
        geom_point(aes(colour = factor(Tissue)))+
        xlab(paste("PC1 - ",pca_var_percent[1],"%",sep=" "))+
        ylab(paste("PC2 - ",pca_var_percent[2],"%",sep=" "))+
        theme_bw()+
        #add mean point and standard error for skin
        geom_segment(aes(x=skin_x_mean, y=skin_y_mean-skin_y_stderr, xend=skin_x_mean,yend=skin_y_mean+skin_y_stderr),color="purple")+
        geom_segment(aes(x=skin_x_mean-skin_x_stderr, y=skin_y_mean, xend=skin_x_mean+skin_x_stderr,yend=skin_y_mean),color="purple")+
        geom_point(aes(x=skin_x_mean,y=skin_y_mean),colour="purple",shape="diamond",size = 5)+
        #add mean point and standard error for blood
        geom_segment(aes(x=blood_x_mean, y=blood_y_mean-blood_y_stderr, xend=blood_x_mean,yend=blood_y_mean+blood_y_stderr),color="red")+
        geom_segment(aes(x=blood_x_mean-blood_x_stderr, y=blood_y_mean, xend=blood_x_mean+blood_x_stderr,yend=blood_y_mean),color="red")+
        geom_point(aes(x=blood_x_mean,y=blood_y_mean),colour="red",shape="diamond",size = 5)+
        #add mean point and standard error for germline
        geom_segment(aes(x=germline_x_mean, y=germline_y_mean-germline_y_stderr, xend=germline_x_mean,yend=germline_y_mean+germline_y_stderr),color="green")+
        geom_segment(aes(x=germline_x_mean-germline_x_stderr, y=germline_y_mean, xend=germline_x_mean+germline_x_stderr,yend=germline_y_mean),color="green")+
        geom_point(aes(x=germline_x_mean,y=germline_y_mean),colour="green",shape="diamond",size = 5)+
        #add mean point and standard error for liver
        geom_segment(aes(x=liver_x_mean, y=liver_y_mean-liver_y_stderr, xend=liver_x_mean,yend=liver_y_mean+liver_y_stderr),color="blue")+
        geom_segment(aes(x=liver_x_mean-liver_x_stderr, y=liver_y_mean, xend=liver_x_mean+liver_x_stderr,yend=liver_y_mean),color="blue")+
        geom_point(aes(x=liver_x_mean,y=liver_y_mean),colour="blue",shape="diamond",size = 5)
    ggsave(paste(tmp_file_path,"analysis/global/plots/",model_name,"/pca_",model_desc,".pdf",sep=""))
    #saving information about the pricipal compoentns to dataframes 
    pca_means_summary <- pca_data %>% 
    group_by(Tissue) %>% 
    summarize(MeanX=mean(X, na.rm=TRUE),StdErrX=sd(X)/sqrt(min_nrow),MeanY=mean(Y, na.rm=TRUE),StdErrY=sd(Y)/sqrt(min_nrow))
    write.csv(pca_means_summary,paste(tmp_file_path,"data/global/pca/",model_name,"/pca_",model_desc,"_measnSummaryDf.csv",sep=""), row.names = FALSE)
    all_components_df <- data.frame(pca$rotation)
    pc1_df <- data.frame(sort((pca$rotation[,1]), decreasing=TRUE))
    pc2_df <- data.frame(sort((pca$rotation[,1]), decreasing=TRUE))
    pcAll_df <- data.frame(pca$rotation)
    write.csv(pc1_df,paste(tmp_file_path,"data/global/pca/",model_name,"/pca_",model_desc,"_pc1.csv",sep=""), row.names = TRUE)
    write.csv(pc2_df,paste(tmp_file_path,"data/global/pca/",model_name,"/pca_",model_desc,"_pc2.csv",sep=""), row.names = TRUE)
    write.csv(pcAll_df,paste(tmp_file_path,"data/global/pca/",model_name,"/pca_",model_desc,"_pcAll.csv",sep=""), row.names = TRUE)
    }
