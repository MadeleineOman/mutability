import re 
import pandas as pd

import sys

model_name = sys.argv[1]
model_desc = sys.argv[2]
# model_name = "model8"
# model_desc = "_fullModel"
bin_size = 400

tmp_path_to_files = ""

MAE_df ={}
MAE_df_onBins = {}
MSE_df = {}
r_squared_df = {}
adj_r_squared_df = {}

tissues = ["blood","germline","liver","skin"]
for tissue in tissues: 
    curtissue_MSEs =[]
    curtissue_MAEs = []
    curtissue_MAEs_onBins = []
    curtissue_r_squareds = []
    curtissue_adj_r_squareds = []
    for tissue_predOn in tissues: 
        for line in open("{f}analysis/{t}/plots/{m}/{t}_on_{p}_createPlot_textOutput{d}.txt".format(f=tmp_path_to_files,m=model_name,t=tissue,p=tissue_predOn,d=model_desc)).readlines():
            if "model mse predicting directly on mutations is" in line: 
                cur_MSE = line.split()[-1].rstrip("\n")
                curtissue_MSEs.append(cur_MSE)
            if "model mae predicting directly on mutations is" in line: 
                cur_MAE = line.split()[-1].rstrip("\n")
                curtissue_MAEs.append(cur_MAE)
            if re.search(r"^r-squared is",line): # need regex for this one so doesnt match the adj r squared line
                cur_r_squared = round(float(line.split()[-1].rstrip("\n")),3)
                curtissue_r_squareds.append(cur_r_squared)
            if "adjusted r-squared" in line: 
                cur_adj_r_squared = round(float(line.split()[-1].rstrip("\n")),3)
                curtissue_adj_r_squareds.append(cur_adj_r_squared)
            if "MAE of the model on the bins" in line: 
                cur_MAE_onBins = round(float(line.split()[-1].rstrip("\n")),3)
                curtissue_MAEs_onBins.append(cur_MAE_onBins)
                
    MSE_df[tissue] = curtissue_MSEs
    MAE_df[tissue] = curtissue_MAEs
    MAE_df_onBins[tissue] = curtissue_MAEs_onBins
    r_squared_df[tissue] = curtissue_r_squareds
    adj_r_squared_df[tissue] = curtissue_adj_r_squareds

    
MSE_df = pd.DataFrame.from_dict(MSE_df,orient="index")
MSE_df.columns = ["on_blood","on_germline","on_liver","on_skin"]
MAE_df_onBins = pd.DataFrame.from_dict(MAE_df_onBins,orient="index")
MAE_df_onBins.columns = ["on_blood","on_germline","on_liver","on_skin"]
MAE_df = pd.DataFrame.from_dict(MAE_df,orient="index")
MAE_df.columns = ["on_blood","on_germline","on_liver","on_skin"]
r_squared_df = pd.DataFrame.from_dict(r_squared_df,orient="index")
r_squared_df.columns = ["on_blood","on_germline","on_liver","on_skin"]
adj_r_squared_df = pd.DataFrame.from_dict(adj_r_squared_df,orient="index")
adj_r_squared_df.columns = ["on_blood","on_germline","on_liver","on_skin"]

adj_r_squared_df.to_csv("{t}data/global/dataframes/{m}/adj_r_squared_df{d}_{b}.csv".format(t=tmp_path_to_files,m=model_name,d=model_desc,b=bin_size))
r_squared_df.to_csv("{t}data/global/dataframes/{m}/r_squared_df{d}_{b}.csv".format(t=tmp_path_to_files,m=model_name,d=model_desc,b=bin_size))
MAE_df.to_csv("{t}data/global/dataframes/{m}/MAE_df{d}_{b}.csv".format(t=tmp_path_to_files,m=model_name,d=model_desc,b=bin_size))
MSE_df.to_csv("{t}data/global/dataframes/{m}/MSE_df{d}_{b}.csv".format(t=tmp_path_to_files,m=model_name,d=model_desc,b=bin_size))
MAE_df_onBins.to_csv("{t}data/global/dataframes/{m}/MAE_df_onBins{d}_{b}.csv".format(t=tmp_path_to_files,m=model_name,d=model_desc,b=bin_size))