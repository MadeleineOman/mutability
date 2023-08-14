import pandas as pd
import numpy as np
import sys


# model_name = "model9"
# model_desc = "_equiv_toLowest"
model_name = sys.argv[1]
model_desc = sys.argv[2]


tmp_path_to_files = ""

MAE_df ={}
r_squared_df = {}

tissues = ["blood","germline","liver","skin"]
comps_toDo = ["blood_on_germline","blood_on_liver","blood_on_skin", "skin_on_liver","skin_on_germline","liver_on_germline"]
for tissue in tissues: 
    curtissue_MAEs = []
    curtissue_r_squareds = []
    for tissue_predOn in tissues: 
        if tissue+"_on_"+tissue_predOn in comps_toDo:  
            for line in open("{f}analysis/global/plots/{m}/tStatScatter_{t}_on_{p}_textOutput{d}.txt".format(f=tmp_path_to_files,m=model_name,t=tissue,p=tissue_predOn,d=model_desc)).readlines()[0:2]:
                if "mean absolute error" in line: 
                    cur_MAE = round(float(line.split()[-1].rstrip("\n")),3)
                    curtissue_MAEs.append(cur_MAE)
                if "r-squared " in line: 
                    cur_r_squared = round(float(line.split()[-1].rstrip("\n")),3)
                    curtissue_r_squareds.append(cur_r_squared)
            
        else: 
            curtissue_MAEs.append(np.nan)
            curtissue_r_squareds.append(np.nan)
    #print(curtissue_MAEs, curtissue_r_squareds)
    MAE_df[tissue] = curtissue_MAEs
    r_squared_df[tissue] = curtissue_r_squareds
    
MAE_df = pd.DataFrame.from_dict(MAE_df,orient="index")
MAE_df.columns = ["on_blood","on_germline","on_liver","on_skin"]
r_squared_df = pd.DataFrame.from_dict(r_squared_df,orient="index")
r_squared_df.columns = ["on_blood","on_germline","on_liver","on_skin"]
    
r_squared_df.to_csv("{t}data/global/dataframes/{m}/tStat_pairwise_r_squared_df{d}.csv".format(t=tmp_path_to_files,m=model_name,d=model_desc))
MAE_df.to_csv("{t}data/global/dataframes/{m}/tStat_pairwise_MAE_df{d}.csv".format(t=tmp_path_to_files,m=model_name,d=model_desc))