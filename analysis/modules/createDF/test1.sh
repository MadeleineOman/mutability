python analysis/modules/createDF/createDF.py germline model7 '[0,100,10000]' _allTissueSpecTracks;
grep 'buffer' data/germline/dataframes/model7/predictorDf_allTissueSpecTracks.txt >> data/germline/dataframes/model7/predictorDf_allTissueSpecTracks_errorlog.txt;
grep 'discord' data/germline/dataframes/model7/predictorDf_allTissueSpecTracks.txt >> data/germline/dataframes/model7/predictorDf_allTissueSpecTracks_errorlog.txt;
grep -v 'discord' data/germline/dataframes/model7/predictorDf_allTissueSpecTracks.txt > data/germline/dataframes/model7/predictorDf_allTissueSpecTracks_noDiscord.txt;
grep -v 'buffer' data/germline/dataframes/model7/predictorDf_allTissueSpecTracks_noDiscord.txt >  data/germline/dataframes/model7/predictorDf_allTissueSpecTracks.txt