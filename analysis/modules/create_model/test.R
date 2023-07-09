args = commandArgs(trailingOnly=TRUE)
tissue = args[1]
model_name = args[2]
if ("_bloodEquiv" %in% args){ equiv_toLowest = TRUE
}else{equiv_toLowest = FALSE }
if ("_noTriplets" %in% args){exclude_triplet = TRUE
}else{ exclude_triplet = FALSE }
if ("_noCpG" %in% args){allTissueSpecTracks = TRUE
}else{allTissueSpecTracks = FALSE }
if ("_noTCX_CCX" %in% args){ exclude_TCX_CCX = TRUE
}else{ exclude_TCX_CCX = FALSE }
if ("_allTissueSpecTracks" %in% args){allTissueSpecTracks = TRUE
}else{allTissueSpecTracks = FALSE }

print(args)
print(equiv_toLowest)
print(allTissueSpecTracks)