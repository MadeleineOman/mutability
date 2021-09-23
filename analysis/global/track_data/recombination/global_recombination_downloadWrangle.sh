### recombination download wrangle 

# - downloaded from: 
#     - manually, from the UCSC table broswer 
#     - hg18 > mapping and seqeuncing > recomb rate > recombRate 
#     - options: 
#         - genome-wide 
#         - output format = all fields from seelcted table  
#         - output file name = UCSCtableBrowser_hg18_mappingAndSequencing_recombinationRate_recombRate_manualDownload.txt 
#         - filetype returned = palian text format 
#     - submit, file ownload window opens and svae it to /desktop/research/projects/mutability/analysis/global/recombination/
#     - upload to data/global/track_data/recombination/
# - download date: 
#     - lol it depends when i finally do the thing 
# - genome build 
#     - hg18 
# - tissue 
#     - NA
# - lab method:  
#     - The deCODE genetic map was created at deCODE Genetics and is based on 5,136 microsatellite markers for 146 families with a total of 1,257 meiotic events. For more information on this map, see Kong, et al., 2002
#     - The recombination rate track represents calculated sex-averaged rates of recombination based on either the deCODE, Marshfield, or Genethon genetic maps."..." Female- and male-specific recombination rates, as well as rates from the Marshfield and Genethon maps, can also be displayed by choosing the appropriate filter option on the track description page.
# - file format & head: 
#     - OG: 
#         - seems in bed format, though the spacing seems suscpicious . also lots of columns, 
#         - #chrom	chromStart	chromEnd	name	decodeAvg	decodeFemale	decodeMale	marshfieldAvg	marshfieldFemale	marshfieldMale	genethonAvg	genethonFemale	genethonMale
#         - chr1	0	1000000	recombRate	1.11892	1.24759	0.990242	0	0	0	0	0	0
# - what the data represents: 
#     - Each base is assigned the recombination rate calculated by assuming a linear genetic distance across the immediately flanking genetic markers. The recombination rate assigned to each 1 Mb window is the average recombination rate of the bases contained within the window.

#rename the file 
cat analysis/global/track_data/recombination/UCSCtableBrowser_hg18_mappingAndSequencing_recombinationRate_recombRate_manualDownload.txt > data/global/track_data/recombination/recombination.bed

#bgzip the manually uploaded file 
bgzip data/global/track_data/recombination/recombination.bed

#tabix indexing 
tabix -p bed data/global/track_data/recombination/recombination.bed.gz