### global repeats download wrangle 

# - downloaded from: 
#     - manually, from the UCSC table broswer 
#     - hg18 > Variation and Repeats > Simple repeats > simple repeat  
#     - options: 
#         - genome-wide 
#         - output format = all fields from seelcted table  
#         - output file name = UCSCtableBrowser_hg18_variationAndREepeats_simpleRepeats_simpleRepeat.txt
#         - filetype returned = gzipped 
#     - submit, file ownload window opens and svae it to /desktop/research/projects/mutability/data/repeats/
#     - upload to analysis/global/track_data/repeats/
# - download date: 
#     - lol it depends when i finally do the thing 
# - genome build 
#     - hg18 
# - tissue 
#     - NA
# - lab method:  
#     - This track displays simple tandem repeats (possibly imperfect repeats) located by Tandem Repeats Finder (TRF) which is specialized for this purpose. These repeats can occur within coding regions of genes and may be quite polymorphic. Repeat expansions are sometimes associated with specific diseases.
# - file format & head: 
#     - tab delimited
#     - #bin	chrom	chromStart	chromEnd	name	period	copyNum	consensusSize	perMatch	perIndel	score	A	C	G	T	entropy	sequence
#     - 585	chr1	0	468	trf	6	77.2	6	95	3	789	33	51	0	15	1.43	TAACCC
# - what the data represents: 
#     - basically yes/no if there is a repeat there (also the repeat itself, bases & copy #) 

#sort the file and put it in the data folder 
zcat analysis/global/track_data/repeats/UCSCtableBrowser_hg18_variationAndREepeats_simpleRepeats_simpleRepeat.txt.gz | sort -k2,2 -k3,3n > data/global/track_data/repeats/repeats.bed

#bgzip and tabix 
bgzip data/global/track_data/repeats/repeats.bed
tabix -s 2 -b 3 -e 4 data/global/track_data/repeats/repeats.bed.gz