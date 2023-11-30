# **Comparing the predictors of mutability among healthy human tissues inferred from mutations in single cell genome data**
### by Madeleine Oman and Rob W Ness

This is the repository that contains all the analyses and code used in the above publication. Snakemake was used as a project manager, and the workflow can be visualized in the following documents: 
- rulegraph.pdf : the basic bones of the whole analysis, from the webscraping all the way to figure generation
- dag.pdf : the whole workflow. basically the "rulegraph" workflow but for every tissue
- filegraph.pdf : the workflow but with the linking files explicitely given between steps. 

the "Snakefile" file dictates the wrokflow, and contains all the rules and requirements explicitely. 

