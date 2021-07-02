
#BASH SCRIPT 

#used for naming 
species = "hsapiens"
project = "mutability"

#mkdirs : general 
mkdir -p /research/projects/$species/$project/{data-analysis-raw_data-conda}
mkdir -p /scratch/projects/$species/$project/{data-analysis-raw_data}
#mkdirs: specific to this project 
mkdir -p /research/projects/$species/$project/raw_data/mutation_datasets/{dsmnc_somatic_muts}
mkdir -p /research/projects/$species/$project/raw_data/predictors/{global-germline-blood}

#make the sym link to the scratch folder that will be holding the data 
cd /research/projects/$species/$project/ #idk why cd here, it's what rbovb dies in his script so why not 
ln -s  /scratch/projects/$species/$project/data /research/projects/$species/$project/data/

#make the conda env 
conda create --prefix /research/projects/$species/$project/conda/$project
conda config --append envs_dirs /research/projects/$species/$project/conda/

#Git stuff
git init
touch .gitignore
echo "./raw_data/" >>.gitignore
echo "./data/" >>.gitignore

git add project_setup.sh
