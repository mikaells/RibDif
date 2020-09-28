# rRNA_counter
Count and compare 16S rRNA genes within a genome

Takes an input genus and returns a large set of analysed files, the main being phylogenetic trees of full length 16S rRNA genes and common amplicons thereof.


## Installation

#Clone from github
git clone https://github.com/mikaells/rRNA_counter/

#Create conda environment form yml file
conda env create -f rRNA_counter/rRNA_counter.yml

#Activate conda environment
conda activate rRNA

#Set permissions of files to be executable
chmod 755 rRNA_counter/*

#Run analysis by as follows
path-to-rRNA_counter/rRNA_counter.sh -g GENUS
