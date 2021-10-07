# RibDif
RibDif evaluates the differences in 16S ribosomal RNA genes within a genus or species.

Takes an input genus and returns a large set of analysed files, the main being heatmaps comparing amplicons of 16S rRNA genes. Default amplicons are V3V4 and V1V9, but 
custom primers can be optionally be supplied.

This can be useful for determining if species resolution is possible within certain genera using common metataxonomic approaches. Often they are not.

Ribdif runs in Linux and installs through conda (https://docs.conda.io/en/latest/) in a couple of minutes. It has been tested on Ubuntu 18.04 and 20.04, other distros 
may or may not work

## Installation

#Clone from github

git clone https://github.com/mikaells/RibDif/

#Create conda environment from yml file

conda env create -f RibDif/RibDif.yml

#Activate conda environment 

conda activate RibDif

#Set permissions of files to be executable 

chmod 755 RibDif/*

#silence the parallel citation message if need be (and do with the message what you wish)

parallel --citation   


Note: in_silico_PCR.pl is borrowed from https://github.com/egonozer/in_silico_pcr


## Usage
#Run analysis by as follows 

path-to-RibDif/RibDif.sh -g $genus

#An example run with Ruegeria as target genus:

~/RibDif/RibDif.sh -g Ruegeria

#Can run with species as well, make sure to use quotes

~/RibDif/RibDif.sh -g "Mycoplasma bovis"


Running Ruegeria will less than a minute, while large genera like Bacillus might take more than an hour as it has 7000+ genomes with 10+ 16S genes.

In some genera, unspecific R-errors are thrown. It is caused by some genomes having no rRNA genes, but causes no problems downstream. How a completed genome can have no 
rRNA genes is of course cause of some concern. Occasionally, an error occors when people have used weird characters in their genome names. I have tried parsing all cases, but people are suprisingly creative when it comes to naming schemes.


## Full set of options:

 -g|--genus     
 #the genus you want
 
[-c|--clobber]  
#Delete previous run

[-p|--primers]   
#Path to custom primer-file, must be a tab seperated file with name, forward and reverse primers. See the file 'v3v4.primers' 

[-a|--ANI]      
#ANI is off by default, turn on if you care about individual genomes

[-f|--frag]     
#off by default, full genomes are required for detecting all 16S-genes, use with caution 

[-m|--msa]     
#off by default, make multiple sequence alignment and trees.  

[-i|--id]       
#1 (100% identity) as default, if the final evaluation should be done at amplicons clustered at another identity, 
like .99 (pretty experimental right now). Does not cluster at the genome level, so beware.

[-t|--threads]  
#default all, number of threads to use


## Output

The program generates a new directory named after the genus in question. Within that is a summary file and three sub-directories. Many files are only present with all options turned on.

The main files of interest for amplicon metataxonomic analysis is the amplicon/ folder, in which you will find analysis for all all primers. The heatmap pdfs is the most important ones, as they show which species overlap, which is also summarized in the 'XX-overlap-summary.txt' files.


amplicons/ #files for amplicons
  
    $genus-$region.amplicons # all 16S amplicons genes  
    $genus-$region.aln # alligned from 16S amplicons (alligned by mafft)  
    $genus-$region.tree # tree from allignment  
    $genus-$region-meta.csv # overview of each amplicon from each genome and the cluster they belong to
    $genus-$region.uc # clustering file from vsearch
    $genus-$region-heatmap.pdf # heatmap of 1) clusters found in each genome and 2) genomes that would overlap with given primers
    $region-clusters/ #amplicons clustered into unique clusters by vsearch, in fasta format
    $genus-$region-overlap-summary.txt # analysis summary
    
full/ # files for full length 16S analysis
  
    $genus.16S # all 16S genes  
    $genus.aln # alligned from 16S genes (alligned by mafft)  
    $genus.tree # tree from allignment  
    $genus-meta.csv # metadata for each sequence  

refseq/bacteria/ #contains a folder for each genome, named after refseq-identifier

    GCF*.fna # the original genome file, headers are renamed for downstream compatibility
    GCF*.fna.fai # indexed fna file  
    GCF*.rRNA # rRNA file generated by barrnap  
    GCF*.16S # only the 16S file from .rRNA file  
    GCF*.16sAln # alligned .16S file (by muscle)  
    GCF*.16sTree # tree from .16sAln file (by fasttree)  
    16S_div.pdf # pdf for 16S tree and diversity across genes  
    indiv_16S_dir/ # folder for individual 16S files, for pyani compatibility  
    x.fna #16s files  
    ani/ # files generated from pyani


$genus-summary.csv # contains summary statistics for each genome. Columns 5-8 are in nucleotide mismatches between the genomes' 16S genes, TotalDiv is the cummulative shannon index
 
    #GCF     Genus   Species #16S    Mean    SD      Min     Max   TotalDiv 


  
The tree files are just newick-files and can be used in any treeviewer or R. The meta.csv file is formatted for CLC, but is just a table and can be used for other treeviewers supporting annotation.

# Example output

V3V4 amplicons generated from the Pseudoalteromonas genus results in 51 unique alleles, and these are identical multiple species such as the allele found in agarivorans/atlantica/carrageenovora/espejiana/issachenkonii/tetraodonis, which hence cannot be seperated using V3V4 sequencing. Moreover, several genomes have multiple distinct alleles, which will artificially overinflate diversity estimates. Additionally, some genomes have alleles shared with different species, such as a genome of P. arctica which has alleles seperately shared with P. translucida and P. nigrifaciens. V3V4 amplicons from this genome would erronously suggest that three different species are present.

![My image](https://github.com/mikaells/RibDif/blob/master/img/Pseudoalteromonas-V3V4_clusterdistri.png)

A more concise view of which genomes overlap. Evidently, multiple species cannot be distinguished from V3V4 sequencing.
![My image](https://github.com/mikaells/RibDif/blob/master/img/Pseudoalteromonas-V3V4_confusionmat.png)


RibDif gives the following overlap summary which may be easier to parse for large genera:

    Summary:
    
    Genomes: 56
    	Named: 35
    	Non-named: 21

    Named species: 24
    
    29 of 56 (51.79%) genomes have multiple alleles.

    19 of 24 (79.17%) species have at least one overlap.

    The following species overlap:
	agarivorans/atlantica/carrageenovora/espejiana/issachenkonii/sp./tetraodonis
	nigrifaciens/sp./translucida
	piscicida/rubra/sp./viridis
	arabiensis/shioyasakiensis
	arctica/sp./translucida
	arctica/translucida
	piratica/spongiae
	donghaensis/sp.
	ruthenica/sp.
	marina/sp.

