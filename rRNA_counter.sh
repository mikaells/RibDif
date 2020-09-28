#!/bin/bash

####
#A program to analyse the 16s variation in each genome of a genus
#Mikael Lenz Strube & Mikkel Lindegaard
#14-09-2020
####

#Input sanitation
if [ "$#" -lt 2 ]; then
	echo -e "\nUsage is\nrRNA_counter \n\t-g|--genus <genus>\n\t[-c|--clobber true/false].\n\n"
	exit;
fi

#setting global variable, will be overwritten i argument is present

while :
do
 case "$1" in
	-h | --help)
		echo -e "\nUsage is\nrRNA_counter \n\t-g|--genus <genus>\n\t[-c|--clobber true/false].\n\n"
		exit 0
		;;
	-g | --genus)
		genus="$2"
		shift 2
		;;
	-c | --clobber)
		clobber=$2
		shift 2
		;;
	--)
		shift
		break
		;;
	-*)
		echo -e "Error: Unexpected option $1"
		exit 4
		;;
	*)
		break
		;;
 esac
done

echo -e "\n***rRNA_counter running on $genus***\n\n"

if [[ $clobber = true  ]]
then
	echo -e "Removing old run of $genus."
	if [[ -d $genus ]]
	then
		rm -r $genus
		echo -e "\n"
	else
		echo -e "\t$genus-folder does not exist, ignoring --clobber\n\n"
	fi
elif [[ -d $genus ]]
then
	echo -e "$genus-folder already exists, run again with --clobber true or find another folder to run in.\n\n"
	exit
fi


echo -e  "Downloading all strains of $genus into $genus/refseq/bacteria/ with ncbi-genome-download.\n\n";
ncbi-genome-download -s 'refseq' -F 'fasta' -l 'complete' -g $genus -o $genus -p 10 bacteria

#checking if download worked
if [ ! -d $suge ]
then
	echo -e "\n\tDownload failed, is $genus a correct genus?\n\n"
	exit
fi

#gunzip all in parallel
echo -e "Gunzipping all files.\n\n"
find $genus/refseq/bacteria/ -name "*gz" | parallel 'gunzip {}'



#renaming fna-files
echo -e "Renaming fastas and adding GCF (for genomes with multiple chromosomes).\n\n"
find $genus/refseq/bacteria/ -name "*fna" | parallel 'sed -i "s/[:,/()]//g; s/[: ]/_/g" {} '
find $genus/refseq/bacteria/ -name "*fna" | parallel ' GCF=$(echo $(basename $(dirname {})));  sed -E -i "s/^>(.*)/>$GCF"_"\1/g" {} '

#run barrnap
echo -e "Finding all rRNA genes longer than 90% of expected length with barrnap.\n\n"
find $genus/refseq/bacteria/ -name "*fna" | parallel ' barrnap --kingdom "bac" --quiet --threads 1 --reject 0.90 -o "{.}.rRNA" {}'  > barrnap.log 2>&1

#fish out 16S
echo -e "Fishing out 16S genes.\n\n"
find $genus/refseq/bacteria/ -name "*rRNA" | parallel 'grep "16S" {} -A1 > {.}.16S'

#renaming headers in 16s
echo -e "Renaming 16S headers.\n\n"
find $genus/refseq/bacteria/ -name "*.16S" | parallel ' sed -E -i "s/^>16S_rRNA::(.*):.*/>\1_/g" {} ; awk -i inplace "/^>/ { \$0=\$0"_"++i }1" {}  '

#run splitting 16S, cant make it work in parallel
echo -e "Splitting 16S.\n\n"
original_PWD=$PWD
for folder in $genus/refseq/bacteria/*; do
	mkdir $folder/indiv_16S_dir/
	cd  $folder/indiv_16S_dir/
	awk '/^>/{f=(++i-1)".fna"}1 {print > f}' ../*.16S
	cd $original_PWD
done

#calculating ANI for each genome
echo -e "Calculating mismatches for each genome with ANI.\n\n"
ls -d $genus/refseq/bacteria/*/indiv_16S_dir/ | parallel 'average_nucleotide_identity.py -i {} -o {}/../ani/'

echo -e "Alligning 16S genes within genomes with muscle and builing trees with fastree.\n\n"
find $genus/refseq/bacteria/ -name "*.16S" | parallel 'muscle -in {} -out {.}.16sAln -quiet; sed -i "s/[ ,]/_/g" {.}.16sAln; fasttree -quiet -nopr -gtr -nt {.}.16sAln > {.}.16sTree '

#Summarizing data
echo -e "Summarizing data into $genus/$genus-summary.csv.\n\n"
echo -e "GCF\tGenus\tSpecies\t#16S\tMean\tSD\tMin\tMaxtTotalDiv" > $genus/$genus-summary.csv
ls -d $genus/refseq/bacteria/* | parallel Rscript ~/rRNA_counter/run16sSummary.R {}/ani/ANIm_similarity_errors.tab {}/*16sAln {}/16S_div.pdf {}/*fna {}/*16sTree >> $genus/$genus-summary.csv


wait;

mkdir $genus/full
find $genus/refseq/bacteria/ -name "*16S" -exec cat {}  \; > $genus/full/$genus.16S

echo -e "Alligning all 16S rRNA genes with mafft and building tree with fasttree.\n\n"
mafft --auto --quiet --thread 20 $genus/full/$genus.16S > $genus/full/$genus.aln
fasttree -quiet -nopr -gtr -nt $genus/full/$genus.aln > $genus/full/$genus.tree

echo -e "Making amplicons with in_silico_pcr.\n\n"
mkdir $genus/amplicons/
/home/common/scripts/in_silico_PCR.pl -s $genus/full/$genus.16S -a CCTACGGGNGGCNGCAG    -b GACTACNNGGGTATCTAATCC -m -i > $genus/amplicons/$genus-V3V4.summary 2> $genus/amplicons/$genus-V3V4.temp.amplicons
/home/common/scripts/in_silico_PCR.pl -s $genus/full/$genus.16S -a AGAGTTTGATCCTGGCTCAG -b CGGTTACCTTGTTACGACTT  -m -i > $genus/amplicons/$genus-V1V9.summary 2> $genus/amplicons/$genus-V1V9.temp.amplicons

#renaming headers
seqkit replace --quiet -p "(.+)" -r '{kv}' -k $genus/amplicons/$genus-V3V4.summary $genus/amplicons/$genus-V3V4.temp.amplicons > $genus/amplicons/$genus-V3V4.amplicons
seqkit replace --quiet -p "(.+)" -r '{kv}' -k $genus/amplicons/$genus-V1V9.summary $genus/amplicons/$genus-V1V9.temp.amplicons > $genus/amplicons/$genus-V1V9.amplicons


#deleting old amplicon files
rm $genus/amplicons/$genus-V3V4.temp.amplicons
rm $genus/amplicons/$genus-V1V9.temp.amplicons

echo -e "Alligning all amplicons with mafft and building tree with fasttree.\n\n"
mafft --auto --quiet --thread 20 $genus/amplicons/$genus-V3V4.amplicons > $genus/amplicons/$genus-V3V4.aln
fasttree -quiet -nopr -gtr -nt $genus/amplicons/$genus-V3V4.aln > $genus/amplicons/$genus-V3V4.tree

echo -e "Making gene summary file for CLC import.\n\n"
Rscript ~/rRNA_counter/Format16STreesForCLC.R $genus/full/$genus.tree $genus/full/$genus-CLC.csv

echo -e "Making amplicon summary file for CLC import.\n\n"
Rscript ~/rRNA_counter/Format16STreesForCLC.R $genus/amplicons/$genus-V3V4.tree $genus/amplicons/$genus-V3V4-CLC.csv

echo -e "Done.\n\n"
