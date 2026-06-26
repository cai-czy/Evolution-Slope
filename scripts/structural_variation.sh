#!/bin/bash

############################################################
# Plant-microbiome-genome coupling promotes
# incipient ecological speciation
#
# Structural variation analyses
#
# Software
# MUMmer v4
# SyRI v1.7.1
# plotsr v1.1.1
############################################################


###############################
# Input genome assemblies
###############################

REF=B59
QRY=T166

THREADS=16

mkdir -p syri_results
mkdir -p coords


############################################################
# Structural variation detection
############################################################

for chr in {1..7}
do

echo "Processing chromosome ${chr}H ..."

########################################
# Whole-genome alignment
########################################

nucmer --maxmatch -c 500 -b 500 -l 100 -t ${THREADS} -p chr${chr} ${REF}/${REF}.chr${chr}.fa ${QRY}/${QRY}.chr${chr}.fa


########################################
# Filter alignments
########################################

delta-filter -m -i 90 -l 150 chr${chr}.delta > chr${chr}.filtered.delta


########################################
# Generate coordinate file
########################################

show-coords -THrd chr${chr}.filtered.delta > coords/chr${chr}.coords


########################################
# Detect structural variants
########################################

syri -c coords/chr${chr}.coords -d chr${chr}.filtered.delta -r ${REF}/${REF}.chr${chr}.fa -q ${QRY}/${QRY}.chr${chr}.fa -k --prefix syri_results/chr${chr}

done


############################################################
# Visualization
############################################################

plotsr --sr syri_results/chr1syri.out syri_results/chr2syri.out syri_results/chr3syri.out syri_results/chr4syri.out syri_results/chr5syri.out syri_results/chr6syri.out syri_results/chr7syri.out --genomes genomes.txt -o Structural_variation.pdf


############################################################
# End
############################################################
