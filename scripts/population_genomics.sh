#!/bin/bash

############################################################
# Plant-microbiome-genome coupling promotes
# incipient ecological speciation
#
# Population genomic analyses
#
# Software
# VCFtools v0.1.17
# pixy v2.0.0
# PopLDdecay v3.43
# SMC++ v1.15.5
# fastsimcoal2 v2.6.0.3
############################################################


############################################################
# Input files
############################################################

VCF=barley.snp.vcf.gz
POP1=Basalt.txt
POP2=Terrarossa.txt
THREADS=16

############################################################
# FST
# Sliding window = 1 Mb
############################################################

vcftools --gzvcf ${VCF} --weir-fst-pop ${POP1} --weir-fst-pop ${POP2} --fst-window-size 1000000 --out FST_1Mb


############################################################
# Tajima's D
# Sliding window = 250 kb
############################################################

vcftools --gzvcf ${VCF} --TajimaD 250000 --out TajimaD_250kb

############################################################
# Nucleotide divergence (dXY)
# Sliding window = 1 Mb
############################################################

pixy --vcf ${VCF} --populations populations.txt --window_size 1000000 --output_folder pixy_output --output_prefix dXY


############################################################
# Linkage disequilibrium (LD) decay
############################################################

PopLDdecay -InVCF ${VCF} -OutStat LD_decay -MaxDist 500


############################################################
# Historical effective population size
############################################################

smc++ estimate --cores ${THREADS} 6.5e-9 SMC_output sample.smc.gz


############################################################
# Demographic inference
############################################################

easySFS.py -i ${VCF} -p population.txt -a -o easySFS

fsc26 -t ES.tpl -e ES.est -n 100000 -L 40 -m -M -q -c 20 -MUTRATE 6.5e-9 -B 20 


############################################################
# End
############################################################
