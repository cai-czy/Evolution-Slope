#!/bin/bash

############################################################
# Plant-microbiome-genome coupling promotes
# incipient ecological speciation
#
# Genome annotation pipeline
#
# Software
# RepeatModeler v2.0.7
# RepeatMasker v4.2.2
# Infernal v1.1.5
# Augustus v3.5.0
# Miniprot v0.13
# HISAT2 v2.2.1
# StringTie v2.1.7
# TransDecoder v5.5.0
# EvidenceModeler
# eggNOG-mapper v2.1.13
############################################################


############################################################
# Repeat annotation
############################################################

BuildDatabase -name B59_db B59.reorder.fa

RepeatModeler -database B59_db -threads 48 -LTRStruct

RepeatMasker -e rmblast -pa 64 -gff -xsmall -s -a -lib B59-families.fa -dir repeatmasker_B59 B59.reorder.fa


############################################################
# ncRNA annotation
############################################################

cmscan --cpu 48 --cut_ga --rfam --nohmmonly --fmt 2 --tblout B59.tblout Rfam.cm B59.reorder.fa

perl infernal-tblout2gff.pl --cmscan --fmt2 B59.tblout > B59.ncRNA.gff3


############################################################
# Ab initio gene prediction
############################################################

autoAugTrain.pl --genome MorexV3.fasta --trainingset MorexV3.gff3 --species barley --optrounds 5 --cpus 64

bam2hints --in merged.bam --out B59.hints.gff

augustus --species=barley --extrinsicCfgFile=extrinsic.M.RM.E.W.cfg --hintsfile=B59.hints.gff --genemodel=complete --softmasking=1 --codingseq=on --exonnames=on

perl augustus_GTF_to_EVM_GFF3.pl augustus.gtf > B59.augustus.gff3


############################################################
# Homology-based annotation
############################################################

cat Morex.faa Aet.faa Rice.faa ZWY.faa Atlit.faa > total_pep.raw.faa

seqkit rmdup total_pep.raw.faa -o total_pep.rmdup.faa

cd-hit -i total_pep.rmdup.faa -o total_pep.clean.faa -c 0.95 -n 5

miniprot -t 48 -d B59.ref.mpi B59.final.fa.masked

miniprot -t 48 -I --gff --outc 0.8 --outn 1 B59.ref.mpi total_pep.clean.faa > B59.miniprot.gff3

python miniprot_GFF_2_EVM_GFF3.py B59.miniprot.gff3 > B59.pep.gff3


############################################################
# Transcript-based annotation
############################################################

hisat2-build B59.final.fa.masked B59

hisat2 -x B59 -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz --dta -S sample.sam

samtools sort -o sample.sort.bam sample.sam

samtools merge merged.bam *.sort.bam

stringtie -p 48 -o stringtie_merged.gtf merged.bam -G MorexV3.gff3

gtf_genome_to_cdna_fasta.pl stringtie_merged.gtf B59.final.fa.masked > transcripts.fasta

gtf_to_alignment_gff3.pl stringtie_merged.gtf > stringtie.gff3

TransDecoder.LongOrfs -t transcripts.fasta

diamond blastp -d uniprot_sprot.fasta -q longest_orfs.pep --max-target-seqs 1 --outfmt 6 > blastp.out

TransDecoder.Predict -t transcripts.fasta --retain_blastp_hits blastp.out

cdna_alignment_orf_to_genome_orf.pl transcripts.transdecoder.gff3 stringtie.gff3 transcripts.fasta > B59.transdecoder.gff3


############################################################
# EvidenceModeler
############################################################

partition_EVM_inputs.pl --genome B59.final.fa.masked --gene_predictions B59.augustus.gff3 --protein_alignments B59.pep.gff3 --segmentSize 100000 --overlapSize 10000 --partition_listing partitions.out

write_EVM_commands.pl --partitions partitions.out --genome B59.final.fa.masked --gene_predictions B59.augustus.gff3 --protein_alignments B59.pep.gff3 --output_file_name B59.evm.out --weights weights.txt

execute_EVM_commands.pl commands.list

recombine_EVM_partial_outputs.pl --partitions partitions.out --output_file_name B59.evm.out

convert_EVM_outputs_to_GFF3.pl --partitions partitions.out --output_file_name B59.evm.out --genome B59.final.fa.masked


############################################################
# Functional annotation
############################################################

emapper.py -i B59.protein.fa -m diamond --tax_scope 33090 -o B59

############################################################
# BUSCO assessment
############################################################

busco -i B59.protein.fa -l poales_odb12 -m proteins -o B59.busco -c 32 --offline

############################################################
# End
############################################################
