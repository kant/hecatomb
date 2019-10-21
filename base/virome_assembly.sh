#!/bin/bash

# Takes quality controlled reads (though contaminant_remval.sh Step 7) and performe a metagenomic assembly

# Steps:
# 0) Some stats file generation
# 1) Digital normalization
# 2) Assembly
# 3) Assembly decontamination
# ?) Sketch comparison
# ?) Virfinder
# ?) Bracken


# Prep output directories
mkdir -p ./assembly
mkdir -p ./assembly/stats
mkdir -p ./assembly/megahit_contigs

# Set variables
IN=./QC/step_7
OUT=./assembly

# Set file names
for i in $IN/*_R1.s7.out.fastq; do
        F=`basename $i _R1.s7.out.fastq`;

# Step 0: Tabulate some kmer statistics
bbcountunique.sh in=$IN/"$F"_R1.s7.out.fastq in2=$IN/"$F"_R2.s7.out.fastq \
	interval=2500 \
	out=$OUT/stats/"$F"_uniq_kmer_stats.txt;

# Step 1: Digital Normalization

bbnorm.sh in=$IN/"$F"_R1.s7.out.fastq in2=$IN/"$F"_R2.s7.out.fastq extra=$IN/"$F"_singletons.s7.out.fastq \
	out=$OUT/"$F"_R1.norm.out.fastq out2=$OUT/"$F"_R2.norm.out.fastq outt=$OUT/"$F"_tossed.norm.fastq \
	fixspikes=t \
	target=40 mindepth=2 \
	hist=$OUT/"$F"_norm.hist \
	ow=t;

# Step 2: Assembly
megahit -1 $OUT/"$F"_R1.norm.out.fastq -2 $OUT/"$F"_R2.norm.out.fastq -o $OUT/"$F"_megahit_out \
	--out-prefix "$F".mh;
cp $OUT/"$F"_megahit_out/*.mh.contigs.fa $OUT/megahit_contigs;

done