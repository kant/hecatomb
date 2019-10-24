#!/bin/bash bash

# Set bash options for proper failure recognition
set -o errext
set -o pipefail

### DESCRIPTION ###
# Script to remove non-biological sequences (primers, adapters), low-quality bases, host-sequences and obiovus (100% ID)  bacterial sequences from virome sequenced libraries
# References:
# Heavy reliance on:
        # BBtools: https://jgi.doe.gov/data-and-tools/bbtools/

# Set Variables
CONPATH=/mnt/data1/databases/contaminants # adapter/primers
HOSTPATH=/mnt/data1/databases/human_masked # host
BACPATH=/mnt/data1/databases/masked_dbs/genome-collection-combos/reformated/bac_giant # masked bacterial and giant virus genome

#re-establish loop to finish everything
for i in *_R1.fastq.gz; do
        F=`basename $i _R1.fastq.gz`;

# Step 7: Trim low-quality bases
bbduk.sh in=./QC/step_6/"$F"_unmapped.s6.out.fastq \
	out=./QC/step_7/"$F"_R1.s7.out.fastq out2=./QC/step_7/"$F"_R2.s7.out.fastq outs=./QC/step_7/"$F"_singletons.s7.out.fastq \
	stats=./QC/step_7/"$F".s7.stats \
	qtrim=r trimq=30 \
	maxns=2 minlength=50 \
	ordered=t \
	ow=t;

# Step 8: Merge forward (R1) and reverse (R2) reads
bbmerge.sh in1=./QC/step_7/"$F"_R1.s7.out.fastq in2=./QC/step_7/"$F"_R2.s7.out.fastq \
	out=./QC/step_8/"$F"_merged.fastq outu1=./QC/step_8/"$F"_R1.unmerged.fastq outu2=./QC/step_8/"$F"_R2.unmerged.fastq \
	rem k=62 extend2=50 ecct vstrict=t \
	ordered=t \
	-Xmx128g \
	ow=t;

	# Split singletons and combine R1 and R2 files
	grep -A 3 '1:N:' ./QC/step_7/"$F"_singletons.s7.out.fastq | sed '/^--$/d' > ./QC/step_7/"$F"_singletons_R1.out.fastq;
	grep -A	3 '2:N:' ./QC/step_7/"$F"_singletons.s7.out.fastq | sed '/^--$/d' > ./QC/step_7/"$F"_singletons_R2.out.fastq;

	cat ./QC/step_8/"$F"_merged.fastq ./QC/step_8/"$F"_R1.unmerged.fastq ./QC/step_7/"$F"_singletons_R1.out.fastq > ./QC/step_8/"$F"_R1.s8.out.fastq;
	cat ./QC/step_8/"$F"_merged.fastq ./QC/step_8/"$F"_R2.unmerged.fastq ./QC/step_7/"$F"_singletons_R2.out.fastq > ./QC/step_8/"$F"_R2.s8.out.fastq;

# Step 9: Remove bacterial contaminants reserving viral and ambiguous sequences
bbmap.sh in=./QC/step_8/"$F"_R1.s8.out.fastq \
	path=$BACPATH \
	outm=./QC/step_9/"$F"_bacterial.fastq outu=./QC/step_9/"$F"_viral_amb.fastq \
	semiperfectmode=t \
	quickmatch fast \
	refstats=./QC/step_9/"$F"_refstats.txt scafstats=./QC/step_9/"$F"_scafstats.txt \
	ordered=t \
	ow=t;

done
