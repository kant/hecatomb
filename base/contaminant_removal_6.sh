#!/bin/bash bash

# Set bash options for proper failure recognition
#set -o errext
#set -o pipefail

### DESCRIPTION ###
# Script to remove non-biological sequences (primers, adapters), low-quality bases, host-sequences and obiovus (100% ID)  bacterial sequences from virome sequenced libraries
# References:
# Heavy reliance on:
        # BBtools: https://jgi.doe.gov/data-and-tools/bbtools/

##function for joining array, thanks stack!
function join_by { local IFS="$1"; shift; echo "$*"; }

# Set Variables
CONPATH=/mnt/data1/databases/contaminants # adapter/primers
HOSTPATH=/mnt/data1/databases/human_masked # host
BACPATH=/mnt/data1/databases/masked_dbs/genome-collection-combos/reformated/bac_giant # masked bacterial and giant virus genome

# Step 6: Host removal

#gather all file names and input to bbwrap
COUNT=$(ls ./QC/step_5/*R1* | wc -l)
TEMP=0
declare -a NAMES 
declare -a NAMES_R1
declare -a NAMES_R2
declare -a MAPPED
declare -a UNMAPPED

for i in ./QC/step_5/*_R1.s5.out.fastq; do
  NAMES[${TEMP}]=`basename $i _R1.s5.out.fastq`;
  NAMES_R1[${TEMP}]=./QC/step_5/${NAMES[$TEMP]}_R1.s5.out.fastq;
  NAMES_R2[${TEMP}]=./QC/step_5/${NAMES[$TEMP]}_R2.s5.out.fastq;
  MAPPED[${TEMP}]=./QC/step_6/${NAMES[$TEMP]}_hostmapped.s6.out.fastq;
  UNMAPPED[${TEMP}]=./QC/step_6/${NAMES[$TEMP]}_unmapped.s6.out.fastq
  TEMP=$((TEMP+1));
done

bbwrap.sh in1=$(join_by , ${NAMES_R1[*]}) in2=$(join_by , ${NAMES_R2[*]}) \
	outm=$(join_by , ${MAPPED[*]}) outu=$(join_by , ${UNMAPPED[*]}) \
        semiperfectmode=t \
        quickmatch fast \
        ordered=t \
        path=$HOSTPATH \
        ow=t
      

#old single mapping as model for wrap
#bbmap.sh in=./QC/step_5/"$F"_R1.s5.out.fastq in2=./QC/step_5/"$F"_R2.s5.out.fastq \
#	outu=./QC/step_6/"$F"_unmapped.s6.out.fastq outm=./QC/step_6/"$F"_hostmapped.s6.out.fastq \
#	semiperfectmode=t \
#	quickmatch fast \
#	ordered=t \
#	path=$HOSTPATH \
#	ow=t;
