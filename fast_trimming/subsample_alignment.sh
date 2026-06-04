#!/bin/bash



subsample_size=1000

# ---------------------------------
# Subsampling the fastq file
subsampling_seq.o $1 ${1/.fastq/_sub}${subsample_size}".fastq" $subsample_size

# ---------------------------------
# Convert to fasta
fastq_to_fasta.awk ${1/.fastq/_sub}${subsample_size}".fastq" > ${1/.fastq/_sub}${subsample_size}".fa"

# ---------------------------------
# Multiple sequences alignment
# Gap opening penalty
gp_op=2.0
# Gap extension penalty
gp_ex=1.0

# Create a MSA of the subsampled sequences
mafft --thread 24 --maxiterate 100 --op $gp_op --ep $gp_ex ${1/.fastq/_sub}${subsample_size}".fa" > ${1/.fastq/_sub}${subsample_size}".align"



											  # # List of lines to reverse for the alignment
# list=(6 7 9 13 14 19 21 24 25 27 30 31 32 39 40 44 45 48)
# index=1
# ilist=0
# while read -r line; do
#     temp=$(($index % 2))
#     if [[ $temp -eq 0 ]]; then
# 	# echo "seq #"$(($index / 2))
# 	if [[ $(($index / 2)) -eq ${list[$ilist]} ]]; then
# 	    # echo "(REV)"
# 	    seq=$(echo $line | rev)
# 	    comp_dna.sh $seq
# 	    echo ""
# 	    # echo $line
# 	    ((ilist++))
# 	else
# 	    echo $line
# 	fi
#     else
# 	echo $line
#     fi
#     # echo ${list[$ilist]}
#     ((index++))
# done < $1
