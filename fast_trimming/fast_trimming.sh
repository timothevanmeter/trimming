#!/bin/bash

# ARGUMENT PARSING
while [[ "$#" -gt 0 ]]; do
    case $1 in
	-i|--input) input=$2; shift ;;
	-h|--help) help=1 ;;
	-t|--threshold) threshold=$2; shift ;;
	-z|--summary) summary=1 ;;
        -s|--silent) silent=1 ;;
	-k|--keep) keep=1 ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

if [[ -z ${input} ]]; then
    printf "No input was provided!\nExiting program\n\n" > /dev/tty
    exit 1;
fi
printf "\ninput = %s\n" ${input} > /dev/tty

# printf "\nThreshold = %s\n" ${threshold} > ${stream}

if [[ -n ${help} && $help -eq 1 ]]; then 
    cat <<EOF
---------------------------------------------------------------------------
			- TRIMMING CODE -
This script processes a fastq file [ask tommaso what to put here]
There are N steps:
      -Extracting only the sequences between two adapters
      -Ensuring that these sequences have non-ambiguous nucleotides
      -Trimming them to 30 nucleotides
      -Discarding sequences with an expected error greater or equal to a threshold [threshold]
      -Translating nucleotide sequences into amino acids (without any stop codon)
      -Counting the frequency of each sequence in the dataset
---------------------------------------------------------------------------
Usage:
	./fast_trimming.sh -i <file_to_process.fastq> [options]
	./fast_trimming.sh -i <file_to_process.fastq> -o <output_file_name> [options]
Options:
	-i --input	fastq file to process.
	-o --output	name of the file in which the output should be saved. 
	-h --help	Show this information.
	-s --silent	Disable all messages.
	-k --keep	Keep all intermediary files generated during runtime.
	-t --threshold	Define the expected error threshold above which sequences are discarded. By default it is set to 0.01
---------------------------------------------------------------------------
For any questions or issues contact: timothe[dot]vm[at]posteo[dot]net
---------------------------------------------------------------------------
EOF
    exit 0;
fi
		   
if [[ -n ${silent} && $silent -eq 1 ]]; then 
    stream=/dev/null
else
    stream=/dev/tty
fi

####################################################################################
# Description of cutadapt options used in this script
# ---------------------------------------------------
# -cores X : [technical] simply the number of CPU cores to allocate for the parallelised task
# -e X : error rate tolerance
# "If E is an integer >= 1, then E errors in a full-length adapter match are allowed. For each specified adapter, this is converted to a maximum allowed error rate. This allows proportionally fewer errors for shorter (partial) adapter matches."
# -O X : setting the minimum allowed overlap in number of base
# -no-indels : Does not allow insertions and deletions when matching adapters against reads.
# -action=trim : Removes adapters when found.
# -untrimmed-output <filename> : Write the reads in which no adapter was found to a separate file.
# -g <sequence> : Specifying a 5' or linked adapter sequence.
# -o <filename> : Specify the file name for the output.
# --max-n X : 
##################################################

if [[ -n ${summary} && $summary -eq 1 ]]; then
    # # COUNTING THE TOTAL NUMBER OF READS IN THE FILE
    total=$(wc -l ${input} | awk '{print $1}')
    total=$(($total/4))
    # # WRITING THE INFORMATION TO SUMMARY FILE
    printf "native_sequence %d\n" $total > "${input}."summary.out
fi

# TRIMMING BOTH LINKED ADPATERS WITH NO INDELS AND ONLY IF THEY ARE EXACTLY MATCHED (length=10 && -O=10).
cutadapt --cores=4 -e 1 -O 10 --no-indels --action=trim --untrimmed-output untrimmed_.fastq -g GAGAGGCAAC...CGCCAAGCAG -o Fwd_.fastq ${input} > ${stream}

# TAKING THE READS FOR WHICH ADAPTERS SEQUENCES WERE NOT MATCHED AND TRYING TO MATCH THE COMPLEMENTARY SEQUENCES FOR BOTH ADAPTERS.
# UNTRIMMED SEQUENCES ARE DISCARDED.
cutadapt --cores=4 -e 1 -O 10 --no-indels --action=trim --discard-untrimmed -g CTGCTTGGCG...GTTGCCTCTC -o Rev_.fastq ./untrimmed_.fastq > ${stream}

# REVERSES AND THEN OBTAINS THE COMPLEMENTARY
#  SEQUENCE FROM THE INPUT SEQUENCE
# -r : inverses the sequence 5'===>3' ---> 3'===>5'
# -p : obtains the complementary sequence ATGCA ---> TACGT
# THIS WAY THE OUTPUTTED SEQUENCES CAN BE MERGED WITH THE
#  FORWARD ONES WITHOUT ANY PROBLEMS.
# seqkit seq -r -p Rev_.fastq -o Rev_RC.fastq
# Same command, but removes the warning:
seqkit --seq-type dna seq -r -p Rev_.fastq -o Rev_RC.fastq

# MERGES EVERYTHING INTO A SINGLE FILE
cat Fwd_.fastq Rev_RC.fastq > data.fastq

if [[ -n ${summary} && $summary -eq 1 ]]; then
    # # COUNTING NUMBER OF TRIMMED SEQUENCES
    trimmed=$(wc -l data.fastq | awk '{print $1}')
    trimmed=$(($trimmed/4))
    # # WRITING THE INFORMATION TO SUMMARY FILE
    printf "trimmed_sequence %d\n" $trimmed >> "${input}."summary.out
fi

# INSURES TO KEEP ONLY SEQUENCES OF EXACTLY 30 BASES (max=min=30)
cutadapt --cores=4 --minimum-length 30 -o noless30.fastq data.fastq > ${stream}
cutadapt --cores=4 --maximum-length 30 -o data30.fastq noless30.fastq > ${stream}

if [[ -n ${summary} && $summary -eq 1 ]]; then
    # # COUNTING NUMBER OF TRIMMED SEQUENCES
    trente=$(wc -l data30.fastq | awk '{print $1}')
    trente=$(($trente/4))
    # # WRITING THE INFORMATION TO SUMMARY FILE
    printf "30bp_sequence %d\n" $trente >> "${input}."summary.out
fi

# REMOVES ALL THE 'N' NUCLEOTIDE, I.E. THE ONES THAT HAVE AMBIGUOUS
#  IDENTITIES IN THE SEQUENCE
cutadapt --cores=4 --max-n 0 -o data10adapt.fastq data30.fastq > ${stream}


# ------------------------------------------------------------------------------


# OUTPUTTING THE EXPECTED ERROR FOR EACH READ :: DEBUGGING
# vsearch --fastq_filter data_10adapt.fastq --fastq_ascii 33 --eeout --fastaout data_10adapt_errors.fasta > ${stream}

# QUALITY FILTER
# DISCARDING ALL THE READS WITH AN EXPECTED ERROR GREATER OR EQUAL TO 0.01
vsearch --fastq_filter data10adapt.fastq --eeout --fastq_ascii 33 --fastq_maxee $threshold --fastaout data10adapt_filtered.fasta > ${stream}

if [[ -n ${summary} && $summary -eq 1 ]]; then
    quality=$(wc -l data10adapt_filtered.fasta | awk '{print $1}')
    quality=$(($quality/2))
    printf "Qthreshold %d\n" $quality >> "${input}."summary.out
fi

# TRANSLATING THE READS FROM DNA TO AMINO ACIDS
seqkit translate --allow-unknown-codon -f 1 data10adapt_filtered.fasta > AA_filtered.fasta

# REMOVING THE HEADERS IN THE FASTA FILE
# awk '!/>/{print $0}' AA_filtered.fasta > AA_filtered.seq

# DISCARDING ALL THE READS CONTAINING A STOP CODON
awk '{if(NR%2==0){if($0!~/\*/){print header "\n" $0}}else{header=$0}}' AA_filtered.fasta > AA_filtered-nostop.fasta

if [[ -n ${summary} && $summary -eq 1 ]]; then
    dixmer=$(wc -l data10adapt_filtered.fasta | awk '{print $1}')
    dixmer=$(($dixmer/2))
    printf "Qthreshold %d\n" $dixmer >> "${input}."summary.out
fi

# TRIMMING THE PROTEIN SEQUENCE
# REMOVING THE LINKER AMINO ACIDS ON EACH SIDE OF THE SEQUENCE.
awk '{if(NR%2==0){print substr($0,3,7)}else{print $0}}' AA_filtered-nostop.fasta > AA_filtered-nostop-cropped.fasta

# COLLAPSING READS:
# COUNTING THE FREQUENCY OF EACH SEQUENCE AND SAVING THE COUNTS
./colps AA_filtered-nostop-cropped.fasta AA_filtered-nostop-cropped-collapsed.fasta

if [[ -n ${summary} && $summary -eq 1 ]]; then
    septmer=$(awk 'NR>1{sum+=$2}END{print sum}' AA_filtered-nostop-cropped-collapsed.fasta)
    printf "7mer %d\n" $septmer >> "${input}."summary.out
fi

# Graphing the results 
# gnuplot -p plot.scores.gp ; display time_scores.hist.png

if [[ -n ${keep} && $keep -eq 1 ]]; then 
    printf "\nAll intermediary files were kept in the local folder" > ${stream}
else
    rm untrimmed_.fastq 
    rm Fwd_.fastq
    rm Rev_.fastq
    rm Rev_RC.fastq
    rm noless30.fastq
    rm data.fastq
    rm data30.fastq
    rm data10adapt.fastq
    rm data10adapt_filtered.fasta
    rm AA_filtered.fasta
    rm AA_filtered-nostop.fasta
    rm AA_filtered-nostop-cropped.fasta
    printf "\nAll intermediary files were removed" > ${stream}
fi

printf "#----------------------------------------------------#" > ${stream}
printf "" > ${stream}
printf "           DONE!" > ${stream}
printf "" > ${stream}


