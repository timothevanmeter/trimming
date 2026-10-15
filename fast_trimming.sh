#!/bin/bash


# Parallel usage:
# parallel --eta -j 4 fast_trimming.sh -z -t 0.1 -is 21 -a1 <adapter1> -a2 <adapter2> -i {} -jn {#} ::: *.fastq.gz


# ARGUMENT PARSING
while [[ "$#" -gt 0 ]]; do
    case $1 in
	-i|--input) input=$2; shift ;;	
	-t|--threshold) threshold=$2; shift ;;
	-a1|--adapter1) adapter1=$2; shift ;;
	-a2|--adapter2) adapter2=$2; shift ;;
	-is|--insertSize) insertSize=$2; shift ;;
	-jn|--jobNumber) jobNumber=$2; shift ;;
	-h|--help) help=1 ;;
	-z|--summary) summary=1 ;;
        -s|--silent) silent=1 ;;
	-k|--keep) keep=1 ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

# ------------------------------------------------
# Pushing all stdout and stderr to the log file:
# Strips .gz then .fastq so plain and compressed inputs give the same name
LOGFILE=${input%.gz}
LOGFILE=${LOGFILE%.fastq}"-trimming.log"

# Redirect stdout and stderr
exec > >(tee -a "$LOGFILE") \
     2> >(tee -a "$LOGFILE" >&2)
# ------------------------------------------------

# CHECKING THAT THE MANDATORY ARGUMENTS WERE PROVIDED
if [[ -z ${input} ]]; then
    printf "No input was provided!\nExiting program\n\n" > /dev/tty
    exit 1;
fi
if [[ -z ${adapter1} ]]; then
    printf "No sequence was provided for adapter 1!\nExiting program\n\n" > /dev/tty
    exit 1;
fi
if [[ -z ${adapter2} ]]; then
    printf "No sequence was provided for adapter 2!\nExiting program\n\n" > /dev/tty
    exit 1;
fi
if [[ -z ${insertSize} ]]; then
    printf "No insert size was provided!\nExiting program\n\n" > /dev/tty
    exit 1;
fi

# CHECKING IF THE INPUT IS COMPRESSED OR NOT
if [[ ${input} == *".gz"* ]]; then
    output=${input%.fastq.gz}
    compressed=1
else
    output=${input%.fastq}
    compressed=0
fi

printf "\n\tinput = %s\n\tadapter 5 = %s\n\tadapter 3 = %s\n\tinsert size = %d\n" ${input} ${adapter1} ${adapter2} ${insertSize} > /dev/tty

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
	-k --keep	Keep all intermediary files generated during runtime.
	-t --threshold	Define the expected error threshold above which sequences are discarded. By default it is set to 0.01
	-a1 --adapter1	Define the 5` sequence around the sequence of interest 
	-a2 --adapter2	Define the 3` sequence around the sequence of interest 
	-is --insertSize       Define the size of the insert or sequence of interest that you expect: only sequences of this exact size will be kept
---------------------------------------------------------------------------
For any questions or issues contact: timothe[dot]vm[at]posteo[dot]net
---------------------------------------------------------------------------
EOF
    exit 0;
fi
		   
# if [[ -n ${silent} && $silent -eq 1 ]]; then 
#     stream=/dev/null
# else
#     stream=/dev/tty
# fi

############################################################
# EXPORTING all the processing information to a corresponding
# json file to document how the trimming was performed
############################################################
# --------------------------------------
# INITIATING the variables and their values to be passed
# to jq:
DATE=$(date -I)
FOLDER=$(pwd)
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
SUMMARY="${output}.summary"
COUNTS="${output}.counts"
# --------------------------------------
# CALLING jq to create the json file:
jq -n \
   --arg script "$0" \
   --arg jobNumber "$jobNumber" \
   --arg adapter1 "$adapter1" \
   --arg adapter2 "$adapter2" \
   --arg path "$FOLDER" \
   --arg scriptPath "$SCRIPT_DIR" \
   --arg input "$input" \
   --arg summary "$SUMMARY" \
   --arg logFile "$LOGFILE" \
   --arg counts "$COUNTS" \
   --arg date "$DATE" \
   --arg name "vsearch" \
   --arg threshold "$threshold" \
   --arg insertSize "$insertSize" \
   '{ main: { script: $script, jobNumber: $jobNumber, path: $scriptPath, input: $input, path: $path, adapter1: $adapter1, adapter2: $adapter2, insertSize: $insertSize, threshold: $threshold, summary: $summary, counts: $counts, logFile: $logFile, date: $date }}' > ${output}"-trimming.json"
############################################################

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
    if [[ $compressed -eq 1 ]]; then
	total=$(zcat ${input} | wc -l | awk '{print $1}')
    else
	total=$(cat ${input} | wc -l | awk '{print $1}')
    fi
    total=$(($total/4))
    # # WRITING THE INFORMATION TO SUMMARY FILE
    printf "native_sequence %d\n" $total > "${output}.summary"
fi

# TRIMMING BOTH LINKED ADPATERS WITH NO INDELS AND ONLY IF THEY ARE EXACTLY MATCHED (length=10 && -O=10).
# cutadapt --cores=24 --cores=24 -e 0.2 -O 10 --no-indels --action=trim --untrimmed-output untrimmed_.fastq -g GAGCACCGGT...GGCTTAAGTT -o Fwd_.fastq ${input}
# Try sequences further away from the 7mer
cutadapt --cores=24 -e 0.2 -O 10 --no-indels --action=trim --untrimmed-output untrimmed_${jobNumber}.fastq -g ${adapter1}...${adapter2} -o Fwd_${jobNumber}.fastq ${input}

rc_a1=$(comp_dna.sh ${adapter1} | rev)
rc_a2=$(comp_dna.sh ${adapter2} | rev)

# TAKING THE READS FOR WHICH ADAPTERS SEQUENCES WERE NOT MATCHED AND TRYING TO MATCH THE COMPLEMENTARY SEQUENCES FOR BOTH ADAPTERS.
# UNTRIMMED SEQUENCES ARE DISCARDED.
cutadapt --cores=24 -e 0.2 -O 10 --no-indels --action=trim --discard-untrimmed -g ${rc_a2}...${rc_a1} -o Rev_${jobNumber}.fastq ./untrimmed_${jobNumber}.fastq

# REVERSES AND THEN OBTAINS THE COMPLEMENTARY
#  SEQUENCE FROM THE INPUT SEQUENCE
# THIS WAY THE OUTPUTTED SEQUENCES CAN BE MERGED WITH THE
#  FORWARD ONES WITHOUT ANY PROBLEMS.
vsearch --fastx_revcomp Rev_${jobNumber}.fastq -fastqout Rev_RC_${jobNumber}.fastq

# MERGES EVERYTHING INTO A SINGLE FILE
cat Fwd_${jobNumber}.fastq Rev_RC_${jobNumber}.fastq > data_${jobNumber}.fastq

if [[ -n ${summary} && $summary -eq 1 ]]; then
    # # COUNTING NUMBER OF TRIMMED SEQUENCES
    trimmed=$(wc -l data_${jobNumber}.fastq | awk '{print $1}')
    trimmed=$(($trimmed/4))
    # # WRITING THE INFORMATION TO SUMMARY FILE
    printf "trimmed_sequence %d\n" $trimmed >> "${output}.summary"
fi

# INSURES TO KEEP ONLY SEQUENCES OF EXACTLY THE NUMBER OF BASES
# GIVEN BY THE VARIABLE insertSize (max=min=insertSize)
cutadapt --cores=24 --minimum-length ${insertSize} -o minSize_${jobNumber}.fastq data_${jobNumber}.fastq
cutadapt --cores=24 --maximum-length ${insertSize} -o data_insertSize_${jobNumber}.fastq minSize_${jobNumber}.fastq

if [[ -n ${summary} && $summary -eq 1 ]]; then
    # # COUNTING NUMBER OF TRIMMED SEQUENCES
    size=$(wc -l data_insertSize_${jobNumber}.fastq | awk '{print $1}')
    size=$(($size/4))
    # # WRITING THE INFORMATION TO SUMMARY FILE
    printf "%dbp_sequence %d\n" $insertSize $size >> "${output}.summary"
fi

# REMOVES ALL THE SEQUENCES CONTAINING ONE OR MORE 'N' NUCLEOTIDE,
# I.E. THE ONES THAT HAVE AMBIGUOUS IDENTITIES IN THE SEQUENCE
cutadapt --cores=24 --max-n 0 -o data_no-ambiguous-na_${jobNumber}.fastq data_insertSize_${jobNumber}.fastq

# ------------------------------------------------------------------------------
# QUALITY FILTER
# DISCARDING ALL THE READS WITH AN EXPECTED ERROR GREATER OR EQUAL TO 0.1
vsearch --fastq_filter data_no-ambiguous-na_${jobNumber}.fastq --eeout --fastq_ascii 33 --fastq_maxee $threshold --fastaout data_filtered_${jobNumber}.fasta

if [[ -n ${summary} && $summary -eq 1 ]]; then
    quality=$(wc -l data_filtered_${jobNumber}.fasta | awk '{print $1}')
    quality=$(($quality/2))
    printf "Qthreshold %d\n" $quality >> "${output}.summary"
fi

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

# TRANSLATING THE READS FROM DNA TO AMINO ACIDS
# seqkit translate --allow-unknown-codon -f 1 data_filtered_${jobNumber}.fasta > AA_filtered_${jobNumber}.fasta
translate data_filtered_${jobNumber}.fasta AA_filtered_${jobNumber}.fasta

# DISCARDING ALL THE READS CONTAINING A STOP CODON
awk '{if(NR%2==0){if($0!~/\*/){print header "\n" $0}}else{header=$0}}' AA_filtered_${jobNumber}.fasta > AA_filtered-nostop_${jobNumber}.fasta

if [[ -n ${summary} && $summary -eq 1 ]]; then
    dixmer=$(wc -l AA_filtered-nostop_${jobNumber}.fasta | awk '{print $1}')
    dixmer=$(($dixmer/2))
    printf "stop_codon_filter %d\n" $dixmer >> "${output}.summary"
fi

# COLLAPSING READS:
# COUNTING THE FREQUENCY OF EACH SEQUENCE AND SAVING THE COUNTS
hash_collapser AA_filtered-nostop_${jobNumber}.fasta "${output}.counts"

if [[ -n ${summary} && $summary -eq 1 ]]; then
    septmer=$(awk '{sum+=$2}END{print sum}' "${output}.counts")
    printf "7mer %d\n" $septmer >> "${output}.summary"
fi


# REMOVING OR KEEPING THE INTERMEDIARY FILES
if [[ -n ${keep} && $keep -eq 1 ]]; then 
    printf "\nAll intermediary files were kept in the local folder"
else
    rm untrimmed_${jobNumber}.fastq 
    rm Fwd_${jobNumber}.fastq
    rm Rev_${jobNumber}.fastq
    rm Rev_RC_${jobNumber}.fastq
    rm data_${jobNumber}.fastq
    rm data_insertSize_${jobNumber}.fastq
    rm minSize_${jobNumber}.fastq
    rm data_no-ambiguous-na_${jobNumber}.fastq
    rm data_filtered_${jobNumber}.fasta
    rm AA_filtered_${jobNumber}.fasta
    rm AA_filtered-nostop_${jobNumber}.fasta
    printf "\nAll intermediary files were removed"
fi

printf "#----------------------------------------------------#\n"
printf "\n"
printf "\t\tDONE!\n"
printf "\n"
