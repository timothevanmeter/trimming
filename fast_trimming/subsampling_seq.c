#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "kseq.h"
#include "getlines.h"

//--------------------------------------------------------
// USAGE:

// COMPILATION:
// gcc -Wall -pedantic -o subsampling_seq.o subsampling_seq.c getlines.c ht.c -lm

// RUNNING THE EXECUTABLE:
// ./subsampling.o <input_file> <output_file> <probability_of_sampling> OR <number_of_sequences>

// EXAMPLE:
// ./subsampling.o qscore.dat qscore_out.dat 0.1
// OR 
// ./subsampling.o qscore.dat qscore_out.dat 1000
//--------------------------------------------------------

KSEQ_INIT(int, read)

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

void printProgress(double percentage) {
    int val = (int) (percentage * 100);
    int lpad = (int) (percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    fprintf(stdout, "\t\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    fflush(stdout);
}

// -----------------------------------------------------------
int main( int argc, char **argv )
{
  srand(time(NULL));   // Initialization, should only be called once.
  kseq_t *seq;
  FILE * input;
  FILE * out;
  // char ch;
  char *end;
  // char BUFFER[30];
  unsigned int count = 0;
  unsigned int kept = 0;
  long unsigned int lines = getlines(argv[1]) / 4;
  long unsigned int number = 0;
  double probability;
  // printf("\n\tNumber of lines: %ld\n", lines);
  // Checking whether the argument is a probability or a
  // number of sequences to be kept.
  if(strtof(argv[3], &end) > 1.0) {
    number = strtol(argv[3], &end, 10);
    // printf("\n\tNumber of lines to keep: %ld", number);
    probability = number / lines  + 0.01;
  } else {
    probability = strtof(argv[3], &end);
    // printf("\n\tProbability to keep a sequence: %.4f", probability);
  }
  if(number > lines) {
    // Error message
    printf("\nERROR: Argument %s is bigger than the number of sequences in input file %s\nPlease choose a number of sequences to keep smaller or equal to %ld\n", argv[3], argv[1], lines);
    exit(1);
  }  
  double r = 0.0;
  // double fileSize = (double)getlines(argv[1]);
  // double i = 0.0;
    
  fprintf(stdout,"-----------------------------------------------------\n");
  input = fopen(argv[1], "r");
  seq = kseq_init(fileno(input));
  
  fprintf(stdout, "\n\tprocessing file %s\n", argv[1]);
  fflush(stdout);
  out = fopen(argv[2], "w+");
  // --------------------------------------

    while (kseq_read(seq) >= 0) {
      r = (double)rand() / RAND_MAX;
      if(probability >= r) {
	// fprintf(stdout, "%s\n", BUFFER) ;
	// memcpy( temp, seq->seq.s, MAX_SEQ_SIZE );
	fprintf(out, "@%s\n", seq->name.s);
	fprintf(out, "%s\n", seq->seq.s);
	fprintf(out, "+\n%s\n", seq->qual.s) ;
	// fprintf(stdout, "%s\n", seq->seq.s) ;
	kept++;
      }
      count++;
      if(number != 0 && number == kept) {
	break;
      }
      // i++;
      // printProgress(((double)count + 1) / fileSize);
      // printProgress((i + 1) / fileSize);
  }
  fflush(out);

  kseq_destroy(seq);
  fclose(input);
  fclose(out);
  fprintf(stdout,"\t%d lines processed and %d lines kept\n", count, kept);
  fprintf(stdout,"\n-----------------------------------------------------\n");
  return 0;
}
