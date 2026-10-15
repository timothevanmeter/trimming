#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <unistd.h>
#include <limits.h>
#include "ht_array_DYNAMIC.h"

// ##############################################################
// This script expects as input a space-separated file with
// in the first column a key used to identify an entry and
// in the second column an integer value.

// The files listed as input are read and the keys are used to
// populate a hash table. The hash table has any key associated
// with a vector of all the integer values corresponding to that
// key.

// The script resolves the problem of presence / absence of keys
// between different files by including 0 values in the vectors.

// -----------------
// USAGE:

// How to call the executable:
// $ ./conc.o [-o <output_name>] <file_1_to_concatenate> <file_2_to_concatenate> <file_N_to_concatenate>
// The ".csv" extension is appended to <output_name>.
// Without -o the output is saved under concatenation_output.csv

// The makefile has 5 different options
// By default if make is launch without any options it will do make all (option 5)

// 1) COMPILATION and production of an executable:
// $ make compile

// 2) RUNNING the executable:
//     By default this is done on the files T4_retina_100.counts and T6_retina_100.counts located in the data_test folder
// $ make run

// 3) TESTING the executable:
//     By default this is done on the files _m1.counts and _m1.counts located in the data_test folder
// $ make test

// 4) REMOVING the output file and the executable:
// $ make clean

// 5) COMPILING and RUNNING options 1) + 2)
// $ make all
// -----------------

// ##############################################################
#define KEY_SIZE 40
#define DEFAULT_OUTPUT "concatenation_output"
int TOTAL_FILE_NUMBER = 0;

// -----------------------------------------------------------
// -----------------------------------------------------------

void usage(const char *prog) {
  fprintf(stderr, "Usage: %s [-o output_name] file_1 file_2 [file_N ...]\n", prog);
}

// -----------------------------------------------------------
// -----------------------------------------------------------

void exit_nomem(void) {
  fprintf(stderr, "out of memory\n");
  exit(1);
}
// -------------------------------------------------
// -------------------------------------------------
void print_hash(ht *htable) {
  hti iterator = ht_iterator(htable);
  hti *it = &iterator;
  while(ht_next(it)) {
    fprintf(stdout, "%s,", it->key);
    int i = 0;
    while(i < TOTAL_FILE_NUMBER) {
      fprintf(stdout, "%d,", it->counts[i]);
      i++;
    }
    fprintf(stdout, "\n");
  }
}

// -------------------------------------------------
// -------------------------------------------------
void hash_to_file(ht *htable, const char *output, char **files) {
  printf("\nSaving results to file: %s\n", output);
  FILE *file;
  hti iterator = ht_iterator(htable);
  hti *it = &iterator;
  file = fopen(output, "w");
  // WRITING the header to the file
  fprintf(file, "identifier,");
  for(int i=0; i < TOTAL_FILE_NUMBER; i++) {
    fprintf(file, "%s", files[i]);
    if(i != TOTAL_FILE_NUMBER - 1) {
      fprintf(file, ",");
    }
  }
  fprintf(file, "\n");
  // WRITING the data to the file
  while(ht_next(it)) {
    fprintf(file, "%s,", it->key);
    int i = 0;
    while(i < TOTAL_FILE_NUMBER) {
      fprintf(file, "%d", it->counts[i]);
      i++;
      if(i != TOTAL_FILE_NUMBER) {
	fprintf(file, ",");
      }
    }
    fprintf(file, "\n");
  }
  fclose(file);
}

// -------------------------------------------------
void test_input_format(int n_files, char **files) {
  fprintf( stdout, "---------------------------\n");
  char ch;
  int i = 0;
  int READ = 1;
  FILE *fp;
  for(int f = 0; f < n_files; f++) {
    fp = fopen(files[f], "r");
    fprintf(stdout, "Verifying format of file: %s\n", files[f]);
    fflush(stdout);
    ch = fgetc(fp);
    READ = 1;
    i = 0;
    while(ch != '\n' && ch != '\r') {
      if(ch == ' ') {
	READ++;
	ch = fgetc(fp);
	continue;
      }
      if(READ == 1) {
	if(i > KEY_SIZE) {
	  fprintf(stderr, "\nERROR: first field's identifier exceeds KEY_SIZE's length=%d\nEdit the length of KEY_SIZE in concatenate_tables.c:48 to allow for longer keys.", KEY_SIZE);
	  exit(EXIT_FAILURE);
	}
      } else {
	if(!isdigit(ch)) {
	  fprintf(stderr, "\nERROR: second field's count has non-digit characters: %c\n", ch);
	  exit(EXIT_FAILURE);
	}
      }
      ch = fgetc(fp);
      i++;
    }
    if(READ == 1) {
      fprintf(stderr, "\nERROR: the fields are not space separated\n");
      exit(EXIT_FAILURE);
    }
    fclose(fp);
  }
}

// ###########################################################################

int main(int argc, char **argv)
{
  FILE * fp;
  const char * out_base = DEFAULT_OUTPUT;
  char output[PATH_MAX];
  int opt;
  int output_given = 0;
  char * key = malloc(sizeof(char) * KEY_SIZE);
  char * countC = malloc(sizeof(char) * KEY_SIZE);
  int countI;
  int n = 0;
  
  
  fprintf(stdout, "-----------------------------------------------------");
  // NEED TO:
  // TESTING IF THE INPUT FILE IS CORRECTLY FORMATTED !

  // TESTING that the size of the key is not larger
  // than the allocated string size
  // Testing that they are only two fields per line
  // separated by whitespace
  // test_input_format(n_files, files);

  // PARSING THE OPTIONS: -o <output_name>
  while((opt = getopt(argc, argv, "o:h")) != -1) {
    switch(opt) {
    case 'o':
      out_base = optarg;
      output_given = 1;
      break;
    case 'h':
      usage(argv[0]);
      exit(EXIT_SUCCESS);
    default:
      usage(argv[0]);
      exit(EXIT_FAILURE);
    }
  }
  // Every remaining argument is a file to concatenate
  char ** files = argv + optind;
  TOTAL_FILE_NUMBER = argc - optind;

  // TESTING THAT THE MINIMUM ARGUMENTS WERE PROVIDED
  // FOR THE SCRIPT TO FUNCTION
  if(TOTAL_FILE_NUMBER < 2) {
    fprintf(stderr, "\nERROR: Not enough arguments!\nNeed at least two files to concatenate\n");
    usage(argv[0]);
    exit(EXIT_FAILURE);
  }
  if(out_base[0] == '\0') {
    fprintf(stderr, "\nERROR: the output name is empty\n");
    exit(EXIT_FAILURE);
  }
  int len = snprintf(output, sizeof output, "%s.csv", out_base);
  if(len < 0 || (size_t)len >= sizeof output) {
    fprintf(stderr, "\nERROR: the output name is too long\n");
    exit(EXIT_FAILURE);
  }
  if(!output_given) {
    fprintf(stdout, "\nWARNING: No -o option was provided for the output file's name\nBy default it will be saved under %s\n", output);
  }
  test_input_format(TOTAL_FILE_NUMBER, files);

  fprintf( stdout, "---------------------------\n");
  printf("FILE NUMBER = %d\n", TOTAL_FILE_NUMBER);
  fprintf( stdout, "---------------------------\n");
  ht * htable = ht_create();
  if (htable == NULL) {
    exit_nomem();
  }

  // LOOPING OVER ALL THE INPUT FILES CONTAINING THE COUNTS
  for(int f = 0; f < TOTAL_FILE_NUMBER; f++) {
    fp = fopen(files[f], "r");
    fprintf(stdout, "\n\tProcessing file %s\n", files[f]);
    fflush(stdout);
    while (fscanf(fp, "%39s %d", key, &countI) == 2) {
      // ht_set expects a 1-based file number
      ht_set(htable, key, &countI, f + 1);
      n++;
    }
    fclose(fp);
    fprintf(stdout,"\t%d items processed\n\n", n);
    n = 0;
  }
  fprintf( stdout, "---------------------------\n");
  // print_hash(htable);
  hash_to_file(htable, output, files);
  // CLEANING UP ALL THE MESS,
  //  CLOSING DOORS AND TURNING OFF THE LIGHTS!
  ht_destroy(htable);
  free(countC);
  free(key);
  
  printf("-----------------------------------------------------\n");
  return 0;
}
