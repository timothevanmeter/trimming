#include <zlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "kseq.h"
#include "ht.h"
#include "list_counts.h"


/* USAGE: */
/* COMPILATION: */
/* gcc -Wall -pedantic -o colps hash_collapser.c list_counts.c ht.c -lm */
/* COMPILATION FOR PERFORMANCE ANALYSIS: */
/* gcc -Wall -pedantic -o colps hash_collapser.c list_counts.c ht.c -g -pg -lm */
/* EXECUTION: */
/* ./colps <fastq/fasta_file> <output_file_name> */


KSEQ_INIT(int, read)

#define MAX_SEQ_SIZE 10
  
// -----------------------------------------------------------
// -----------------------------------------------------------
// -----------------------------------------------------------
// -----------------------------------------------------------

  void transpose_hash( hti it, list_count **lc, size_t hash_table_size ) {
  int n = 0;
  while (ht_next(&it)) {
    if(  hash_table_size >= 100 &&
	 n % (int) floor((double) hash_table_size/100*2) == 0 ) {

      /* fprintf( stdout, "\t%s\n", it.key); */
      fprintf( stdout, "=" );
      
      fflush( stdout );
    }
    set_seq( lc, it.key, MAX_SEQ_SIZE, *(int*)it.value );
    free(it.value);
    n++;
  }
}
  
// -----------------------------------------------------------
// -----------------------------------------------------------  

void exit_nomem(void) {
    fprintf(stderr, "out of memory\n");
    exit(1);
}

// ###########################################################################

int main(int argc, char **argv)
{
  FILE *fp;
  FILE *output;
  output = fopen(argv[2], "w+");
  kseq_t *seq;
  int n = 0;
  char temp[MAX_SEQ_SIZE+1];
  memset( temp, 0, MAX_SEQ_SIZE+1 );

  size_t hash_table_size = 0;
  ht* counts = ht_create();
  if (counts == NULL) {
    exit_nomem();
  }

  /* char *p; */
  /* long max_iterations = strtol(argv[2], &p, 10); */
  fprintf( stdout, "\n-----------------------------------------------------\n");

  fprintf( stdout, "\tProgram called with: %s\n", argv[1] );
  
  fprintf( stdout, "\n\tCollapsing sequences" );
  fprintf( stdout, "\n\t[" );
  fflush( stdout );
  sleep( 1 );
	
  fp = fopen(argv[1], "r");
  seq = kseq_init(fileno(fp));
	
  while (kseq_read(seq) >= 0) {
    ++n;
    /* if( n == max_iterations ) { */
    /*   break; */
    /* } */
    /* if(  max_iterations >= 100 && */
    /* 	 n % (int) floor((double) max_iterations/100*2) == 0 ) { */
    /*   fprintf( stdout, "=" ); */
    /*   fflush( stdout ); */
    /* } */
    if( n % 1000000 == 0 ) {
      fprintf( stdout, "=" );
      fflush( stdout );
    }
    memcpy( temp, seq->seq.s, MAX_SEQ_SIZE );

    void* value = ht_get(counts, temp);
    if (value != NULL) {
      // Already exists, increment int that value points to.
      int* pcount = (int*)value;
      (*pcount)++;
      continue;
    }
    // Word not found, allocate space for new int and set to 1.
    int* pcount = malloc(sizeof(int));
    if (pcount == NULL) {
      exit_nomem();
    }
    *pcount = 1;
    if (ht_set(counts, temp, pcount) == NULL) {
      exit_nomem();
    }
  }
  fprintf( stdout, "]\n\t%d sequences read\n", n );
  fflush( stdout );

  /* CREATING AN ITERATOR TO TRAVERSE THROUGH */
  /* THE HASH TABLE. */
  hti it = ht_iterator(counts);
  /* INITIALISING A LINKED LIST ; THE ASSOCIATED */
  /* STRUCTURE IS DEFINED IN LIST_COUNTS.C */
  list_count *lc = init_list_counts();
  list_count **lst_cnt = &lc;

  /* WE ARE ONLY COPYING THE CONTENTS OF THE HASH */
  /* TABLE TO THE LINKED LIST. */
  fprintf( stdout, "\n\tTransposing hash table to linked list\n");
  fprintf( stdout, "\t[" );
  fflush( stdout );
  hash_table_size = ht_length( counts );


  transpose_hash( it, lst_cnt, hash_table_size );
  
  fprintf( stdout, "]\n" );
  fflush( stdout );
  /* /\* printf("%d\n", (int)ht_length(counts)); *\/ */

  /* SORTING THE LINKED LIST CONTAINING THE */
  /* SEQUENCES AND THE ASSOCIATED COUNTS. */
  set_head( lst_cnt, list_sort( lc ) );

  /* /\* PRINTING THE SORTED SEQUENCES TO OUTPUT *\/ */
  /* /\* IN FASTA FORMAT. *\/ */
  /* fprintf( stdout, "\n\tPrinting results to fasta\n"); */
  /* fflush( stdout ); */
  /* print_fasta( lst_cnt, output ); */

  /* PRINTING THE SORTED SEQUENCES TO OUTPUT */
  /* IN PLAIN DAT FORMAT. */
  fprintf( stdout, "\n\tPrinting results to dat file\n");
  fflush( stdout );
  print_dat( lst_cnt, output );

  /* FREEING ALL THE MEMORY */
  /* FOR THE LINKED LIST: */
  free_list_content( get_head( lc ));
  free_list( lc );
  /* FOR THE HASH TABLE: */
  ht_destroy(counts);
  /* FOR THE KSEQ STRUCTURE: */
  kseq_destroy(seq);
  fclose(fp);
  printf("\n-----------------------------------------------------\n");
  return 0;
}
