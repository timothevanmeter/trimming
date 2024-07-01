#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "list_counts.h"
  
// -----------------------------------------------------------
// -----------------------------------------------------------
// -----------------------------------------------------------
// -----------------------------------------------------------

typedef struct seq_count {
  char *seq;
  int count;
  seq_count *next;
}seq_count;

seq_count * init_seq_count( char *sequence, int size, int count ) {
  seq_count *sc = malloc( sizeof( seq_count ));
  sc->seq = (char*) malloc( size * sizeof( char ));
  memset(sc->seq, 0, size);
  strcpy( sc->seq, sequence );
  sc->count = count;
  sc->next = NULL;
  return sc;
}

void free_seq_count( seq_count *trash ) {
  if( trash != NULL ) {
    if( trash->seq != NULL ) {
      free( trash->seq );
    }
    free( trash );
    trash = NULL;
  }
}

void print_seq( seq_count *sc ) {
  if( sc != NULL ) {
    printf( " %s %d\n", sc->seq, sc->count );
  }
}

// -----------------------------------------------------------
// -----------------------------------------------------------

typedef struct list_count {
  seq_count *head;
  seq_count *tail;
}list_count;

seq_count * get_head( list_count *lc ) {
  return lc->head;
}

list_count * init_list_counts() {
  list_count *lc;
  lc = malloc( sizeof( list_count ));
  lc->head = NULL;
  lc->tail = NULL;
  return lc;
}

void free_list( list_count *lc ) {
  free( lc );
  lc = NULL;
}

void free_list_content( seq_count *head ) {
  while( head != NULL ){
    seq_count *temp = head;
    head = head->next;    
    free_seq_count( temp );
  }
}


/* USED IN THE TRANSPOSE_HASH FUNCTION */
/* No need to compare the sequences here, which are already */
/*   collapsed in the hash table. */
/*   To be more efficient the head of the list is not past, but */
/*   the last added unit instead.  */
void set_seq( list_count **lc, char *seq, int size, int count ) {
  if( !(*lc)->head ) {
    (*lc)->head= init_seq_count( seq, size, count );
    (*lc)->tail = (*lc)->head;
  } else {
    seq_count *temp = (*lc)->tail;
    seq_count *new = init_seq_count( seq, size, count );
    /* temp = (*lc)->head; */
    while( temp->next != NULL ) {
      temp = temp->next;
    }
    temp->next = new;
    (*lc)->tail = temp->next;
  }
}


void set_head( list_count **lc, seq_count *new_head ) {
  (*lc)->head = new_head;
}

void add_seq( list_count **lc, char *seq, int size ) {
  if( !(*lc)->head ) {
    (*lc)->head= init_seq_count( seq, size, 1 );    
  } else {
    seq_count *temp = (*lc)->head;
    seq_count *new = init_seq_count( seq, size, 1 );
    int add = 0;
    temp = (*lc)->head;
    while( temp->next != NULL ) {
      // THE SEQUENCE IS ALREADY STORED IN THE LIST
      // SO WE JUST INCREMENT THE COUNT.
      if( strcmp( seq, temp->next->seq ) == 0 ) {
	temp->next->count++;
	// NO NEED FOR MORE SPACE IN MEMORY!
	free_seq_count( new );
	add = 1;
	break;
      }
      temp = temp->next;
    }
    if( add == 0 ){
      temp->next = new;
    }
  }
}

void print_list( list_count **lc ) {
  seq_count *temp = (*lc)->head;
  while( temp != NULL ) {
    if( temp->count != 1 ) {
      print_seq( temp );
    }
    temp=temp->next;
  }
}

void print_fasta( list_count **lc, FILE *output ) {
  seq_count *temp = (*lc)->head;
  while( temp != NULL ) {
    fprintf( output, "\n> %d", temp->count );
    fprintf( output, "\n%s", temp->seq );
    temp=temp->next;
  }
  fclose( output );
}

void print_dat( list_count **lc, FILE *output ) {
  seq_count *temp = (*lc)->head;
  while( temp != NULL ) {
    fprintf( output, "\n%s %d", temp->seq, temp->count );
    temp=temp->next;
  }
  fclose( output );
}

// -----------------------------------------------------------
// -----------------------------------------------------------

int cmp( seq_count *a, seq_count *b) {
  if (a->count < b->count)
    return -1;
  if (a->count > b->count)
    return +1;
  return 0;
}

seq_count * list_sort( list_count *lc ) {
  /* list_count *list = *lc; */
  seq_count *list_head = lc->head;
  seq_count *p, *q, *e, *tail; //, *oldhead;
  int insize, nmerges, psize, qsize, i;

  /*
   * Silly special case: if `list' was passed in as NULL, return
   * NULL immediately.
   */
  if( !list_head )
    return NULL;

  insize = 1;

  while( 1 ) {
    p = list_head;
    list_head = NULL;
    tail = NULL;

    nmerges = 0;  /* count number of merges we do in this pass */

    while( p ) {
      nmerges++;  /* there exists a merge to be done */
      /* step `insize' places along from p */
      q = p;
      psize = 0;
      for( i = 0; i < insize; i++ ) {
	psize++;
	q = q->next;
	if( !q ) break;
      }

      /* if q hasn't fallen off end, we have two list_heads to merge */
      qsize = insize;

      /* now we have two list_heads; merge them */
      while( psize > 0 || (qsize > 0 && q) ) {

	/* decide whether next element of merge comes from p or q */
	if( psize == 0 ) {
	  /* p is empty; e must come from q. */
	  e = q; q = q->next; qsize--;
	} else if (qsize == 0 || !q) {
	  /* q is empty; e must come from p. */
	  e = p; p = p->next; psize--;
	} else if (cmp(p,q) <= 0) {
	  /* First element of p is lower (or same);
	   * e must come from p. */
	  e = q; q = q->next; qsize--;
	} else {
	  /* First element of q is lower; e must come from q. */
	  e = p; p = p->next; psize--;
	}

	/* add the next element to the merged list_head */
	if (tail) {
	  tail->next = e;
	} else {
	  list_head = e;
	}
	tail = e;
      }

      /* now p has stepped `insize' places along, and q has too */
      p = q;
    }
    tail->next = NULL;

    /* If we have done only one merge, we're finished. */
    if (nmerges <= 1)   /* allow for nmerges==0, the empty list_head case */
      return list_head;

    /* Otherwise repeat, merging list_heads twice the size */
    insize *= 2;
  }
}



/* ########################################################################### */
/* ########################################################################### */
/* ########################################################################### */
/* ########################################################################### */



/* int main(int argc, char **argv) */
/* { */
/*   FILE *fp;	 */
/*   kseq_t *seq; */
/*   int n = 0; //, slen = 0; //, qlen = 0, init = 0; */
/*   list_count *list_count = init_list_counts(); */
/*   list_count **lst_cnt = &list_count; */
/*   char temp[31]; */
/*   memset( temp, 0, 31 ); */

/*   char *p; */
/*   long max_iterations = strtol(argv[2], &p, 10); */
/*   /\* long total_seqs = strtol(argv[3], &p, 10); *\/ */
/*   /\* printf( "\nTOTAL SEQS = %ld\n", total_seqs ); *\/ */
/*   fprintf( stdout, "\n-----------------------------------------------------\n"); */

/*   fprintf( stdout, "\tProgram called with: %s\n", argv[1] ); */
  
/*   fprintf( stdout, "\n\tCollapsing sequences" ); */
/*   fprintf( stdout, "\n\t[" ); */
/*   fflush( stdout ); */
/*   sleep( 1 ); */
	
/*   fp = fopen(argv[1], "r"); */
/*   seq = kseq_init(fileno(fp)); */
	
/*   while (kseq_read(seq) >= 0) { */
/*     ++n; */
/*     /\* if(  n % (int) floor((((double) total_seqs/100))*2) == 0 ) { *\/ */
/*     /\*   fprintf( stdout, "=" ); *\/ */
/*     /\*   fflush( stdout ); *\/ */
/*     /\* } *\/ */
/*     /\* if( n % 1000 == 0 ) { *\/ */
/*     /\*   list_count->head = list_sort( lst_cnt ); *\/ */
/*     /\* } *\/ */
/*     // slen = seq->seq.l; //, qlen = seq->qual.l; */
/*     /\* seqstr = seq->seq.s; *\/ */
/*     //    REMINDER!!!!! */
/*     // THE FOLLOWING TRIMS THE SEQUENCE TO THE FIRST TEN */
/*     // NUCLOETIDES JUST FOR CLARITY DURING TESTING. */
/*     memcpy( temp, seq->seq.s, 30 ); */
/*     add_seq( lst_cnt, temp, 31 ); */
/*     if( n == max_iterations ) { */
/*       break; */
/*     } */
/*   } */
/*   fprintf( stdout, "]\n" ); */
/*   fprintf(stdout, "\n-----------------------------------------------------\n"); */
/*   fflush( stdout ); */

/*   list_count->head = list_sort( lst_cnt ); */
/*   print_list( lst_cnt ); */
/*   /\* print_fasta( lst_cnt ); *\/ */

/*   free_list_content( list_count->head ); */
/*   free_list( list_count ); */

/*   kseq_destroy(seq); */
/*   fclose(fp); */
/*   printf("\n-----------------------------------------------------\n"); */
/*   return 0; */
/* } */
