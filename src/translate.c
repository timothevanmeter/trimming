#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "kseq.h"


/* ---------------------------------------------------------------- */
/* USAGE: */
/* TRANSLATES EACH NUCLEOTIDE SEQUENCE OF A FASTQ/FASTA FILE */
/* (PLAIN OR GZIPPED) INTO AMINO ACIDS, IN FRAME 1 ONLY. */
/* STOP CODONS ARE WRITTEN AS '*', CODONS WITH AN AMBIGUOUS */
/* BASE AS 'X', AND THE 1-2 TRAILING NUCLEOTIDES ARE DROPPED. */
/* ---------------------------------------------------------------- */
/* THE GENETIC CODE USED HERE IS THE STANDARD CODE LISTED IN NCBI: */
/* https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=cgencodes#SG1 */
/* ---------------------------------------------------------------- */
/* COMPILATION: */
/* gcc -O2 -Wall -pedantic -o translate translate.c -lz */
/* ---------------------------------------------------------------- */
/* EXECUTION: */
/* ./translate [-p] <fastq/fasta_file> [output_file] */
/*   -p : PLAIN OUTPUT, ONE PEPTIDE PER LINE (DEFAULT IS FASTA) */
/*   WITHOUT output_file THE RESULT IS WRITTEN TO STDOUT */
/* ---------------------------------------------------------------- */


KSEQ_INIT(gzFile, gzread)

#define INVALID 4

// -----------------------------------------------------------
// -----------------------------------------------------------

/* NUCLEOTIDE -> 2-BIT CODE: A=0, C=1, G=2, T/U=3, ANYTHING ELSE=4 */
static uint8_t nt2code[256];

/* STANDARD GENETIC CODE, INDEXED BY (c0<<4)|(c1<<2)|c2 IN ACGT ORDER */
static const char codon2aa[65] =
  "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF";

void init_nt2code(void) {
  memset( nt2code, INVALID, sizeof(nt2code) );
  nt2code['A'] = nt2code['a'] = 0;
  nt2code['C'] = nt2code['c'] = 1;
  nt2code['G'] = nt2code['g'] = 2;
  nt2code['T'] = nt2code['t'] = 3;
  nt2code['U'] = nt2code['u'] = 3;
}

// -----------------------------------------------------------
// -----------------------------------------------------------

/* TRANSLATES n_codons CODONS OF nt INTO aa, WHICH MUST HOLD */
/* AT LEAST n_codons+1 BYTES. */
void translate_seq( const char *nt, size_t n_codons, char *aa ) {
  size_t k;
  for( k = 0; k < n_codons; k++ ) {
    const unsigned char *p = (const unsigned char*) nt + 3*k;
    uint8_t c0 = nt2code[p[0]];
    uint8_t c1 = nt2code[p[1]];
    uint8_t c2 = nt2code[p[2]];
    if( (c0 | c1 | c2) & INVALID )
      aa[k] = 'X';
    else
      aa[k] = codon2aa[(c0 << 4) | (c1 << 2) | c2];
  }
  aa[n_codons] = '\0';
}

// -----------------------------------------------------------
// -----------------------------------------------------------

void exit_nomem(void) {
  fprintf(stderr, "out of memory\n");
  exit(1);
}

void usage(const char *prog) {
  fprintf(stderr, "usage: %s [-p] <fastq/fasta_file> [output_file]\n", prog);
  exit(1);
}

// ###########################################################################

int main(int argc, char **argv)
{
  int plain = 0;
  int argi = 1;
  if( argi < argc && strcmp(argv[argi], "-p") == 0 ) {
    plain = 1;
    argi++;
  }
  if( argi >= argc || argc - argi > 2 )
    usage(argv[0]);

  gzFile fp = gzopen(argv[argi], "r");
  if( fp == NULL ) {
    fprintf(stderr, "cannot open %s\n", argv[argi]);
    return 1;
  }
  FILE *output = stdout;
  if( argc - argi == 2 ) {
    output = fopen(argv[argi+1], "w");
    if( output == NULL ) {
      fprintf(stderr, "cannot open %s\n", argv[argi+1]);
      return 1;
    }
  }

  init_nt2code();

  /* ONE OUTPUT BUFFER, ONLY GROWN WHEN A LONGER READ APPEARS. */
  char *aa = NULL;
  size_t aa_cap = 0;
  long n = 0;
  int ret;

  kseq_t *seq = kseq_init(fp);
  while ((ret = kseq_read(seq)) >= 0) {
    size_t n_codons = seq->seq.l / 3;
    if( n_codons + 1 > aa_cap ) {
      aa_cap = n_codons + 1;
      aa = realloc(aa, aa_cap);
      if( aa == NULL )
	exit_nomem();
    }
    translate_seq( seq->seq.s, n_codons, aa );

    if( plain )
      fprintf( output, "%s\n", aa );
    else if( seq->comment.l )
      fprintf( output, ">%s %s\n%s\n", seq->name.s, seq->comment.s, aa );
    else
      fprintf( output, ">%s\n%s\n", seq->name.s, aa );
    n++;
  }
  if( ret < -1 )
    fprintf(stderr, "warning: truncated or malformed record after %ld sequences\n", n);
  fprintf( stderr, "%ld sequences translated\n", n );

  /* FREEING ALL THE MEMORY */
  free(aa);
  kseq_destroy(seq);
  gzclose(fp);
  if( output != stdout )
    fclose(output);
  return 0;
}
