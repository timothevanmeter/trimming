// -----------------------------------------------------------
// -----------------------------------------------------------
// -----------------------------------------------------------
// -----------------------------------------------------------
  
typedef struct seq_count seq_count;

seq_count * init_seq_count( char *sequence, int size, int count );

void free_seq_count( seq_count *trash );

void print_seq( seq_count *sc );

// -----------------------------------------------------------
// -----------------------------------------------------------

typedef struct list_count list_count;

seq_count * get_head( list_count *lc );

list_count * init_list_counts(void);

void free_list( list_count *lc );

void free_list_content( seq_count *head );

void set_seq( list_count **lc, char *seq, int size, int count );

void set_head( list_count **lc, seq_count *unit );

void add_seq( list_count **lc, char *seq, int size );

void print_list( list_count **lc );

void print_fasta( list_count **lc, FILE *output );

void print_dat( list_count **lc, FILE *output );

// -----------------------------------------------------------
// -----------------------------------------------------------

seq_count * list_sort( list_count *lc );

/* ########################################################################### */
