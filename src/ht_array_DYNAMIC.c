// Based on ht by Ben Hoyt (https://github.com/benhoyt/ht),
// MIT License, Copyright (c) 2021 Ben Hoyt: see src/LICENSE-ht.txt
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <assert.h>
#include "ht_array_DYNAMIC.h"

#define INITIAL_CAPACITY 16

extern int TOTAL_FILE_NUMBER;

typedef struct {
  char* key;  // key is NULL if this slot is empty
  int* counts;
  size_t size; // size of the array containing the counts
} ht_entry;

// Hash table structure: create with ht_create, free with ht_destroy.
struct ht {
  ht_entry* entries;  // hash slots
  size_t capacity;    // size of entries array
  size_t length;      // number of items in hash table
};

ht* ht_create(void) {
  // Allocate space for hash table struct.
  ht* table = malloc(sizeof(ht));
  if (table == NULL) {
    return NULL;
  }
  table->length = 0;
  table->capacity = INITIAL_CAPACITY;

  // Allocate (zero'd) space for entry buckets.
  table->entries = calloc(table->capacity, sizeof(ht_entry));
  if (table->entries == NULL) {
    free(table); // error, free table before we return!
    return NULL;
  }
  return table;
}

void ht_destroy(ht* table) {
  // First free allocated keys.
  for (size_t i = 0; i < table->capacity; i++) {
    free(table->entries[i].key);
    free(table->entries[i].counts);
  }
  // Then free entries array and table itself.
  free(table->entries);
  free(table);
}

#define FNV_OFFSET 14695981039346656037UL
#define FNV_PRIME 1099511628211UL

// Return 64-bit FNV-1a hash for key (NUL-terminated). See description:
// https://en.wikipedia.org/wiki/Fowler–Noll–Vo_hash_function
static uint64_t hash_key(char* key) {
  uint64_t hash = FNV_OFFSET;
  for (char* p = key; *p; p++) {
    hash ^= (uint64_t)(unsigned char)(*p);
    hash *= FNV_PRIME;
  }
  return hash;
}

// ------------------------------------------------- 
void print_counts(int * array, size_t size) {
  int i = 0;
  while(i < size+1) {
    printf("%d,", array[i]);
    i++;
  }
  printf("\n");
}
// ------------------------------------------------- 

void add_count(int** counts, int new_count, int file_number) {
  (*counts)[file_number] = new_count;
}

// ------------------------------------------------- 
// ------------------------------------------------- 
int * create_array(int file_number, int *count) {
  int *array = malloc(TOTAL_FILE_NUMBER * sizeof(int));
  if (array == NULL) {
    fprintf(stderr, "Memory allocation failed\n");
    exit(1);
  }
  for (int i = 0; i < TOTAL_FILE_NUMBER; i++) {
    array[i] = 0;
  }
  array[file_number] = *count;
  return array;
}

// ------------------------------------------------- 
// ------------------------------------------------- 
// Internal function to set an entry (without expanding table). 
char* ht_set_entry(ht_entry* entries, size_t capacity,
		   char* key, int* count, size_t* plength, int file_number, int copy) {
  // ----------------------------------------------
  // AND hash with capacity-1 to ensure it's within entries array.
  uint64_t hash = hash_key(key);
  size_t index = (size_t)(hash & (uint64_t)(capacity - 1));
  // ----------------------------------------------
  if(copy == 1) {
    while(entries[index].key != NULL) {
      index++;
      if(index >= capacity)
	index = 0;
    }
    entries[index].key = (char*)key;
    entries[index].size = file_number;
    entries[index].counts = count;
    return key;
  }
  // ----------------------------------------------
  // Loop till we find an empty entry.
  while (entries[index].key != NULL) {
    if (strcmp(key, entries[index].key) == 0) {
      // Found key (it already exists), update value.
      add_count(&entries[index].counts, *count, file_number-1);
      entries[index].size = TOTAL_FILE_NUMBER;
      return entries[index].key;
    }
    // Key wasn't in this slot, move to next (linear probing).
    index++;
    if (index >= capacity) {
      // At end of entries array, wrap around.
      index = 0;
    }
  }
  // Didn't find key, allocate+copy if needed, then insert it.
  if (plength != NULL) {
    key = strdup(key);
    if (key == NULL) {
      return NULL;
    }
    (*plength)++;
  }
  entries[index].key = (char*)key;
  entries[index].counts = create_array(file_number-1, count);
  entries[index].size = TOTAL_FILE_NUMBER;
  return key;
}

// Expand hash table to twice its current size. Return true on success,
// false if out of memory.
static bool ht_expand(ht* table) {
  // fprintf(stdout, "\tEXPANDING THE HASH TABLE\n");
  // fflush(stdout);
  // Allocate new entries array.
  // fprintf(stdout, "\n\tEXPANDING HASH TABLE\n");
  size_t new_capacity = table->capacity * 2;
  if (new_capacity < table->capacity) {
    return false;  // overflow (capacity would be too big)
  }
  ht_entry* new_entries = calloc(new_capacity, sizeof(ht_entry));
  if (new_entries == NULL) {
    return false;
  }
  // Iterate entries, move all non-empty ones to new table's entries.
  for (size_t i = 0; i < table->capacity; i++) {
    ht_entry entry = table->entries[i];
    if (entry.key != NULL) {
      ht_set_entry(new_entries, new_capacity, entry.key,
		   entry.counts, NULL, entry.size, 1);
    }
  }
  // Free old entries array and update this table's details.
  free(table->entries);
  table->entries = new_entries;
  table->capacity = new_capacity;
  return true;
}

char* ht_set(ht* table, char* key, int* count, int file_number) {
  assert(count != NULL);
  if (count == NULL) {
    return NULL;
  }
  // If length will exceed half of current capacity, expand it.
  if (table->length >= table->capacity / 2) {
    if (!ht_expand(table)) {
      return NULL;
    }
  }
  // Set entry and update length.
  return ht_set_entry(table->entries, table->capacity, key, count,
		      &table->length, file_number, 0);
}

size_t ht_length(ht* table) {
  return table->length;
}

hti ht_iterator(ht* table) {
  hti it;
  it._table = table;
  it._index = 0;
  return it;
}

bool ht_next(hti* it) {
  // Loop till we've hit end of entries array.
  ht* table = it->_table;
  while (it->_index < table->capacity) {
    size_t i = it->_index;
    it->_index++;
    if (table->entries[i].key != NULL) {
      // Found next non-empty item, update iterator key and value.
      ht_entry entry = table->entries[i];
      it->key = entry.key;
      it->counts = entry.counts;
      it->size = entry.size;
      return true;
    }
  }
  return false;
}

// ------------------------------------------------- 
// ------------------------------------------------- 
void add(int ** array, int new) {
  int * temp = (*array);
  while((*temp) != 0) {
    temp++;
  }
  (*temp) = new;
  array = &temp;
}

// ###########################################################################
// int main(int argc, char **argv) 
// { 
//   int * array = calloc(10, sizeof(long int)); 
//   int ** ar = &array; 
//   // ----------------------------------------------- 
//   array[0] = 1; 
//   for(int i=1; i < 10; i++) { 
//     add(ar, array[i-1] * i); 
//   } 
//   printf( "\n--------------------\n" ); 
//   print_counts(array, 10); 
//   // ----------------------------------------------- 
//   // ----------------------------------------------- 
//   ht *table = ht_create(); 
//   char *k1 = "AA"; 
//   char *k2 = "BB"; 
//   int **a = NULL; 
//   for(int i=0; i < 10; i++) { 
//     int *temp = &array[i]; 
//     printf("%d ", i+1); 
//     ht_set(table, k1, temp, i+1); 
//   } 
//   printf("\n"); 
//   // ------------------------ 
//   int temp = 14; 
//   int *t = &temp; 
//   ht_set(table, k2, t, 2); 
//   // ------------------------   
//   int file_number = 4; 
//   int temp2 = 14; 
//   int *t2 = &temp2; 
//   a = create_array(file_number, t2); 
//   int *a_ptr = (*a); 
//   for(int i=0; i < file_number; i++) { 
//     printf("%d,", (*a)[i]); 
//   } 
//   printf("\n"); 
//   // ------------------------ 
//   hti iterator = ht_iterator(table); 
//   hti *it = &iterator; 
//   while(ht_next(it)) { 
//     printf("%s,", it->key); 
//     print_counts(it->counts, it->size); 
//   } 
//   // ------------------------ 
//   free(a_ptr); 
//   a = NULL; 
//   free(array); 
//   ar = NULL; 
//   ht_destroy(table); 
//   // ------------------------ 
//   return 0; 
// } 
