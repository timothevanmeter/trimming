#include <stdio.h>

int getlines(char * filepath) {
  FILE *file;
  char ch;
  int lines = 0;
  file = fopen(filepath, "r");
  while(!feof(file))
    {
      ch = fgetc(file);
      if(ch == '\n')
	{
	  lines++;
	}
    }
  fclose(file);
  return lines;
}
