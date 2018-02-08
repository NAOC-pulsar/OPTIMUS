// Routine to display the header information from a simulated binary data file
// gcc -O3 -lm -o inspectBinaryFile inspectBinaryFile.c simulate.c T2toolkit.c

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "simulate.h"

int main(int argc,char *argv[])
{
  char fname[1024];
  int i;
  header *head;
  FILE *fin;
   
  head = (header *)malloc(sizeof(header));
  for (i=0;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	strcpy(fname,argv[++i]);
    }
  fin = fopen(fname,"rb");
  simulateReadHeaderParameters(head,fin);
  simulateDisplayHeaderParameters(head);
  fclose(fin);
  
  
  free(head);
  
}
