//  Copyright (C) 2015,2016 George Hobbs
// This file is part of the pfits software package
//

/* pfits is free software: you can redistribute it and/or modify 
 * it under the terms of the GNU General Public License as published by 
 * the Free Software Foundation, either version 3 of the License, or 
 * (at your option) any later version. 
 * pfits is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
 * GNU General Public License for more details. 
 * You should have received a copy of the GNU General Public License 
 * along with pfits.  If not, see <http://www.gnu.org/licenses/>. 
*/

//gcc -lm -o pfits_describe pfits_describe.c pfits_setup.c -I/Users/hob044/hob044/software/cfitsio/include/ -L/Users/hob044/hob044/software/cfitsio/lib -lcfitsio 

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "pfits.h"

void basicDescription(dSetStruct **dSet,int nFiles);

int main(int argc,char *argv[])
{
  dSetStruct **dSet;
  int debug=0;
  int i;
  int nFiles=0,nFiles_t=0;

  // Count number of files
  for (i=0;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	nFiles++;
    }

  // Allocate memory for these files
  dSet = (dSetStruct **)malloc(sizeof(dSetStruct *)*nFiles);
  
  // Initialise everything
  for (i=0;i<nFiles;i++) initialise(&(dSet[i]),debug);
  
  // Read inputs
  for (i=0;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	setFilename(argv[++i],dSet[nFiles_t++],debug);
    }
  // Open the file
  for (i=0;i<nFiles;i++) pfitsOpenFile(dSet[i],debug);

  // Load header information
  for (i=0;i<nFiles;i++) pfitsLoadHeader(dSet[i],debug);

  // Do the description
  basicDescription(dSet,nFiles);

  
  // Close the file
  //  pfitsCloseFile(dSet,debug);
  
  // De-allocate the memory
  for (i=0;i<nFiles;i++) deallocateMemory(&dSet[i],debug);
  free(dSet);
}


void basicDescription(dSetStruct **dSet,int nFiles)
{
  int i;
  
  printf("Filename:                ");
  for (i=0;i<nFiles;i++) printf("%-20.20s",dSet[i]->fileName);
  printf("\n");

  printf("Source:                  ");
  for (i=0;i<nFiles;i++) printf("%-20.20s",dSet[i]->head->source);
  printf("\n");

  printf("Obs. freq (MHz):         ");
  for (i=0;i<nFiles;i++) printf("%-20.3f",dSet[i]->head->freq);
  printf("\n");

  printf("Obs. bandwidth (MHz):    ");
  for (i=0;i<nFiles;i++) printf("%-20.3f",dSet[i]->head->bw);
  printf("\n");

  printf("Start time IMJD:         ");
  for (i=0;i<nFiles;i++) printf("%-20d",dSet[i]->head->stt_imjd);
  printf("\n");

  printf("Start time seconds (s):  ");
  for (i=0;i<nFiles;i++) printf("%-20.3f",dSet[i]->head->stt_smjd);
  printf("\n");

  printf("Number of subints:       ");
  for (i=0;i<nFiles;i++) printf("%-20d",dSet[i]->head->nsub);
  printf("\n");

  printf("Number of channels:      ");
  for (i=0;i<nFiles;i++) printf("%-20d",dSet[i]->head->nchan);
  printf("\n");

  printf("Channels bw (MHz):       ");
  for (i=0;i<nFiles;i++) printf("%-20.3f",dSet[i]->head->chanbw);
  printf("\n");

  printf("Number of bits:          ");
  for (i=0;i<nFiles;i++) printf("%-20d",dSet[i]->head->nbits);
  printf("\n");

  printf("Number of polarisations: ");
  for (i=0;i<nFiles;i++) printf("%-20d",dSet[i]->head->npol);
  printf("\n");

  printf("Number of samples/block: ");
  for (i=0;i<nFiles;i++) printf("%-20d",dSet[i]->head->nsblk);
  printf("\n");

  printf("Sample time (s):         ");
  for (i=0;i<nFiles;i++) printf("%-20.7f",dSet[i]->head->tsamp);
  printf("\n");

  printf("Observation length (s)   ");
  for (i=0;i<nFiles;i++) printf("%-20.3f",dSet[i]->head->tsamp*dSet[i]->head->nsblk*dSet[i]->head->nsub);
  printf("\n");
} 
