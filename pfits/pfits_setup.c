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

// Standard routines for pfits

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "pfits.h"
#include "fitsio.h"

int initialise(dSetStruct **dSet,int debugFlag)
{
  if (debugFlag==1) printf("Entering initialise\n");
  
  // Allocate memory for the main structure
  if (!(*dSet = (dSetStruct *)malloc(sizeof(dSetStruct))))
    errorStop("Unable to allocate enough memory to allocate the data structure\n",*dSet,debugFlag);

  // Allocate memory for the header information
  //
  (*dSet)->headerMemorySet=0;
  if (!((*dSet)->head = (headerStruct *)malloc(sizeof(headerStruct))))
    errorStop("Unable to allocate enough memory to allocate the header structure\n",(*dSet),debugFlag);
  (*dSet)->headerMemorySet=1;
  (*dSet)->head->headerSet=0;

  // Set some defaults
  strcpy((*dSet)->fileName,"UNSET");
  (*dSet)->fileOpen = 0;
  (*dSet)->fileSet = 0;
  if (debugFlag==1) printf("Complete initialise\n");
  return 0;
}

int deallocateMemory(dSetStruct **dSet,int debugFlag)
{
  if (debugFlag==1) printf("Entering deallocateMemory\n");
  if ((*dSet)->headerMemorySet==1)
    free((*dSet)->head);
  free(*dSet);
  if (debugFlag==1) printf("Complete deallocateMemory\n");

  return 0;
}

void errorStop(char *err,dSetStruct *dSet,int debugFlag)
{
  printf("%s",err);
  deallocateMemory(&dSet,debugFlag);
  exit(1);
}

void pfitsLoadHeader(dSetStruct *dSet,int debug)
{
  int status=0;
  int colnum=0;
  float nval=0;
  int initFlag=0;
  
  if (debug==1)
    printf("Entering loadHeader\n");

  if (dSet->headerMemorySet==0)
    errorStop("Trying to load the header information before allocating memory\n",dSet,debug);

  // Read primary header information
    if (status) {printf("pos 0\n"); fits_report_error(stderr,status); exit(1);}
  fits_read_key(dSet->fp,TFLOAT,"OBSFREQ",&(dSet->head->freq),NULL,&status);
  if (status) {printf("pos 1 STT_OBSFREQ\n"); fits_report_error(stderr,status); exit(1);}
  fits_read_key(dSet->fp,TFLOAT,"OBSBW",&(dSet->head->bw),NULL,&status);
  if (status) {printf("pos 1 STT_OBSBW\n"); fits_report_error(stderr,status); exit(1);}
  fits_read_key(dSet->fp,TINT,"STT_IMJD",&(dSet->head->stt_imjd),NULL,&status);
  if (status) {printf("pos 1 STT_IMJD\n"); fits_report_error(stderr,status); exit(1);}
  fits_read_key(dSet->fp,TFLOAT,"STT_SMJD",&(dSet->head->stt_smjd),NULL,&status);
  if (status) {printf("pos 1 STT_SMJD\n"); fits_report_error(stderr,status); exit(1);}
  fits_read_key(dSet->fp,TFLOAT,"STT_OFFS",&(dSet->head->stt_offs),NULL,&status);
  if (status) {printf("pos 1 STT_OFFS\n"); fits_report_error(stderr,status); exit(1);}
  fits_read_key(dSet->fp,TSTRING,"SRC_NAME",(dSet->head->source),NULL,&status);      
  fits_read_key(dSet->fp,TSTRING,"RA",(dSet->head->ra),NULL,&status);      
  if (status) {printf("pos 1 RA\n"); fits_report_error(stderr,status); exit(1);}
  fits_read_key(dSet->fp,TSTRING,"DEC",(dSet->head->dec),NULL,&status);      
  if (status) {printf("pos 1 DEC\n"); fits_report_error(stderr,status); exit(1);}
  
  // Read SUBINT header information
  fits_movnam_hdu(dSet->fp,BINARY_TBL,"SUBINT",1,&status);
  if (status) {printf("pos 1.1\n");  fits_report_error(stderr,status); exit(1);}
  fits_read_key(dSet->fp,TINT,"NAXIS2",&(dSet->head->nsub),NULL,&status);
  if (status) {printf("pos 1.2\n");  fits_report_error(stderr,status); exit(1);}
  fits_read_key(dSet->fp,TINT,"NCHAN",&(dSet->head->nchan),NULL,&status);
  if (status) {printf("pos 1.3\n");  fits_report_error(stderr,status); exit(1);}
  fits_read_key(dSet->fp,TINT,"NBITS",&(dSet->head->nbits),NULL,&status);
  if (status) {printf("pos 1.4\n");  fits_report_error(stderr,status); exit(1);}
  fits_read_key(dSet->fp,TINT,"NPOL",&(dSet->head->npol),NULL,&status);
  if (status) {printf("pos 1.5\n");  fits_report_error(stderr,status); exit(1);}
  fits_read_key(dSet->fp,TINT,"NSBLK",&(dSet->head->nsblk),NULL,&status);
  if (status) {printf("pos 1.6\n");  fits_report_error(stderr,status); exit(1);}
  fits_read_key(dSet->fp,TFLOAT,"TBIN",&(dSet->head->tsamp),NULL,&status);
  if (status) {printf("pos 1.7\n");  fits_report_error(stderr,status); exit(1);}
  fits_read_key(dSet->fp,TFLOAT,"CHAN_BW",&(dSet->head->chanbw),NULL,&status);
  if (status) {printf("pos 2\n");  fits_report_error(stderr,status); exit(1);}
      
  // Read the frequency channels
  fits_get_colnum(dSet->fp,CASEINSEN,"DAT_FREQ",&colnum,&status);  
  printf("Reading\n");
  fits_read_col_flt(dSet->fp,colnum,1,1,dSet->head->nchan,nval,dSet->head->chanFreq,&initFlag,&status);
  if (status) {printf("pos 3\n");  fits_report_error(stderr,status); exit(1);}
	
  if (debug==1)
    printf("Leaving loadHeader\n");
}

void pfitsOpenFile(dSetStruct *dSet,int debug)
{
  int status=0;
  
  if (debug==1)
    printf("Entering pfitsOpenFile\n");

  if (dSet->fileSet==0)
    errorStop("Trying to open a file before setting the filename\n",dSet,debug);

  fits_open_file(&(dSet->fp),dSet->fileName,READONLY,&status);
  fits_report_error(stderr,status);

  dSet->fileOpen = 1;
  
  if (debug==1)
    printf("Leaving pfitsOpenFile\n");
}

void setFilename(char *str,dSetStruct *dSet,int debug)
{
  if (debug==1) printf("Entering setFilename\n");
  strcpy(dSet->fileName,str);
  dSet->fileSet = 1;
  if (debug==1) printf("Leaving setFilename\n");
}
