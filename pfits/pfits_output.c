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

//gcc -lm -o pfits_statistcs pfits_statistics.c pfits_statistics.c -I/Users/hob044/hob044/software/cfitsio/include/ -L/Users/hob044/hob044/software/cfitsio/lib -lcfitsio 

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "pfits.h"
#include "fitsio.h"
#include "header.h"


void writeFilterbankHeader(FILE *fout,dSetStruct *dSet);
void swapThebytes(unsigned char *data,int n);
void send_string(char *string,FILE *output); 
void send_int(char *name, int integer,FILE *output); 
void send_double (char *name, double double_precision,FILE *output); 

int main(int argc,char *argv[])
{
  dSetStruct *dSet;
  int debug=0;
  int status=0;
  int i,j;
  unsigned char *data;
  char outName[128];
  long isub0,isub1;
  FILE *fout;
  int colnum;
  int samplesperbyte;
  unsigned char nval = '0';
  int initflag=0;
  int swapbytes=0;
  
  // Initialise everything
  initialise(&dSet,debug);
  
  // Read inputs
  for (i=0;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	setFilename(argv[++i],dSet,debug);
      else if (strcmp(argv[i],"-o")==0)
	strcpy(outName,argv[++i]);
      else if (strcmp(argv[i],"-swapbytes")==0)
	swapbytes=1;
    }

  pfitsOpenFile(dSet,debug);
  pfitsLoadHeader(dSet,debug);
  printf("RA = %s\n",dSet->head->ra);
  isub0 = 0;
  isub1 = dSet->head->nsub;
  samplesperbyte = 8/dSet->head->nbits;
  
  data = (unsigned char *)malloc(sizeof(unsigned char)*dSet->head->nsblk*dSet->head->nchan*dSet->head->npol/samplesperbyte);

  // Open output file
  if (!(fout = fopen(outName,"wb")))
    {
      printf("Unable to open output file: %s\n",outName);
      exit(1);
    }

  // Write header for filterbank file
  writeFilterbankHeader(fout,dSet);


  fits_movnam_hdu(dSet->fp,BINARY_TBL,"SUBINT",1,&status);
  if (status){printf("Unable to find the SUBINT table\n"); exit(1);}
  fits_get_colnum(dSet->fp,CASEINSEN,"DATA",&colnum,&status);

  for (i=isub0;i<isub1;i++)
    {
      fits_read_col_byt(dSet->fp,colnum,i+1,1,dSet->head->npol*dSet->head->nsblk*dSet->head->nchan/samplesperbyte,nval,data,&initflag,&status);
      if (swapbytes==1)
	swapThebytes(data,dSet->head->npol*dSet->head->nsblk*dSet->head->nchan/samplesperbyte);
      fwrite(data, 1, dSet->head->npol*dSet->head->nsblk*dSet->head->nchan/samplesperbyte, fout);
    }
  fclose(fout);
  


  free(data);

  // De-allocate the memory
  deallocateMemory(&dSet,debug);
}

void swapThebytes(unsigned char *data,int n)
{
  int i;

  for (i=0;i<n;i++)
    {
      data[i] = (data[i] & 0xF0) >> 4 | (data[i] & 0x0F) << 4;
      data[i] = (data[i] & 0xCC) >> 2 | (data[i] & 0x33) << 2;
      data[i] = (data[i] & 0xAA) >> 1 | (data[i] & 0x55) << 1;
    }
}

void writeFilterbankHeader(FILE *fout,dSetStruct *dSet)
{
  strcpy(source_name,"simulated");
  machine_id = 0;
  telescope_id = 7;
  fch1 = dSet->head->chanFreq[0];
  foff = dSet->head->chanbw;
  nbits = dSet->head->nbits;
  nchans = dSet->head->nchan;
  tstart = 53000.0;
  tsamp = dSet->head->tsamp;
  nifs = dSet->head->npol;

    // header
  send_string("HEADER_START",fout);
  send_string("source_name",fout);
  send_string(source_name,fout);
  send_int("machine_id",machine_id,fout);
  send_int("telescope_id",telescope_id,fout);
  if (dSet->head->nchan==1) // If nchan = 1 (i.e., dedispersed)
    {
      send_int("data_type",2,fout);
      send_double("refdm",0,fout); // SHOULD SET THIS PROPERLY
    }
  else
    send_int("data_type",1,fout);


  {
    int deg,hr,min;
    float sec;
    char temp[1024];

    sscanf(dSet->head->ra,"%d:%d:%f",&hr,&min,&sec);
    sprintf(temp,"%02d%02d%.1f",hr,min,sec);
    sscanf(temp,"%lf",&src_raj);

    sscanf(dSet->head->dec,"%d:%d:%f",&deg,&min,&sec);
    sprintf(temp,"%02d%02d%.1f",deg,min,sec);
    sscanf(temp,"%lf",&src_dej);
  }

  send_double("src_raj",src_raj,fout);
  send_double("src_dej",src_dej,fout);

  send_double("fch1",fch1,fout);
  send_double("foff",foff,fout);
  send_int("nchans",nchans,fout);
  send_int("nbits",nbits,fout);
  send_double("tstart",tstart,fout);
  send_double("tsamp",tsamp,fout);
  send_int("nifs",nifs,fout);
  send_string("HEADER_END",fout);
}

void send_string(char *string,FILE *output) 
{
  int len;
  len=strlen(string);
  fwrite(&len, sizeof(int), 1, output);
  fwrite(string, sizeof(char), len, output);
  /*fprintf(stderr,"%s\n",string);*/
}

void send_int(char *name, int integer,FILE *output) 
{
  send_string(name,output);
  fwrite(&integer,sizeof(int),1,output);
  /*fprintf(stderr,"%d\n",integer);*/
}

void send_double (char *name, double double_precision,FILE *output) 
{
  send_string(name,output);
  fwrite(&double_precision,sizeof(double),1,output);
  /*fprintf(stderr,"%f\n",double_precision);*/
}

