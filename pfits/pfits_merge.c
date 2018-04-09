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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fitsio.h"

// In order to merge multiple psrfits files we:
// 1) make a copy of the first file

int main(int argc,char *argv[])
{
  char outfile[1024],oname[1024];
  int i;
  fitsfile *outfptr,*infptr;
  char fname[100][1024];
  int status=0;
  int nFiles=0;
  int ii=1;
  long nSubint_out,nSubint_in;
  int colnum;
  int tableWidth;
  unsigned char *loadBytes;
  double *offsSub;
  double nulldouble = 0.0;
  int initflag = 0;
  double stepSize;
  
  printf("Starting\n");
  // Set up some defaults
  strcpy(outfile,"pfitsMerge.sf");

  for (i=1;i<argc;i++)
    {
      if (strcmp(argv[i],"-o")==0)
	strcpy(outfile,argv[++i]);
      else
	  strcpy(fname[nFiles++],argv[i]);
    }
  
  printf("Starting\n");
  sprintf(oname,"!%s",outfile);
  printf("Creating new PSRFITS file\n");
  if ( !fits_create_file(&outfptr, oname, &status) )
    {
      // Copy the first file
      if ( !fits_open_file(&infptr, fname[0], READONLY, &status) )
	{
	  /* Copy every HDU until we get an error */
	  while( !fits_movabs_hdu(infptr, ii++, NULL, &status) )
	    {
	      fits_copy_hdu(infptr, outfptr, 0, &status);
	    }
	  /* Reset status after normal error */
	  if (status == END_OF_FILE) status = 0;	  
	  fits_close_file(infptr,  &status);	  
	}
      else
	{
	  printf("Unable to open input file %s\n",fname[0]);
	  exit(1);
	}
    }
  else
    printf("Unable to open the output file: %s\n",outfile);

  // Now combine this file with the other files
  printf("combining: nFiles = %d\n",nFiles);
  
  for (i=1;i<nFiles;i++)
    {
      printf("Merging file: %s\n",fname[i]);
      if ( !fits_open_file(&infptr, fname[i], READONLY, &status) )
	{
	  fits_movnam_hdu(outfptr,BINARY_TBL,"SUBINT",1,&status);
	  fits_get_num_rows(outfptr,&nSubint_out,&status);
	  printf("nSubint_out = %ld\n",nSubint_out);

	  fits_movnam_hdu(infptr,BINARY_TBL,"SUBINT",1,&status);
	  fits_get_num_rows(infptr,&nSubint_in,&status);

	  // Make more space in the output file
	  fits_insert_rows(outfptr,nSubint_out,nSubint_in,&status);

	  // Determine the width of the SUBINT table in bytes
	  fits_read_key(infptr,TINT,"NAXIS1",&tableWidth,NULL,&status);
	  
	  printf("nSubint_in = %ld %d\n",nSubint_in,tableWidth);
	  // Allocate enough memory
	  loadBytes = (unsigned char *)malloc(sizeof(unsigned char)*tableWidth*nSubint_in);
	  fits_read_tblbytes(infptr,1,1,tableWidth*nSubint_in,&loadBytes[0],&status);
	  fits_write_tblbytes(outfptr,nSubint_out+1,1,tableWidth*nSubint_in,&loadBytes[0],&status);
	  free(loadBytes);
	  fits_close_file(infptr,&status);
	}
      else
	{
	  printf("Unable to open file %s for merging\n",fname[i]);
	  exit(1);
	}
    }
  // Must now fix the OFFS_SUB column to ensure continuity throughout the PSRFITS file
  printf("Updating OFFS_SUB column\n");
  fits_get_num_rows(outfptr,&nSubint_out,&status);
  fits_get_colnum(outfptr,CASEINSEN,"OFFS_SUB",&colnum,&status);
  offsSub = (double *)malloc(sizeof(double)*nSubint_out);
  fits_read_col(outfptr,TDOUBLE,colnum,1,1,nSubint_out,&nulldouble,&(offsSub[0]),&initflag,&status);
  stepSize = offsSub[0]*2;
  for (i=1;i<nSubint_out;i++)
    offsSub[i] = offsSub[i-1]+stepSize;
  fits_write_col(outfptr,TDOUBLE,colnum,1,1,nSubint_out,&(offsSub[0]),&status);
  
  free(offsSub);
  printf("Closing file\n");
  fits_close_file(outfptr,  &status);	  
  printf("Ending\n");
}
