
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
#include <math.h>
#include "fitsio.h"
#include <cpgplot.h>

void processHDU(fitsfile *fp,int hduNum);
void viewFitsTable(fitsfile *fp,int hduNum,int viewTable);

int main(int argc,char *argv[])
{
  int i;
  char fname[128];
  fitsfile *fp;
  int status=0;
  int nkey=-1;
  int morekeys=-1;
  int keysexist=-1;
  char keyname[128],val[128],comment[128];
  int nhdus;
  int hdutype;
  int hduNum;
  int quick=0;
  int nFile=0;
  int fileID[1024];

  for (i=1;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	{
	  strcpy(fname,argv[++i]);
	  fileID[nFile] = i;
	  nFile++;
	}
      else if (strcmp(argv[i],"-q")==0)
	quick=1;
      else
	{
	  fileID[nFile] = i;
	  nFile++;
	}
    }

  if (quick==1)
    {
      char obsMode[1024];
      char projID[1024];
      char obsFreq[1024];
      char obsBW[1024];
      char obsNCHAN[1024];
      char srcName[1024];
      int i;

      for (i=0;i<nFile;i++)
	{
	  strcpy(fname,argv[fileID[i]]);
	  fits_open_file(&fp,fname,READONLY,&status);
	  fits_get_num_hdus(fp,&nhdus,&status);
	  	  
	  fits_read_key(fp,TSTRING,"OBS_MODE",obsMode,NULL,&status);
	  fits_read_key(fp,TSTRING,"SRC_NAME",srcName,NULL,&status);
	  fits_read_key(fp,TSTRING,"PROJID",projID,NULL,&status);
	  fits_read_key(fp,TSTRING,"OBSFREQ",obsFreq,NULL,&status);
	  fits_read_key(fp,TSTRING,"OBSBW",obsBW,NULL,&status);
	  fits_read_key(fp,TSTRING,"OBSNCHAN",obsNCHAN,NULL,&status);
	  printf("%-20.20s %-15.15s %-8.8s %-5.5s %-5.5s %-5.5s %-5.5s ",fname,srcName,obsMode,projID,obsFreq,obsBW,obsNCHAN);
	  if (strcmp(obsMode,"SEARCH")==0)
	    {
	      char nbits[1024];
	      char npol[1024];
	      char tsamp[1024];

	      fits_movnam_hdu(fp,BINARY_TBL,"SUBINT",1,&status);
	      fits_read_key(fp,TSTRING,"NBITS",nbits,NULL,&status);
	      fits_read_key(fp,TSTRING,"NPOL",npol,NULL,&status);
	      fits_read_key(fp,TSTRING,"TBIN",tsamp,NULL,&status);
	      printf("nbits = %s npol = %s tbin = %s",nbits,npol,tsamp);
	    }

	  printf("\n");
	  fits_close_file(fp,&status);
	}
    }
  else
    {
      fits_open_file(&fp,fname,READONLY,&status);
      fits_get_num_hdus(fp,&nhdus,&status);

      do 
	{
	  printf("------------------------------------------\n");
	  printf("File name: %s\n",fname);
	  printf("Number of header data units = %d\n",nhdus);
	  printf("------------------------------------------\n");
	  for (i=0;i<nhdus;i++)
	    {
	      printf("[%2.2d] ",i+1);
	      fits_movabs_hdu(fp,i+1,&hdutype,&status);
	      //      fits_get_hdu_type(fp,&hdutype,&status);
	      if (i==0)
		printf("PRIMARY ");
	      else
		{
		  fits_read_key(fp,TSTRING,"EXTNAME",keyname,comment,&status);
		  printf("%s ",keyname);
		}
	      if (hdutype == IMAGE_HDU)
		printf("(IMAGE_HDU)");
	      else if (hdutype == ASCII_TBL)
		printf("(ASCII_TBL)");
	      else if (hdutype == BINARY_TBL)
		printf("(BINARY_TBL)");
	      printf("\n");
	    }
	  printf("\n\n");
	  printf("Enter value of HDU to view: ");
	  scanf("%d",&hduNum);
	  if (hduNum > 0) processHDU(fp,hduNum);
	} while (hduNum > 0);
      printf("\n\n Goodbye\n");
      fits_close_file(fp,&status);
    }

}

void processHDU(fitsfile *fp,int hduNum)
{
  int keysexist=-1;
  int morekeys=-1;
  int status=0;
  int i;
  int hdutype;
  char keyname[128],val[128],comment[128];
  int nrows,ncol;
  int typecode;
  long repeat,width;
  char card[1024];
  int colnum;
  int sstat=0;
  int viewTable;

  printf("\n\n");
  fits_movabs_hdu(fp,hduNum,&hdutype,&status);
  fits_get_hdrspace(fp,&keysexist,&morekeys,&status);
  printf("Header information\n\n");
  for (i=0;i<keysexist;i++)
    {
      fits_read_keyn(fp,i+1,keyname,val,comment,&status);
      printf("%-15.15s %-24.24s %s\n",keyname,val,comment);
      if (strcmp(keyname,"NAXIS2")==0)
	sscanf(val,"%d",&nrows);
      else if (strcmp(keyname,"TFIELDS")==0)
	sscanf(val,"%d",&ncol);
    }
  printf("End header information\n");
  if (hduNum > 1)
    {
      do
	{
	  sstat=0;
	  fits_movabs_hdu(fp,hduNum,&hdutype,&status);
	  printf("\n\n");
	  printf("Number of rows: %d\n",nrows);
	  printf("Number of columns: %d\n",ncol);
	  printf("\n\n");
	  for (i=0;i<ncol;i++)
	    {
	      printf("[%2.2d] ",i+1);
	      //      printf(" %d %d %d",(int)typecode,(int)repeat,(int)width);
	      status = sstat;
	      fits_get_colname(fp,0,"*",card,&colnum,&status);
	      sstat = status;
	      status = 0; // As there are more than one matching colnames
	      printf("%s",card);
	      
	      /*	      fits_get_coltype(fp,i+1,&typecode,&repeat,&width,&status);
	      if (typecode == TSTRING)
		printf(" (TSTRING %d)",(int)width);
	      else if (typecode == TSHORT)
		printf(" (TSHORT %d)",(int)width);
	      else if (typecode == TLONG)
		printf(" (TLONG %d)",(int)width);
	      else if (typecode == TFLOAT)
		printf(" (TFLOAT %d)",(int)width);
	      else if (typecode == TDOUBLE)
		printf(" (TDOUBLE %d)",(int)width);
	      else
	      printf(" (UNKNOWN %d)",(int)width); */
	      
	      printf("\n");
	      
	    }
	  
	  printf("\n\n");
	  printf("Select data to view (-1 to return) ");
	  scanf("%d",&viewTable);
	  if (viewTable > 0)
	    viewFitsTable(fp,hduNum,viewTable);
	}
      while (viewTable > 0);
    }
  printf("\n\n");
}

void viewFitsTable(fitsfile *fp,int hduNum,int viewTable)
{
  int i;
  char **line;
  char nval[128]="UNKNOWN";
  int anynul=0;
  int status=0;
  int r1,r2,c1,c2;
  int row,col;
  long width,repeat;
  int typecode;
  int maxdim = 5;
  int naxis;
  long naxes[maxdim];
  int type;
  long nrows;
  long maxsize;

  fits_get_coltype(fp,viewTable,&typecode,&repeat,&width,&status);

  fits_get_num_rows(fp,&nrows,&status);
  fits_read_tdim(fp,viewTable,maxdim,&naxis,naxes,&status);
  maxsize=1;
  for (i=0;i<naxis;i++)
    {
      printf("Dimension %d = %ld\n",i+1,naxes[i]);
      maxsize *= naxes[i];
    }

  //
  //  Should check whether the array has strings
  //  then should change the possible element numbers
  //



  printf("Enter start row number (min: 1, max: %ld) ",nrows); scanf("%d",&r1);
  printf("Enter end row number (min: %d, max: %ld) ",r1,nrows); scanf("%d",&r2);
  if (typecode == TSTRING)
    {
      printf("This is a columns of 'strings'\n");
      c1=c2=1;
    }
  else
    {
      printf("Enter start element number (min: 1, max: %ld) ",maxsize); scanf("%d",&c1);
      printf("Enter end element number  (min: %d, max: %ld) ",c1,maxsize); scanf("%d",&c2);
    }
  if (r1 < 1)    {      printf("Setting start row to 1\n"); r1=1;    }
  if (r2 < 1)    {      printf("Setting end row to 1\n"); r2=1;    }
  if (c1 < 1)    {      printf("Setting start element to 1\n"); c1=1;    }
  if (c2 < 1)    {      printf("Setting end element to 1\n"); c2=1;    }

  printf("\n\n");
  if (typecode == TSTRING)
    printf("(1) List on screen, (2) Print to file: ");
  else
    printf("(1) List on screen, (2) Print to file, (3) Plot: ");
  scanf("%d",&type);


  if (type==1 || type==2)
    {
      FILE *fout;

      if (type==2)
	fout = fopen("output.dat","w");
      else
	fout = stdout;

      line = (char **)malloc(sizeof(char *));
      line[0] = (char *)malloc(sizeof(char)*1024);     

      for (row = r1;row <= r2; row++)
	{
	  for (col=c1;col<=c2;col++)
	    {
	      fits_read_col_str(fp,viewTable,row,col,1,nval,line,&anynul,&status);
	      fprintf(fout,"%-6d %-6d %s\n",row,col,line[0]);
	      if (status)
		{
		  fits_report_error(stderr,status);
		  status=0;
		}
	    }
	}
      if (type==2)
	fclose(fout);
      free(line[0]);
      free(line);
    }
  else if (type==3)
    {
      float fy[(c2-c1+1)*(r2-r1+1)];
      float fx[(c2-c1+1)*(r2-r1+1)];
      //      float nval = 0;
      //      int anynul=0;
      int n=0;
      float minx,maxx,miny,maxy;
      float ominx,omaxx,ominy,omaxy;
      float mx,my;
      char key;
      int drawpt=-1;
      int drawline=1;

      line = (char **)malloc(sizeof(char *));
      line[0] = (char *)malloc(sizeof(char)*1024);     


      cpgbeg(0,"/xs",1,1);
      cpgask(0);
      for (row = r1;row <= r2; row++)
	{
	  for (col=c1;col<=c2;col++)
	    {
	      //	      printf("With %d %d\n",row,col);
	      fits_read_col_str(fp,viewTable,row,col,1,nval,line,&anynul,&status);
	      //	      printf("Have %s\n",line[0]);
	      sscanf(line[0],"%f",&fy[n]);
	      //	      printf("Have %s %g\n",line[0],fy[n]);
	      fx[n] = n;
	      n++;
	    }
	  //	  fits_read_col_flt(fp,viewTable,row,c1,c2-c1,nval,fy,&anynul,&status);
	  //	  for (i=0;i<c2-c1;i++)
	  //	    fx[i] = i;
	  //	  n+=c2-c1;
	}
      free(line[0]);
      free(line);

      minx = maxx = fx[0];
      miny = maxy = fy[0];
      for (i=0;i<n;i++)
	{
	  if (minx > fx[i]) minx = fx[i];
	  if (maxx < fx[i]) maxx = fx[i];
	  if (miny > fy[i]) miny = fy[i];
	  if (maxy < fy[i]) maxy = fy[i];
	}
      ominx = minx;
      omaxx = maxx;
      ominy = miny;
      omaxy = maxy;

      printf("Have %g %g %g %g\n",minx,maxx,miny,maxy);
      printf("Press 'q' to quit plot\n");
      printf("'p' toggle plotting points\n");
      printf("'l' toggle plotting lines\n");
      printf("'z' zoom\n");
      printf("'u' unzoom\n");
      printf("left click: display mouse position\n");
      do {
	cpgenv(minx,maxx,miny,maxy,0,1);
	if (drawpt==1) cpgpt(n,fx,fy,9);
	if (drawline==1) cpgline(n,fx,fy);
	cpgcurs(&mx,&my,&key);
	if (key=='l') drawline*=-1;
	if (key=='p') drawpt*=-1;
	if (key=='z')
	  {
	    float mx2,my2;
	    cpgband(2,0,mx,my,&mx2,&my2,&key);
	    if (mx > mx2) {maxx = mx; minx = mx2;}
	    else {maxx = mx2; minx = mx;}
	    
	    if (my > my2) {maxy = my; miny = my2;}
	    else {maxy = my2; miny = my;}	   
	  }
	else if (key=='u')
	  {
	    minx = ominx;
	    maxx = omaxx;
	    miny = ominy;
	    maxy = omaxy;
	  }
	else if (key=='A')
	  printf("mouseX = %g, mouseY = %g\n",mx,my);
      } while (key != 'q');

      cpgend();
    }
}
