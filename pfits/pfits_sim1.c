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
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "fitsio.h"
#include "T2toolkit.h"

#define MAX_BURSTS 500
#define MAX_NB_RFI 10

typedef struct simDataStruct {
  char fname[128];
  int npol;
  int nbits;
  int nsblk;
  int nsub;
  int nchan;
  float tsamp;
  float fchan0;
  float chanbw;
} simDataStruct;

typedef struct burstStruct {
  float time;
  float dm;
  float flux;
  float width;
} burstStruct;

typedef struct nb_rfi_struct {
  float t1;
  float t2;
  float f1;
  float f2;
  float flux;
} nb_rfi_struct;
  

void createSimulation(simDataStruct *simData,burstStruct *bursts,int nBursts,nb_rfi_struct *nb_rfi,int num_nb_rfi);
void loadInputFile(char *fname,simDataStruct *simData,burstStruct *bursts,int *nBursts,nb_rfi_struct *nb_rfi,int *num_nb_rfi);
void deleteUnwantedTables(fitsfile *fp);
void createPrimaryHeader(fitsfile *fp,simDataStruct *simData);
void createSubintHeader(fitsfile *fp,simDataStruct *simData);
void createData(fitsfile *fp,simDataStruct *simData,burstStruct *bursts,int nBursts,nb_rfi_struct *nb_rfi,int num_nb_rfi);
void digitise(float *data,unsigned char *cdata,simDataStruct *simData);

int main(int argc,char *argv[])
{
  int i;
  char fname[128];
  simDataStruct simData;
  burstStruct *bursts;
  nb_rfi_struct *nb_rfi;
  int nBursts=0;
  int num_nb_rfi = 0;
  
  for (i=0;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	strcpy(fname,argv[++i]);
    }

  bursts = (burstStruct *)malloc(sizeof(burstStruct)*MAX_BURSTS);
  nb_rfi = (nb_rfi_struct *)malloc(sizeof(nb_rfi_struct)*MAX_NB_RFI);
  
  loadInputFile(fname,&simData,bursts,&nBursts,nb_rfi,&num_nb_rfi);
  createSimulation(&simData,bursts,nBursts,nb_rfi,num_nb_rfi);

  free(bursts);
  free(nb_rfi);
}

void loadInputFile(char *fname,simDataStruct *simData,burstStruct *bursts,int *nBursts,nb_rfi_struct *nb_rfi,int *num_nb_rfi)
{
  FILE *fin;
  char line[2048];
  char first[128],second[128];
  
  fin = fopen(fname,"r");
  while (!feof(fin))
    {
      if (fgets(line,2048,fin)!=NULL)
	{
	  sscanf(line,"%s %s",first,second);
	  if (strcmp(first,"name:")==0)
	    strcpy(simData->fname,second);
	  else if (strcmp(first,"nbits:")==0)
	    sscanf(second,"%d",&(simData->nbits));
	  else if (strcmp(first,"npol:")==0)
	    sscanf(second,"%d",&(simData->npol));
	  else if (strcmp(first,"nchan:")==0)
	    sscanf(second,"%d",&(simData->nchan));
	  else if (strcmp(first,"nsblk:")==0)
	    sscanf(second,"%d",&(simData->nsblk));
	  else if (strcmp(first,"nsub:")==0)
	    sscanf(second,"%d",&(simData->nsub));
	  else if (strcmp(first,"tsamp:")==0)
	    sscanf(second,"%f",&(simData->tsamp));
	  else if (strcmp(first,"fchan0:")==0)
	    sscanf(second,"%f",&(simData->fchan0));
	  else if (strcmp(first,"chanbw:")==0)
	    sscanf(second,"%f",&(simData->chanbw));
	  else if (strcmp(first,"nb_rfi:")==0)
	    {
	      char *tok,*tok2,*tok3;
	      char line2[2048];
	      char *s1,*s2;
	      /* get the first token */
	      tok = strtok_r(line, " ",&s1);
	      
	      tok = strtok_r(NULL, ", ",&s1);
	      /* walk through other tokens */
	      while( s1 != NULL )
		{
		  strcpy(line2,tok);
		  printf( " %s >%s<\n", tok,line2 );
		  printf("---\n");
		  tok = strtok_r(NULL, " ,\n",&s1);
		  tok2 = strtok_r(line2,"=",&s2);
		  tok3 = strtok_r(NULL," ",&s2);
		  printf("tok2 = %s %s\n",tok2,tok3);
		  
		  if (strcmp(tok2,"t1")==0)
		    sscanf(tok3,"%f",&(nb_rfi[*num_nb_rfi].t1));
		  else if (strcmp(tok2,"t2")==0)
		    sscanf(tok3,"%f",&(nb_rfi[*num_nb_rfi].t2));
		  else if (strcmp(tok2,"f1")==0)
		    sscanf(tok3,"%f",&(nb_rfi[*num_nb_rfi].f1));
		  else if (strcmp(tok2,"f2")==0)
		    sscanf(tok3,"%f",&(nb_rfi[*num_nb_rfi].f2));
		  else if (strcmp(tok2,"s")==0)
		    sscanf(tok3,"%f",&(nb_rfi[*num_nb_rfi].flux));
		  else
		    {printf("ERROR: Unknown command: %s\n",tok2); exit(1);}
		}
	      (*num_nb_rfi)++;
	    }	  
	  else if (strcmp(first,"burst:")==0)
	    {
	      char *tok,*tok2,*tok3;
	      char line2[2048];
	      char *s1,*s2;
	      /* get the first token */
	      tok = strtok_r(line, " ",&s1);
	      
	      tok = strtok_r(NULL, ", ",&s1);
	      printf("At this point %s\n",tok);
	      /* walk through other tokens */
	      while( tok != NULL)
		{
		  printf("Again >%s<\n",line2);
		  strcpy(line2,tok);
		  printf("This point\n");
		  printf( " %s >%s<\n", tok,line2 );
		  printf("---\n");
		  tok = strtok_r(NULL, " ,\n",&s1);
		  printf("Tok here is >%s<\n",tok);
		  
		  tok2 = strtok_r(line2,"=",&s2);
		  tok3 = strtok_r(NULL," ",&s2);
		  printf("tok2 = %s %s\n",tok2,tok3);
		  if (strcmp(tok2,"dm")==0)
		    sscanf(tok3,"%f",&(bursts[*nBursts].dm));
		  else if (strcmp(tok2,"t")==0)
		    sscanf(tok3,"%f",&(bursts[*nBursts].time));
 		  else if (strcmp(tok2,"s")==0)
		    sscanf(tok3,"%f",&(bursts[*nBursts].flux));
 		  else if (strcmp(tok2,"w")==0)
		    sscanf(tok3,"%f",&(bursts[*nBursts].width));
		  else
		    {printf("ERROR: Unknown command: %s\n",tok2); exit(1);}
		  printf("Got here\n");
		}
	      printf("Leaving\n");
	      (*nBursts)++;
	    }
	}
    }
  fclose(fin);
  

}


void createSimulation(simDataStruct *simData,burstStruct *bursts,int nBursts,nb_rfi_struct *nb_rfi,int num_nb_rfi)
{
  fitsfile *fp;
  char fname[1024];
  int status=0;
  printf("In here with %s\n",simData->fname);
  sprintf(fname,"!%s(psrheader.fits)",simData->fname);
  fits_create_file(&fp,fname,&status);
  if (status)
    {
      fits_report_error(stderr,status);
      exit(1);
    }

  // Create primary header
  createPrimaryHeader(fp,simData);
  
  // Delete unwanted tables
  deleteUnwantedTables(fp);

  // Create subintegration header
  createSubintHeader(fp,simData);

  // Now simulate and digitise the data
  createData(fp,simData,bursts,nBursts,nb_rfi,num_nb_rfi);
    
  fits_close_file(fp,&status);

  
}

void createData(fitsfile *fp,simDataStruct *simData,burstStruct *bursts,int nBursts,nb_rfi_struct *nb_rfi,int num_nb_rfi)
{
  float *data;
  unsigned char *cdata;
  int sub;
  int i,j,p,b;
  float val,tval;
  long seed = TKsetSeed();
  int status=0;
  int colnum;
  int samplesperbyte = 8/simData->nbits;
  long naxes[4];
  int ival;
  char name[128];
  float datFreq[simData->nchan];
  float datWts[simData->nchan];
  float datOffs[simData->nchan];
  float datScl[simData->nchan];

  float fref = simData->fchan0+simData->chanbw*simData->nchan/2.0;
  float toff;
  
  
  strcpy(name,"SUBINT");
  fits_movnam_hdu(fp,BINARY_TBL,name,0,&status);fits_report_error(stderr,status);
  
  data = (float *)malloc(sizeof(float)*simData->npol*simData->nchan*simData->nsblk);
  cdata = (unsigned char *)malloc(sizeof(unsigned char)*simData->npol*simData->nchan*simData->nsblk);

  naxes[0] = simData->nchan;
  naxes[1] = simData->npol;
  naxes[2] = simData->nsblk;
  
  for (i=0;i<simData->nchan;i++)
    {
      datFreq[i] = simData->fchan0+simData->chanbw*i;
      datWts[i] = datOffs[i] = datScl[i] = 1.0;
    }
  
  for (sub=0;sub<simData->nsub;sub++)
    {
      printf("Processing subint: %d/%d\n",sub+1,simData->nsub);
      for (i=0;i<simData->nsblk;i++)
	{
	  tval = (sub*simData->nsblk+i)*simData->tsamp;
	  for (p=0;p<simData->npol;p++)
	    {
	      for (j=0;j<simData->nchan;j++)
		{
		  val = TKgaussDev(&seed);
		  for (b=0;b<nBursts;b++)
		    {
		      toff = 4.15e-3*bursts[b].dm*(pow(fref/1000.0,-2)-pow((simData->fchan0+simData->chanbw*j)/1000.0,-2));
		      if (bursts[b].time-bursts[b].width/2.0 <= tval+toff && bursts[b].time+bursts[b].width/2.0 > tval+toff)
			{
			  //			  printf("val = %g\n",val);
			  val+=bursts[b].flux;
			}
		    }
		  for (b=0;b<num_nb_rfi;b++)
		    {
		      if (datFreq[j]>=nb_rfi[b].f1 && datFreq[j]<nb_rfi[b].f2 &&
			  (tval >=nb_rfi[b].t1 && tval < nb_rfi[b].t2))
			val+=nb_rfi[b].flux;
		    }
		  data[i*simData->npol*simData->nchan+p*simData->nchan+j]=val;
		  
		}
	    }
	}
      digitise(data,cdata,simData);

      fits_insert_rows(fp,sub,1,&status);
      if (status) { fits_report_error(stderr,status); exit(1);}
      fits_get_colnum(fp,CASEINSEN,(char *)"INDEXVAL",&colnum,&status); fits_report_error(stderr,status);
      ival = sub+1;  fits_write_col(fp,TINT,colnum,sub+1,1,1,&ival,&status);

      fits_get_colnum(fp,CASEINSEN,(char *)"DAT_FREQ",&colnum,&status); fits_report_error(stderr,status);
      fits_modify_vector_len(fp,colnum,simData->nchan,&status); fits_report_error(stderr,status);
      fits_write_col(fp,TFLOAT,colnum,sub+1,1,simData->nchan,datFreq,&status);

      fits_get_colnum(fp,CASEINSEN,(char *)"DAT_WTS",&colnum,&status); fits_report_error(stderr,status);
      fits_modify_vector_len(fp,colnum,simData->nchan,&status); fits_report_error(stderr,status);
      fits_write_col(fp,TFLOAT,colnum,sub+1,1,simData->nchan,datWts,&status);

      fits_get_colnum(fp,CASEINSEN,(char *)"DAT_OFFS",&colnum,&status); fits_report_error(stderr,status);
      fits_modify_vector_len(fp,colnum,simData->nchan,&status); fits_report_error(stderr,status);
      fits_write_col(fp,TFLOAT,colnum,sub+1,1,simData->nchan,datOffs,&status);

      fits_get_colnum(fp,CASEINSEN,(char *)"DAT_SCL",&colnum,&status); fits_report_error(stderr,status);
      fits_modify_vector_len(fp,colnum,simData->nchan,&status); fits_report_error(stderr,status);
      fits_write_col(fp,TFLOAT,colnum,sub+1,1,simData->nchan,datScl,&status);

      
      // Write the data
      if (status) { fits_report_error(stderr,status); exit(1);}
      printf("writing the data %d %d %d %d %d\n",samplesperbyte,simData->nsblk,simData->npol,simData->nchan,simData->nsblk*simData->npol*simData->nchan/samplesperbyte);
      fits_get_colnum(fp,CASEINSEN,(char *)"DATA",&colnum,&status); fits_report_error(stderr,status); 
      fits_modify_vector_len (fp, colnum, (simData->nchan*simData->npol*simData->nsblk), &status); fits_report_error(stderr,status);
      fits_delete_key(fp, (char *)"TDIM18", &status); // THIS SHOULD NOT BE HARDCODED
      fits_write_tdim(fp,colnum,3,naxes,&status);fits_report_error(stderr,status);
	      
      fits_write_col_byt(fp,colnum,sub+1,1,simData->nsblk*simData->npol*simData->nchan/samplesperbyte,cdata,&status);
      if (status) { fits_report_error(stderr,status); exit(1);}
      printf("Done writing\n"); 
    }
  
  
  free(data);
  free(cdata);
}

void digitise(float *data,unsigned char *cdata,simDataStruct *simData)
{
  int i,j;
  unsigned char tc;
  double bit_level=0;
  
  if (simData->nbits==1)
    {
      long n=0;
      for (i=0;i<simData->nchan*simData->nsblk*simData->npol/8;i++)
	{
	  tc=0;
	  for (j=0;j<8;j++)
	    {
	      if (data[n] > bit_level)
		tc = tc | (1 << (7-j));
	      n++;
	    }
	  cdata[i] = tc;
	}
    }
  else
    {
      printf("Error: do not know how to process this number of bits\n");
    }
}

void deleteUnwantedTables(fitsfile *fp)
{
  int status=0;
  int newHDUtype;
  
  fits_movnam_hdu(fp, BINARY_TBL, (char *)"BANDPASS", 0, &status);  fits_delete_hdu(fp, &newHDUtype, &status);   fits_report_error(stdout,status);
  fits_movnam_hdu(fp, BINARY_TBL, (char *)"COHDDISP", 0, &status);  fits_delete_hdu(fp, &newHDUtype, &status);   fits_report_error(stdout,status);
  fits_movnam_hdu(fp, BINARY_TBL, (char *)"PSRPARAM", 0, &status);  fits_delete_hdu(fp, &newHDUtype, &status);   fits_report_error(stdout,status);
  fits_movnam_hdu(fp, BINARY_TBL, (char *)"POLYCO", 0, &status);  fits_delete_hdu(fp, &newHDUtype, &status);   fits_report_error(stdout,status);
  fits_movnam_hdu(fp, BINARY_TBL, (char *)"T2PREDICT", 0, &status);  fits_delete_hdu(fp, &newHDUtype, &status);   fits_report_error(stdout,status);
  fits_movnam_hdu(fp, BINARY_TBL, (char *)"FLUX_CAL", 0, &status);  fits_delete_hdu(fp, &newHDUtype, &status);   fits_report_error(stdout,status);
  fits_movnam_hdu(fp, BINARY_TBL, (char *)"CAL_POLN", 0, &status);  fits_delete_hdu(fp, &newHDUtype, &status);   fits_report_error(stdout,status);
  fits_movnam_hdu(fp, BINARY_TBL, (char *)"FEEDPAR", 0, &status);  fits_delete_hdu(fp, &newHDUtype, &status);   fits_report_error(stdout,status);

}

void createPrimaryHeader(fitsfile *fp,simDataStruct *simData)
{
  int status=0;
  char str[128];
  float freqc = simData->fchan0+simData->chanbw*simData->nchan/2.0;
  float bw = fabs(simData->chanbw*simData->nchan);
  int stt_imjd=55000;
  float stt_smjd=123;
  float stt_offs=0.0;
  
  fits_movabs_hdu(fp,1,NULL,&status);
  fits_write_date(fp,&status);
  
  fits_update_key(fp,TFLOAT,"OBSFREQ",&freqc,NULL,&status);
  fits_update_key(fp,TFLOAT,"OBSBW",&bw,NULL,&status);
  fits_update_key(fp,TINT,"OBSNCHAN",&(simData->nchan),NULL,&status);
  fits_update_key(fp,TINT,"STT_IMJD",&stt_imjd,NULL,&status);
  fits_update_key(fp,TFLOAT,"STT_SMJD",&stt_smjd,NULL,&status);
  fits_update_key(fp,TFLOAT,"STT_OFFS",&stt_offs,NULL,&status);
}

void createSubintHeader(fitsfile *fp,simDataStruct *simData)
{
  int status=0;
  char name[128];
  char str[128];
  char cval[128];

  strcpy(name,"SUBINT");
  fits_movnam_hdu(fp,BINARY_TBL,name,0,&status);fits_report_error(stderr,status);

  if (simData->nbits==1)
    {
      sprintf(cval,"16X");//,nchan*nsblk);
      fits_update_key(fp, TSTRING, (char *)"TFORM18", cval, NULL, &status); // Was 20
      fits_report_error(stderr,status);
    }

    fits_update_key_str(fp,(char *)"INT_TYPE",(char *)"TIME",(char *)"",&status);  fits_report_error(stderr,status);
  fits_update_key_str(fp,(char *)"INT_UNIT",(char *)"SEC",(char *)"",&status);  fits_report_error(stderr,status);
  fits_update_key_str(fp,(char *)"SCALE",(char *)"FluxDen",(char *)"",&status);  fits_report_error(stderr,status);
  fits_update_key_str(fp,(char *)"NBIN",(char *)"1",(char *)"",&status);  fits_report_error(stderr,status);
  fits_update_key_str(fp,(char *)"ZERO_OFF",(char *)"0",(char *)"",&status);  fits_report_error(stderr,status);
  fits_update_key_str(fp,(char *)"SIGNINT",(char *)"0",(char *)"",&status);  fits_report_error(stderr,status);
  fits_update_key_str(fp,(char *)"NSUBOFFS",(char *)"0",(char *)"",&status);  fits_report_error(stderr,status);


  fits_update_key(fp,TINT,"NPOL",&(simData->npol),NULL,&status);
  fits_update_key(fp,TINT,"NBITS",&(simData->nbits),NULL,&status);
  fits_update_key(fp,TINT,"NSBLK",&(simData->nsblk),NULL,&status);
  fits_update_key(fp,TINT,"NCHAN",&(simData->nchan),NULL,&status);
  fits_update_key(fp,TFLOAT,"TBIN",&(simData->tsamp),NULL,&status);
  fits_update_key(fp,TFLOAT,"CHAN_BW",&(simData->chanbw),NULL,&status);

}
