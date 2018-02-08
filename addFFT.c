#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

int main(int argc,char *argv[])
{
  FILE *fin;
  FILE *fout;
  int i;
  
  long maxVals = 600000;
  float *vals;
  float *freqs;
  long n=0;
  float xval,yval;
  
  if (!(vals = (float *)malloc(sizeof(float)*maxVals)))
    {
      printf("Unable to allocate memory\n");
      exit(1);
    }
  if (!(freqs = (float *)malloc(sizeof(float)*maxVals)))
    {
      printf("Unable to allocate memory\n");
      exit(1);
    }

  
  for (i=1;i<argc;i++)
    {
      fin = fopen(argv[i],"r");
      n=0;
      while (!feof(fin))
	{
	  if (fscanf(fin,"%f %f",&xval,&yval)==2)
	    {
	      if (i==1)
		{
		  freqs[n] = xval;
		  vals[n] = yval;
		}
	      else
		vals[n] += yval;
	      n++;
	    }
	}
      fclose(fin);
    }
  fout = fopen("fftAdd.dat","w");
  for (i=0;i<n;i++)
    fprintf(fout,"%f %g\n",freqs[i],vals[i]);
  fclose(fout);
  free(vals);
free(freqs);
}
