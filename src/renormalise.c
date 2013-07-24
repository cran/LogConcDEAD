/* this file should allow us to compute the value of sigma */

#include <math.h>
#include <R.h>
#include <Rmath.h>

int *convhullnmlc(double *x_in, int *nrow_in, int* ncol_in, int *nf, char *opts);

void renormalise(double y[]) 
  {
  extern int npoints;
  extern int dim;
  extern double *xdata;
  extern int truepoints;
  extern double Jtol;
  extern char* chopts;
  double ymin (double y[],int n);
  double min(double a, double b);
  double absdet(double *a, int n, int useLog);
  double JAD(double *y, int dim, double eps);
  double *allpoints;
  int i, j, k;
  double yminimum;
  int *outpoints; 
  double integral = 0.0;
  int totaldim, totalpoints, nf;
  double *A;
  double *ytmp;
  double tmp;
  int inhull;
  double absdetA;
  double logintegral;

  allpoints = Calloc((npoints)*(dim + 1),double);
  A = Calloc((dim)*(dim),double);
  ytmp = Calloc(dim+1,double);
       
  yminimum = ymin(y,truepoints) - 1.0; 
  totaldim=dim+1; 
 
  for (i=0; i<truepoints; i++) 
    {
    for (j=0; j<dim; j++) 
      {

      allpoints[totaldim*i + j] = xdata[i + j*truepoints];    
  
      }
    allpoints[totaldim*i + dim] = y[i];  
    }
  for (i=truepoints; i<npoints; i++)
    {
      for (j=0; j<dim; j++)
	{
	  allpoints[totaldim*i + j] = xdata[(i-truepoints)+j*truepoints];
	}
      allpoints[totaldim*i + dim] = yminimum;
    }
    
    totalpoints = npoints;

    outpoints = convhullnmlc(allpoints, &totalpoints, &totaldim, &nf, chopts);
   
 /* NOTE: use convhullnmlc.c so that the index does not need to be shifted by 1 */

  for (i=0; i<nf; i++) 
    {
    inhull = 0; 
    for (j=0; j<totaldim; j++) { inhull += (outpoints[i+nf*j]>=truepoints);
  
 }

    /* Remove from the list all facets which don't give a contribution */
    if (inhull==0) 
      { 

/* calculate the contribution to the integral */
  
/* Find the relevant A, ytmp */
      for (j=1; j<=dim; j++) 
        {
        for (k=0; k<dim; k++) 
          {
	  A[(j-1)+k*dim] = allpoints[(outpoints[i+nf*j])*totaldim + k] - allpoints[(outpoints[i])*totaldim + k];
	  }
	}
      for (j=0; j<=dim; j++) {
        ytmp[j] = allpoints[(outpoints[i+nf*j])*totaldim + dim];
	}

      /* Find absdetA */
      absdetA = absdet(A,dim,0);
      
      /* Find JAD */
      tmp = JAD(ytmp, dim, Jtol);

      /* Add this contribution to the integral */
      integral += absdetA*tmp; 
    
      }
    }

   Free(allpoints);
   Free(A);
   Free(ytmp);
   Free(outpoints);
  
   logintegral = log(integral);
   for (i=0; i<truepoints; i++) y[i]-=logintegral ;
  }
