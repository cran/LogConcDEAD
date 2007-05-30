/* this file should allow us to compute the value of sigma */

#include <math.h>
#include <R.h>
#include <Rmath.h>

int *convhullnmlc(double *x_in, int *nrow_in, int* ncol_in, int *nf);

double sigmaeff(double y[]) 
  {
  extern int npoints;
  extern int dim;
  extern double *xdata;
  /* extern int nouter; */
  extern int truepoints;
  double ymin (double y[],int n);
  double min(double a, double b);
  double absdet(double *a, int n, int useLog);
  double mean(double y[], int len);
  double *allpoints;
  int i, j, k;
  double yminimum;
  int *outpoints; 
  double integral = 0.0;
  int totaldim, totalpoints, nf;
  double *zprod;
  double *A;
  double *z;
  double sum;
  int inhull;
  double absdetA;
  int *supporting;
  /* int count=0; */

  allpoints = Calloc((npoints)*(dim + 1),double);
  A = Calloc((dim)*(dim),double);
  z = Calloc(dim,double);
  zprod = Calloc(dim,double);
  supporting = Calloc(npoints,int);
  
   
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
  
    outpoints = convhullnmlc(allpoints, &totalpoints, &totaldim, &nf);
   
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
  
/* Find the relevant A, z */
      for (j=1; j<=dim; j++) 
        {
        for (k=0; k<dim; k++) 
          {
	  A[(j-1)+k*dim] = allpoints[(outpoints[i+nf*j])*totaldim + k] - allpoints[(outpoints[i])*totaldim + k];
	  }
        z[(j-1)] = allpoints[(outpoints[i+nf*j])*totaldim + dim] - allpoints[outpoints[i]*totaldim + dim];
	}

      /* Find absdetA */
      absdetA = absdet(A,dim,0);
      
      /* Find zprod */
      for (j=0; j<dim; j++) 
        {
        zprod[j] = 1.0;
        for (k=0; k<dim; k++) 
          if(k != j)  zprod[j] *= (z[j]-z[k]);
        }


      sum=0.0;
      /* Compute the strange sum which gives us our answer */
      for (j=0; j<dim; j++) 
        sum += (exp(z[j]) - 1)/(z[j]*zprod[j]);  

      /* Add this contribution to the integral */
      integral += absdetA*exp(allpoints[totaldim*(outpoints[i])+dim])*sum; 

      }
    }
  
   Free(allpoints);
   Free(A);
   Free(z);
   Free(zprod);
   Free(outpoints);
   Free(supporting);

   return(integral - mean(y,truepoints));
    }
