/* this file should allow us to compute the subgradient */

#include <math.h>
#include <R.h>
#include <Rmath.h>

int *convhullnmlc(double *x_in, int *nrow_in, int* ncol_in, int *nf);

void subgradeff(double y[], double g[])
{
  extern int npoints;
  extern int dim;
  extern double *xdata;
  /*  extern int nouter;*/
  extern int truepoints;
  double ymin (double y[],int n);
  double min(double a, double b);
  double absdet(double *a, int n, int log);
  double *allpoints;
  int i, j, k;
  double yminimum;
  int *outpoints; 
  /* double integral; */
  int totaldim, totalpoints, nf;
  double *zprod; 
  double prod;
  double *A;
  double *z;
  double sum, tot;
  int inhull;
  double absdeta;

  allpoints = Calloc((npoints)*(dim + 1),double);
  A = Calloc((dim)*(dim),double);
  z = Calloc(dim,double);
  zprod = Calloc(dim,double);

  /* initialise the subgradient vector     */
  for (j=0; j<truepoints;j++) 
    g[j] = (-1.0/(truepoints)); 

  yminimum = ymin(y,truepoints) - 1.0; 
  totaldim=dim+1; 

  /* just using the data points */
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
  
  /* Find the convex hull! */
 
  outpoints = convhullnmlc(allpoints, &totalpoints, &totaldim, &nf);

  /* For each facet of the convex hull, find if it is relevant */
 
for (i=0; i<nf; i++) 
    {
    inhull = 0; 
    for (j=0; j<totaldim; j++) inhull += (outpoints[i+nf*j]>=truepoints);
   
    if (inhull==0) /* i.e. if a face in relevant part of convhull */
      {
 
   /* calculate the contribution to the integral */
    
   /* First find the relevant A, z */
    for (j=1; j<=dim; j++) 
      {
      for (k=0; k<dim; k++)
        {
	A[(j-1)+k*dim] = allpoints[outpoints[i+nf*j]*totaldim + k] - allpoints[outpoints[i]*totaldim + k];
        }
        z[(j-1)] = allpoints[outpoints[i+nf*j]*totaldim + dim] - allpoints[outpoints[i]*totaldim + dim];
      }

    /* Find the absolute value of det A */
    absdeta = absdet(A,dim,0); 

    /* Find the zprod_j = Product_{k!=j} (z_j - z_k)  */
    for (j=0; j<dim; j++) 
      {
      zprod[j] = 1.0;
      for (k=0; k<dim; k++) 
        if(k != j) zprod[j] *= (z[j]-z[k]);
      }

    /* Find the product of all the z[j] */
    prod = 1.0;
    for (j=0;j<dim;j++)  prod *= z[j];
  
    /* Now we'll add the relavant parts on to the subgradient */
    tot = 0.0;
    for (j=0; j<dim;j++) 
      {
      sum = 0.0;
      for (k=0; k<dim; k++) 
        {
        if (k != j) 
          {
	  sum += exp(z[k])/(z[k]*(z[k]-z[j])*zprod[k]) - exp(z[j])/(z[k]*(z[k]-z[j])*zprod[k]); 
          } 
        }
      sum += pow(-1,dim)*(exp(z[j])-1)/(z[j]*prod) + exp(z[j])/(z[j]*zprod[j]); 
        
      g[outpoints[i+nf*(j+1)]] += absdeta*exp(allpoints[totaldim*(outpoints[i])+dim])*sum; 
      tot += sum;
      }

    /* Find the integral for the first vertex  */
    sum = 0.0;
    for(k=0;k<dim;k++) 
      sum += (exp(z[k]) - 1)/(z[k]*zprod[k]);

    g[outpoints[i]] += absdeta*exp(allpoints[totaldim*(outpoints[i])+dim])*(sum-tot);
      }
    }
  /* Free the allocated memory*/
    Free(allpoints);
    Free(A);
    Free(z);
    Free(zprod);
    Free(outpoints);
}
