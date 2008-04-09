/* this file should allow us to compute the subgradient with weights
** 02-apr-08 added numerically stable way to compute J
** due to Lutz Duembgen */

#include <math.h>
#include <R.h>
#include <Rmath.h>

int *convhullnmlc(double *x_in, int *nrow_in, int* ncol_in, int *nf);

void subgradeffw(double y[], double g[])
{
  extern int npoints;
  extern int dim;
  extern double *xdata;
  extern int truepoints;
  extern double Jtol;
  extern double *weights;
  double ymin (double y[],int n);
  double min(double a, double b);
  double absdet(double *a, int n, int log);
  double JiAD(double *y, int i, int dim, double eps);
  double *allpoints;
  int i, j, k;
  double yminimum;
  int *outpoints; 
  int totaldim, totalpoints, nf;
  double *A;
  double *ytmp;
  double Ji;
  int inhull;
  double absdeta;

  allpoints = Calloc((npoints)*(dim + 1),double);
  A = Calloc((dim)*(dim),double);
  ytmp = Calloc((dim+1),double);

  /* initialise the subgradient vector     */
  for (j=0; j<truepoints;j++) 
    g[j] = -weights[j]; 

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
      if (inhull==0) /* i.e. if a point on the surface */
	{
 
	  /* calculate the contribution to the integral */
    
	  /* First find the relevant A, ytmp */
	  for (j=1; j<=dim; j++) 
	    {
	      for (k=0; k<dim; k++)
		{
		  A[(j-1)+k*dim] = allpoints[outpoints[i+nf*j]*totaldim + k] - allpoints[outpoints[i]*totaldim + k];
		}
	    }
	  for (j=0; j<=dim; j++) {
	    ytmp[j] = allpoints[(outpoints[i+nf*j])*totaldim + dim];
	  }
	  /* Find the absolute value of det A */
	  absdeta = absdet(A,dim,0); 
   
	  /* Now we'll add the relavant parts on to the subgradient */
   	  for (j=0; j<=dim;j++) { 
	    Ji = JiAD(ytmp,j,dim,Jtol);
	    g[outpoints[i+nf*j]] += absdeta*Ji;
	  }
	}
    }
  /* Free the allocated memory*/
  Free(allpoints);
  Free(A);
  Free(ytmp);
  Free(outpoints);
}
