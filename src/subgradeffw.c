/* this file should allow us to compute the subgradient with weights
** 02-apr-08 added numerically stable way to compute J
** due to Lutz Duembgen */

#include <math.h>
#include <R.h>
#include <Rmath.h>
#include "logcon.h"

void subgradeffw(double y[], double g[])
{
  extern int npoints;
  extern int dim;
  extern double *xdata;
  extern int truepoints;
  extern double Jtol;
  extern char* chopts;
  extern double *weights;
  double *allpoints;
  int i, j, k;
  double yminimum, ymaximum, yrange;
  int *outpoints; 
  int totaldim, totalpoints, nf;
  double *A;
  double *ytmp;
  int inhull;
  double absdetA;

  allpoints = Calloc((npoints)*(dim + 1),double);
  A = Calloc((dim)*(dim),double);
  ytmp = Calloc((dim+1),double);

  /* initialise the subgradient vector     */
  for (j=0; j<truepoints;j++) 
    g[j] = -weights[j]; 

  yminimum = ymin( y, truepoints ); 
  ymaximum = ymax( y, truepoints );
  ymaximum = totaldim = dim + 1; 
  yrange = ymaximum - yminimum;
  /* just using the data points */
  for (i=0; i<truepoints; i++)  {
      for (j=0; j<dim; j++)  {
	allpoints[totaldim*i + j] = xdata[ i + j * truepoints];
      }
      allpoints[totaldim*i + dim] = y[ i ] / yrange;
  }
  for (i=truepoints; i<npoints; i++) {
    for (j=0; j<dim; j++) {
      allpoints[totaldim*i + j] = xdata[(i-truepoints)+j*truepoints];
    }
    allpoints[totaldim*i + dim] = yminimum / yrange - 0.1;
  }
  
  totalpoints = npoints;
  
  /* Find the convex hull! */
  outpoints = convhullnmlc(allpoints, &totalpoints, &totaldim, &nf, chopts);
  
  /* For each facet of the convex hull, find if it is relevant */
  for (i=0; i<nf; i++) 
    {
      inhull = 0; 
      for (j=0; j<totaldim; j++) inhull += (outpoints[i+nf*j]>=truepoints);
      
      if (inhull==0) /* i.e. if a point on the surface */
	{
	  /* calculate the contribution to the integral */
	  
	  /* First find the relevant A, ytmp */
	  for (j=1; j<=dim; j++)  {
	    for (k=0; k<dim; k++) {
	      //  A[(j-1)+k*dim] = allpoints[outpoints[i+nf*j]*totaldim + k] - allpoints[outpoints[i]*totaldim + k];
	      A[ (j - 1 ) + k * dim ] = xdata[ ( outpoints[ i + nf * j ] ) + k * truepoints ] - xdata[ outpoints[ i ] + k * truepoints ];
	    }
	  }
	  for (j=0; j<=dim; j++) {
	    //  ytmp[j] = allpoints[(outpoints[i+nf*j])*totaldim + dim];
	    ytmp[ j ] = y[ outpoints[ i + nf * j ] ];
	  }
	  /* Find the absolute value of det A */
	  absdetA = absdet(A,dim,0); 
	  
	  /* Now we'll add the relavant parts on to the subgradient */
   	  for (j=0; j<=dim;j++) { 
	    g[outpoints[i+nf*j]] += absdetA*JiAD( ytmp, j, dim, Jtol );
	  }
	}
    }
  /* Free the allocated memory*/
  Free(allpoints);
  Free(A);
  Free(ytmp);
  Free(outpoints);
}
