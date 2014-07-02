/* this file should allow us to compute the value of 
** sigmaw(y) = - sum w_i y_i + \int \exp \hbary(x) dx
** with nonnegative weights summing to one
*/

#include <math.h>
#include <R.h>
#include <Rmath.h>
#include "logcon.h"

double sigmaeffw( double y[] ) {
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
  double integral = 0.0;
  int totaldim, totalpoints, nf;
  double *A;
  double *ytmp;
  int inhull;
 
  /* Setup */
  allpoints = Calloc( ( npoints ) * ( dim + 1 ), double );
  A = Calloc( ( dim ) * ( dim ), double );
  ytmp = Calloc( (dim + 1), double );
  yminimum = ymin( y, truepoints );
  ymaximum = ymax( y, truepoints );
  totaldim = dim + 1; 
  yrange = ymaximum - yminimum;

  /* What convhull do we want? */
  for ( i = 0; i < truepoints; i++ ) {
    for ( j = 0; j < dim; j++ ) {
      allpoints[ totaldim * i + j] = xdata[ i + j * truepoints ];    
    }
    allpoints[ totaldim * i + dim ] = y[ i ] / yrange;  
  }
  for ( i = truepoints; i < npoints; i++ ) {
    for ( j = 0; j < dim; j++ ) {
      allpoints[ totaldim * i + j ] = xdata[ ( i - truepoints ) + j * truepoints ];
    }
    allpoints[ totaldim * i + dim ] = yminimum / yrange - 0.1;
  }
  
  /* NOTE: use convhullnmlc.c so that the index does not need to be shifted by 1 */
  totalpoints = npoints;
  outpoints = convhullnmlc( allpoints, &totalpoints, &totaldim, &nf, chopts );
  
  for ( i = 0; i < nf; i++ ) {
    inhull = 0; 
    for ( j = 0; j < totaldim; j++ ) { 
      inhull += ( outpoints[ i + nf * j ] >= truepoints );
    }
    /* Remove from the list all facets which don't give a contribution */
    if ( inhull == 0 )  { 

      /* calculate the contribution to the integral */
      
      /* Find the relevant A, ytmp */
      for ( j = 1; j <= dim; j++ ) {
	for ( k = 0; k < dim; k++ ) {
	  //	  A[ ( j - 1 ) + k * dim ] = allpoints[ ( outpoints[ i + nf * j ]) * totaldim + k ] - allpoints[ ( outpoints[ i ] ) * totaldim + k ];
	  A[ (j - 1 ) + k * dim ] = xdata[ ( outpoints[ i + nf * j ] ) + k * truepoints ] - xdata[ outpoints[ i ] + k * truepoints ];
	}
      }
      for (j=0; j<=dim; j++) {
	//ytmp[ j ] = allpoints[ ( outpoints[ i + nf * j ] ) * totaldim + dim ];
	ytmp[ j ] = y[ outpoints[ i + nf * j ] ];
      }
      
      /* Find absdetA */
      /*Add this contribution to the integral */
      integral += absdet( A, dim, 0 ) * JAD( ytmp, dim, Jtol );     
    }
  }
  
  /* Free everything */
  Free(allpoints);
  Free(A);
  Free(outpoints);
  Free(ytmp);
  return( integral - dotprod(y,weights,truepoints));
}
