/* Computes \argmax_{y \in \R^n} -\sum w_i y_i + \int \exp \hbary(x) dx */
/* for nonnegative weights summing to one             */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>  
#include <R.h>
#include <Rmath.h>

extern int npoints;
extern int dim;
extern double *xdata;
double *weights;
extern int nouter;
extern int truepoints;
extern double Jtol;

/* Function to be minimized: */
double sigmaeffw(double *y);
double dnull_entry();

/* Subgradient: */
void subgradeffw(double y[], double g[]);
void null_entry();

double solvoptweights(int npoints, double *y_in, double sigma_ralg2(),void subgrad_ralg2(),double *opt_out, double fun() ,void fun2(), double *parameters );

void renormalise(double *y);

void logconestw (double *y_in, double *xdata_in, int *d_in, int *n_in, double *weights_in, double *opt_out, double *sigmavalue_out, double *parameters, double *Jtol_in, int *nouter_in)
{

  /* Initialise */
  truepoints = *n_in;
  dim = *d_in;
  xdata = xdata_in; 
  nouter = *nouter_in;
  npoints = truepoints + nouter;
  weights = weights_in;
  Jtol = *Jtol_in;

  /* Use the solvoptweights */
    *sigmavalue_out=solvoptweights(truepoints,y_in,&sigmaeffw,&subgradeffw,opt_out,&dnull_entry,&null_entry,parameters);
 
    renormalise(y_in);
    /*That's all!*/
}
