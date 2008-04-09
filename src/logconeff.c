/* This function computes the log-concave MLE of X_1, ..., X_n */
/* using the "solvopt" implementation of Shor's r-algorithm */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>  
#include "solvopt.h"
/*#include "qhull.h"*/
#include <R.h>
#include <Rmath.h>

int npoints;
int dim;
double *xdata;
int nouter;
int truepoints;
double Jtol;
double trouble;

/*  Function to be minimized: */
double sigmaeff(double *y);
double dnull_entry();

/* Subgradient: */
void subgradeff(double y[], double g[]);
void null_entry();

/* SolvOpt routine, modified to take extra parameters as arguments: */
double solvopt2(int npoints, double *y_in, double sigma_ralg2(), void subgrad_ralg2(), double *opt_out, double fun(), void fun2(),double *parameters );

void renormalise(double *y);

void logconesteff (double *y_in, double *xdata_in, int *d_in, int *n_in, double *opt_out, double *sigmavalue_out, double *parameters, double *Jtol_in, int *nouter_in)
{
  /* Initialise */
  truepoints = *n_in;
  dim = *d_in;
  xdata = xdata_in; 
  nouter = *nouter_in;
  npoints = truepoints+nouter;
  Jtol = *Jtol_in;

  /* Use the solvopt algorithm */
    *sigmavalue_out=solvopt2(truepoints,y_in,&sigmaeff,&subgradeff,opt_out,&dnull_entry,&null_entry,parameters);

  /* That's all! */

  renormalise(y_in);

}
