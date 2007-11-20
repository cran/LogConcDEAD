#include <R.h>
#include <Rmath.h>

#define dgetrf dgetrf_
void dgetrf(int*, int*, double*, int*, int*, int*);

/*
 * determinant:
 * 
 * returns the log determinant of the n x n
 * replacing it with its LU decomposition
 */


double absdet( double *M, int n, int useLog) 
{
  double det, modulus;
  int i, info;
  int *p;

  p = (int *) malloc(sizeof(int) * n);
  
  /* LU decopmpose M */
  dgetrf(&n, &n, M, &n, p, &info);
  if(info != 0) {
#ifdef DEBUG
    warning("bad chol decomp in log_determinant");
#endif
    return -1e300*1e300;
  }  

 
  /* copied from R source and removing all reference to sign */
  if (useLog) {
    modulus = 0.0;
    for (i = 0; i < n; i++) {
      double dii = M[i*(n + 1)]; /* ith diagonal element */
      modulus += log(dii < 0 ? -dii : dii);
      }
    det = exp(modulus);
  } else {
    modulus = 1.0;
    for (i = 0; i < n; i++)
      modulus *= M[i*(n + 1)];
    if (modulus < 0) {
      modulus = -modulus;
          }
    det = modulus;
  }

  free(p);

  return det;
}






double determinant(M, n, useLog)
int n, useLog;
double *M;
{
  double det, modulus;
  int i, info, sign;
  int *p;

  p = (int *) malloc(sizeof(int) * n);
  
  /* LU decopmpose M */
  dgetrf(&n, &n, M, &n, p, &info);
  if(info != 0) {
#ifdef DEBUG
    warning("bad chol decomp in log_determinant");
#endif
    return -1e300*1e300;
  }  

 
  /* copied from R source to get the sign right */

  sign = 1;
  for (i = 0; i < n; i++) 
    if (p[i] != (i + 1))
      sign = -sign;
  if (useLog) {
    modulus = 0.0;
    for (i = 0; i < n; i++) {
      double dii = M[i*(n + 1)]; /* ith diagonal element */
      modulus += log(dii < 0 ? -dii : dii);
      if (dii < 0) sign = -sign;
    }
    det = sign * exp(modulus);
  } else {
    modulus = 1.0;
    for (i = 0; i < n; i++)
      modulus *= M[i*(n + 1)];
    if (modulus < 0) {
      modulus = -modulus;
      sign = -sign;
    }
    det = sign * modulus;
  }

  free(p);

  return det;
}


void det(double *M_in, int *n_in, int *useLog_in, double *det_out)
{
  *det_out = determinant(M_in, *n_in, *useLog_in);
}

void absdeterminant(double *M_in, int *n_in, int *useLog_in, double *det_out)
{ 
  *det_out = absdet(M_in, *n_in, *useLog_in); 
}


/* A few useful things, probably already available, but here we go: */
double ymin( double y[], int n) 
  {
    int i;
    double tmp = y[0];
    for (i=1; i<n; i++) 
      if (tmp > y[i]) tmp = y[i]; 
    return (tmp); 
  }

double mean(double *y, int len) 
{
  int i;
  double tmp = 0.0;
  for (i=0; i<len; i++)  tmp += y[i];
  tmp = tmp/((double)(len));
  return tmp; 
}



double max(double a,double b) {
	if(a>=b) return(a);
	else return(b);
}

double min(double a,double b) {
	if(a<=b) return(a);
	else return(b);
}

double dotprod(double y[], double w[], int n) {
  int i;
  double tmp = 0.0;
  for (i=0; i<n; i++) {
    tmp+= y[i]*w[i]; 
  }

  return tmp; }
