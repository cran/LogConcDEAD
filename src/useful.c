/* This contains various useful functions */
/*                                        */
/* J* are based on ideas of Lutz Duembgen */
/* at the University of Bern              */

#include <R.h>
#include <Rmath.h>
#include <stdlib.h>
#include <stdio.h>

#define dgetrf dgetrf_
void dgetrf(int*, int*, double*, int*, int*, int*);

/*
 * determinant:
 * 
 * returns the log determinant of the n x n
 * replacing it with its LU decomposition
 */

/* Compute the determinant */
double absdet(double *M, int n, int useLog) 
{
  double det, modulus;
  int i, info;
  int *p;

  p = (int *) malloc(sizeof(int) * n);
  
  /* LU decopmpose M */
  dgetrf(&n, &n, M, &n, p, &info);
  if(info != 0) {
    warning("bad chol decomp in log_determinant");
    return 0;
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
     /* Computes the determinant */
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
    return tmp; 
  }

double ymax( double y[], int n ) {
  int i;
  double tmp = y[0];
  for( i=1; i<n; i++ )
    if( tmp < y[i] ) tmp = y[i];
  return tmp;
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
  return tmp;
}

/* Here are the auxiliary functions for the computation of 
   J(y[0], ..., y[d]) and its derivatives */

int cmp_double (const void *X, const void *Y)
/* A compare function for use in the sorting stage */
{
  if (*((double *)X) > *((double *)Y)) {
    return 1;
  }
  else{
    if (*((double *)X) < *((double *)Y)){
      return -1;
    }
    else	   {
      return 0;
    }
  }
}

double JAD_appr(double *y, int d) {
  /* Uses a Taylor expansion about ymean to the third order */
  double ymean = mean(y, (d+1));
  int i;
  double tmp;
  double *z;
  z = malloc((d+1)*sizeof(double));
  for (i=0; i<=d; i++) {
    z[i] = y[i] - ymean;
  }
   tmp = 1.0;
  /*add the z^2 and z^3 terms*/
  for (i=0; i<=d; i++) {
    tmp += (z[i]*z[i]/(2*(d+1)*(d+2)) + z[i]*z[i]*z[i]/(3*(d+1)*(d+2)*(d+3)));
   }
  /* divide by d! */
  for (i=1; i<=d; i++) {
    tmp /= i;
  }
  tmp *= exp(ymean);
  free(z);
  return(tmp);
}

double JAD_ord(double *y, int d, double eps) {
  double tmp;
  double e1, e2;
  /* Computes the function J in d dimensions given ordered inputs */ 
 /* if y[d]-y[0] < eps, use taylor expansion */
  if(y[d]-y[0] < eps) {
    return(JAD_appr(y, d));
  }
  else {
    if(d==1) {
      tmp = (exp(y[0]) - exp(y[1]))/(y[0]-y[1]);
      return(tmp);
    }
    else {
      e1 = JAD_ord(&y[1], d-1, eps);
      e2 = JAD_ord(&y[0], d-1, eps);
      return((e1 - e2)/(y[d] - y[0])); 
    }  
  }
}

double JAD(double *y, int d, double eps) {
  /* First sort the y, then apply JAD_ord */
  double *z;
  int k;
  double tmp;
  z=Calloc( d+1, double);
  for (k=0; k<=d; k++) {
    z[k] = y[k];
  }
  qsort(z,d+1,sizeof(double),cmp_double);

  tmp = JAD_ord(z,d,eps);
  
  Free(z);
  return(tmp);
}

double JiAD(double *y, int i, int d, double eps) {
  /* Compute the ith partial derivative as
     J(y[0], ..., y[i], y[i], ..., y[d])*/
  double *z;
  int k;
  double tmp;
  z = malloc((d+2)*sizeof(double));
  for (k=0; k<=d; k++) {
    z[k] = y[k];
  }
  z[d+1] = y[i];
  
  qsort(z,d+2,sizeof(double),cmp_double);
  
  tmp = JAD_ord(z,d+1,eps);
  
  free(z);
  return(tmp);
}

double JijAD(double *y, int i, int j, int d, double eps) {
  /* Compute the second partial derivatives as
     J(y[0], ..., y[i], y[i], ...,y[j],y[j], ..., y[d]) */
  double *z;
  int k;
  double tmp;
  z = malloc((d+3)*sizeof(double));
  for (k=0; k<=d; k++) {
    z[k] = y[k];
  }
  z[d+1] = y[i];
  z[d+2] = y[j];
  
  qsort(z,d+3,sizeof(double),cmp_double);
  
  tmp = JAD_ord(z,d+2,eps);
  
  /* If we want second derivative we need to multiply by 2 */
  if (i == j) tmp *= 2;
  
  free(z);
  return(tmp);
}
