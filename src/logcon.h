/*  A header file */

double ymin( double y[],int n );
double ymax( double y[], int n );
double min( double a, double b );
double absdet( double *a, int n, int useLog );
double mean( double y[], int len );
double dotprod( double y[], double weights[], int n );
double JAD( double *y, int dim, double eps );
double JiAD(double *y, int i, int d, double eps );

int *convhullnmlc(double *x_in, int *nrow_in, int *ncol_in, int *nf,
		  char *opts );

