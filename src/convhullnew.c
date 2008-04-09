/* Copyright (C) 2000  Kai Habel
** Copyright R-version (c) 2005 Raoul Grasman
**
** This program is free software; you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation; either version 2 of the License, or
** (at your option) any later version.
**
** This program is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
**
** You should have received a copy of the GNU General Public License
** along with this program; if not, write to the Free Software
** Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
*/

/*
29. July 2000 - Kai Habel: first release
2002-04-22 Paul Kienzle
* Use warning(...) function rather than writing to cerr

23. May 2005 - Raoul Grasman: ported to R
* Changed the interface for R

01. August 2007 R. B. Gramacy
* modified for inclusion in the LogConcDEAD package 
* with silent output
*/

#include <R.h>
#include <Rdefines.h>
#include "qhull_a.h"

/*
DEFUN_DLD (convhulln, args, ,
"-*- texinfo -*-\n\
@deftypefn {Loadable Function} {@var{H} =} convhulln (@var{p}[, @var{opt}])\n\
Returns an index vector to the points of the enclosing convex hull.\n\
The input matrix of size [n, dim] contains n points of dimension dim.\n\n\
If a second optional argument is given, it must be a string containing\n\
extra options for the underlying qhull command.  (See the Qhull\n\
documentation for the available options.)\n\n\
@seealso{convhull, delaunayn}\n\
@end deftypefn")
*/
SEXP convhullnew(const SEXP p, const SEXP options)
{
	SEXP retval;
	int curlong, totlong, i, j;
	unsigned int dim, n;
	int exitcode; 
	boolT ismalloc;
	char flags[250];             /* option flags for qhull, see qh_opt.htm */
	char *opts;
	int *idx;
	double *pt_array;


	FILE *outfile = stdout;      /* output from qh_produce_output() use NULL to skip qh_produce_output() */
	FILE *errfile = stderr;      /* error messages from qhull code */

	retval = R_NilValue;

	if(!isString(options) || length(options) != 1){
		error("Second argument must be a single string.");
	}
	if(!isMatrix(p) || !isReal(p)){
		error("First argument should be a real matrix.");
	}

	i=LENGTH(STRING_ELT(options,0));
	opts = (char *) R_alloc( ((i>1)?i:1), sizeof(char) );
	strcpy(opts, " ");
	if(i>1) strcpy(opts, CHAR(STRING_ELT(options,0)));

	dim = ncols(p);
	n   = nrows(p);
	if(dim <= 0 || n <= 0){
		error("Invalid input matrix.");
	}

	j=0;
	pt_array = (double *) R_alloc(n*dim, sizeof(double)); 
	for(i=0; i < n; i++)
		for(j=0; j < dim; j++)
			pt_array[dim*i+j] = REAL(p)[i+n*j]; /* could have been pt_array = REAL(p) if p had been transposed */

	ismalloc = False;   /* True if qhull should free points in qh_freeqhull() or reallocation */

	/* hmm  lot's of options for qhull here */
	sprintf(flags,"qhull Qt Tcv %s",opts);
	/* outfile = NULL; this produces nonsensical output ... needs to be fixed to get convexhulln quiet */
	exitcode = qh_new_qhull (dim,n,pt_array,ismalloc,flags,outfile,errfile);
	/* If you want some debugging information replace the NULL
	// pointer with stdout
	*/

	if (!exitcode) {  /* 0 if no error from qhull */

		facetT *facet;                  /* set by FORALLfacets */
		vertexT *vertex, **vertexp;		/* set by FORALLfacets */
		unsigned int n = qh num_facets;

		PROTECT(retval = allocMatrix(INTSXP, n, dim));
		idx = (int *) R_alloc(n*dim,sizeof(int));

		qh_vertexneighbors();

		i=0;
		FORALLfacets {
			j=0;
			/*std::cout << "Current index " << i << "," << j << std::endl << std::flush;
			// qh_printfacet(stdout,facet);
			*/
			FOREACHvertex_ (facet->vertices) {
				/* qh_printvertex(stdout,vertex); */
				if (j >= dim)
					warning("extra vertex %d of facet %d = %d",
					j++,i,1+qh_pointid(vertex->point));
				else
					idx[i+n*j++] = 1 + qh_pointid(vertex->point);
			}
			if (j < dim) warning("facet %d only has %d vertices",i,j);
			i++;
		}
		j=0;
		for(i=0;i<nrows(retval);i++)
			for(j=0;j<ncols(retval);j++)
				INTEGER(retval)[i+nrows(retval)*j] = idx[i+n*j];
		UNPROTECT(1);
	}
	qh_freeqhull(!qh_ALL);					/*free long memory */
	qh_memfreeshort (&curlong, &totlong);	/* free short memory and memory allocator */

	if (curlong || totlong) {
		warning("convhulln: did not free %d bytes of long memory (%d pieces)",
			totlong, curlong);
	}
	return retval;
}
