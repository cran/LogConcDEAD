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


/* going to try something different, namely returning the pointer to the array. */

int *convhullnmlc(double *x_in, int *nrow_in, int* ncol_in, int *nf)
{
  int curlong, totlong, i, j;
  unsigned int dim, n;
  int exitcode; 
  boolT ismalloc;
  char flags[250];             /* option flags for qhull, see qh_opt.htm */
  /* char *opts; */
  int *idx;
  double *pt_array;

  FILE *outfile = stdout;      
  FILE *errfile = stderr;      /* error messages from qhull code */

  idx = NULL;
  dim = *ncol_in;
  n = *nrow_in;
  
  if (dim <= 0 || n <= 0)
    error("Invalid input matrix.");

  j=0;
  pt_array = (double *)malloc(n*dim*sizeof(double)); 

  for (i=0; i < n; i++)
    for (j=0; j < dim; j++)
      pt_array[dim*i+j] = x_in[dim*i+j];

  ismalloc = False;   /* True if qhull should free points in qh_freeqhull() or reallocation */

   /* hmm lots of options for qhull here */

  sprintf(flags,"qhull Qt7");

  exitcode = qh_new_qhull (dim,n,pt_array,ismalloc,flags,outfile,errfile);

  if (!exitcode)  /* 0 if no error from qhull */
    {  
    facetT *facet;                  /* set by FORALLfacets */
    vertexT *vertex, **vertexp;	/* set by FORALLfacets */
	
    *nf = qh num_facets;

    /* allocate the memory for the output */
    idx = (int *)malloc((*nf)*dim*sizeof(int));
	
    qh_vertexneighbors();
		
    i=0;
    FORALLfacets 
      {
      j=0;
		
      FOREACHvertex_ (facet->vertices) 
        {
	if (j >= dim)
	  warning("extra vertex %d of facet %d = %d", j++,i,1+qh_pointid(vertex->point));
	else
          idx[i+(*nf)*j++] = qh_pointid(vertex->point);
	}
      
      if (j < dim) warning("facet %d only has %d vertices",i,j);
	
		i++;
       }
		j=0;
	      
     }
	qh_freeqhull(!qh_ALL);	/*free long memory */
	qh_memfreeshort (&curlong, &totlong);	/* free short memory and memory allocator */

    if (curlong || totlong) {
	  warning("convhulln: did not free %d bytes of long memory (%d pieces)",totlong, curlong);
	}
	free(pt_array);
	return(idx);
}
