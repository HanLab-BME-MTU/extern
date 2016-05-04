//===========================================================================
// SISL - SINTEF Spline Library, version 4.5.0.
// Definition and interrogation of NURBS curves and surfaces. 
//
// Copyright (C) 2000-2005, 2010 SINTEF ICT, Applied Mathematics, Norway.
//
// This program is free software; you can redistribute it and/or          
// modify it under the terms of the GNU General Public License            
// as published by the Free Software Foundation version 2 of the License. 
//
// This program is distributed in the hope that it will be useful,        
// but WITHOUT ANY WARRANTY; without even the implied warranty of         
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          
// GNU General Public License for more details.                           
//
// You should have received a copy of the GNU General Public License      
// along with this program; if not, write to the Free Software            
// Foundation, Inc.,                                                      
// 59 Temple Place - Suite 330,                                           
// Boston, MA  02111-1307, USA.                                           
//
// Contact information: E-mail: tor.dokken@sintef.no                      
// SINTEF ICT, Department of Applied Mathematics,                         
// P.O. Box 124 Blindern,                                                 
// 0314 Oslo, Norway.                                                     
//
// Other licenses are also available for this software, notably licenses
// for:
// - Building commercial software.                                        
// - Building software whose source code you wish to keep private.        
//===========================================================================

#include "sisl-copyright.h"

/*
 *
 * $Id: s6lusolp.c,v 1.2 2001-03-19 15:59:02 afr Exp $
 *
 */


#define S6LUSOLP

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s6lusolp(double ea[],double eb[],int nl[],int im,int *jstat)
#else
void s6lusolp(ea,eb,nl,im,jstat)
     double ea[];
     double eb[];
     int    nl[];
     int    im;
     int    *jstat;
#endif
/*
************************************************************************
*
***********************************************************************
*
*   PURPOSE : Solve the equationsystem LU=eb by forth- and back-
*             substitution. L and U are stored in ea.
*
*
*   INPUT   : ea   - The factorized coeffecient-matrix.
*             im   - The number of equations.
*             nl   - Ordering of lines in the matrix.
*
*
*   INPUT/OUTPUT : eb - the right side of the equationsystem and
*                       the found values for the unknowns.
*
*   OUTPUT  : jstat  - Status variable.
*                        < 0 : Error
*                        = 0 : ok
*                        > 0 : Warning
*
*                                                                       
*   METHOD  : Solve on the equation-system LU=eb.
*
*
*   REFERENCES : Cheney & Kincaid : Numerical Mathematics and
*                                   Computing.
*
*-
*   CALLS      :
*
*   WRITTEN BY : Vibeke Skytt, SI, 86-10.
*
************************************************************************
*/
{
  int kpos = 0;      /* Position of error.                             */
  int ki,kj;         /* Counters.                                      */
  double *sx = SISL_NULL; /* Array used to keep solution of equation system
			internally.                                    */
  double tdiv;       /* Dividend in expression.                        */
  
  /* Allocate space for local array.  */
  
  if ((sx = newarray(im,double)) == SISL_NULL) goto err101;
  
  for (ki=0; ki<im-1; ki++)
    {
      /*  Gauss on right side of equation  */
      
      for (kj=ki+1; kj<im; kj++)      
	eb[nl[kj]] -= eb[nl[ki]]*ea[ki+nl[kj]*im];
    }
  
  tdiv = ea[im-1+nl[im-1]*im];
  if (DEQUAL(tdiv,DZERO)) goto warn1;
  sx[im-1] = eb[nl[im-1]]/tdiv;
  
  for (ki=im-2; ki>=0; ki--)
    {
      /*  Backwards substitution.   */
      
      for (kj=ki+1; kj<im; kj++)
	eb[nl[ki]] -= sx[kj]*ea[kj+nl[ki]*im];
      
      tdiv = ea[ki+nl[ki]*im];
      if (DEQUAL(tdiv,DZERO)) goto warn1;
      sx[ki] = eb[nl[ki]]/tdiv;
    }   
  for (ki=0; ki<im; ki++) eb[ki] = sx[ki];
  
  /* Equation system solved.  */
  
  *jstat = 0; 
  goto out;

/* Singular equation system.  */

warn1 : *jstat = 1;
        goto out;

/* Error in space allocation.  */

err101: *jstat = -101;
        s6err("s6lusolp",*jstat,kpos);
        goto out;

out:

/* Free space occupied by local array.  */

if (sx != SISL_NULL) freearray(sx);

return;
}
