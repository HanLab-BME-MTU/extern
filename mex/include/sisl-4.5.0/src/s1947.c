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


#define S1947

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
   s1947(double ea[], int nfirst[], int nlast[], int ik, int im, 
	 double etau[], int in, int incont, double ew[], int inlr, 
	 int *jnred, double efac[], int *jstat)
#else
   void s1947(ea, nfirst, nlast, ik, im, etau, in, incont, ew, inlr, 
	      jnred, efac, jstat)
      double ea[];
      int    nfirst[];
      int    nlast[];
      int    ik;
      int    im;
      double etau[];
      int    in;
      int    incont;
      double ew[];
      int    inlr;
      int    *jnred;
      double efac[];
      int    *jstat;
#endif
/*
*********************************************************************
* 
* PURPOSE    : To adjust a least squares problem originating from 
*              approximating a spline curve from a subspace of the
*              original spline space, in order to satisfy given
*              continuity constraints at the start- and endpoint of
*              the spline. 
* 
* INPUT      : ea     - Real array of dimension (im*ik) containing 
*                       the band part of the coefficient matrix of the 
*                       problem. This matrix has
*                       dimension im*in but since at most
*                       ik entries are nonzero in each row, it can
*                       be stored in a im*ik array together
*                       with two integer arrays indicating the position
*                       of the first and last nonzero elements in each
*                       row. In addition it is known that the last nonzero
*                       element of row ki is stored in ea[ki*ik+ik-1].
*              nfirst - Integer array of dimension (im) containing 
*                       pointers to the first nonzero element of each row 
*                       of the original matrix.
*              nlast  - Integer array of dimension (im) containing 
*                       pointers to the last nonzero element of each row 
*                       of the original matrix.
*	       ik     - The order of the spline space in the underlying least
*                       squares problem.
*              im     - The dimension of the original spline space.
*              etau   - The knot vector of the reduced spline space..
*              in     - The number of unknowns in the original least squares
*                       problem, i.e. the dimension of the reduced spline space.
*              incont - Number of continuity constraints. Number of columns
*                       of ew.
*              inlr   - Number of rows of ew.
*              
*
* 
* OUTPUT     : ew     - Corner element of the coefficient matrix originating
*                       from periodicity. The size of ew1 is incont*inlr.
*              jnred  - Number of lines removed from the band part of the 
*                       equation system due to removal of the latest 
*                       coefficients because of continuity across the seem.
*              efac   - Factors used to express the incont last coefficients
*                       by the incont first.
*              jstat      - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*             
* 
* METHOD     : 
*
*
* REFERENCES : 
*              
*
* USE        : The knot vectors of the current spline space has to have
*              ik-multiple knots at the start- and endpoint.
*
*-
* CALLS      :   
*
* WRITTEN BY : Vibeke Skytt, SINTEF Oslo, 01.95.
*
*********************************************************************
*/
{
   int ki, kj, kr, kih;       /* Counters.       */
   int khindx;                /* Help in index computing. */
   double tdiv, th;           /* Help variables. */
   double *sr1 = SISL_NULL;   /* Help array used to store factors in equality
			        constraints at the start of the spline curve.  */
   double *sr2 = SISL_NULL;   /* To store factors at the end of the curve.     */
   double *stw1 = SISL_NULL;  /* Knot intervals at the start of the curve.     */
   double *stw2 = SISL_NULL;  /* Knot intervals at the end of the curve.       */
   double *shelp = SISL_NULL; /* Help array to avoid overwriting.              */
      
   /* Test input. */
   
   if (DEQUAL(etau[ik-1], etau[ik]) || DEQUAL(etau[in-1], etau[in])) 
      goto err112;
   
   /* Allocate scratch for local arrays. */
   
   if ((sr1 = new0array(2*incont*incont+3*incont, DOUBLE)) == SISL_NULL) 
      goto err101;
   sr2 = sr1 + incont*incont;
   stw1 = sr2 + incont*incont;
   stw2 = stw1 + incont;
   shelp = stw2 + incont;
   
   /* Set up equality constraints for coefficients at the ends of the curve. */
     
   /* Store first and last knot interval. */
   
   stw1[0] = etau[ik] - etau[ik-1];
   stw2[0] = etau[in] - etau[in-1];
      
   /* Compute factors in the constraint equations. */
   
   sr1[0] = sr2[0] = (double)1;
   for (ki=1; ki<incont; ki++)
   {
      stw1[ki] = etau[ik+ki] - etau[ik-1];
      stw2[ki] = etau[in] - etau[in-ki-1];
      
      sr1[ki*incont] = -sr1[(ki-1)*incont]/stw1[0];
      sr2[ki*incont] = sr2[(ki-1)*incont]/stw2[0];
      for (kj=1; kj<=ki; kj++)
      {
	 sr1[ki*incont+kj] = sr1[(ki-1)*incont+kj-1]/stw1[kj-1] -
	    sr1[(ki-1)*incont+kj]/stw1[kj];
	 sr2[ki*incont+kj] = -sr2[(ki-1)*incont+kj-1]/stw2[kj-1] +
	    sr2[(ki-1)*incont+kj]/stw2[kj];
      }
   }

   /* Express the incont coefficients at the end of the spline by a combination
      of the coefficients at the start, starting with the last coefficient. 
      The result is stored in sr2. */
   
   for (ki=1; ki<incont; ki++)
   {
      tdiv = sr2[ki*incont+ki];
      for (kj=ki; kj>=0; kj--)
      {
	 th = sr1[ki*incont+kj];
	 for (kr = ki-1; kr>= kj; kr--)
	    th -= sr2[ki*incont+kr]*sr2[kr*incont+kj];
	 shelp[kj] = th/tdiv;
      }
      memcopy(sr2+ki*incont, shelp, ki+1, DOUBLE);
   }
   
   /* Use the fact that the incont last coefficients may be expressed as
      a known combination of the first incont coefficients to reduce the
      number of unknown in the equation system. */
   
   for (kj=im-1; kj>=0; kj--)
   {
      if (nlast[kj] < in-incont) break;  /* All continuity requirement treated. */
      
      for (kih=0, ki=0; ki<incont; ki++)
      {
	 if (nlast[kj] < in-ki-1) continue;  /* No contribution. */
	 
	 for (kr=0; kr<=ki; kr++)
	    ew[(kj-im+inlr)*incont+kr] += sr2[ki*incont+kr]*
	       ea[(kj+1)*ik-kih-1];
	    
	 kih++;
      }
	 
      
      /* Justify size of equation system. */
      
      if (nlast[kj] >= in-incont) 
      {
	 khindx = nlast[kj] - in + incont + 1;
	 
	 /* Move line of equation system to the left in order to have 
	    the last non-zero element at position kj*ik+ik-1 in ea. */
	 
	 for (ki=ik-1; ki>=ik-in+incont+nfirst[kj]; ki--)
	    ea[kj*ik+ki] = ea[kj*ik+ki-khindx];
	 for (; ki>=0; ki--)
	    ea[kj*ik+ki] = DZERO;
      }
      nlast[kj] = in-incont-1;
      if (nlast[kj] < nfirst[kj]) (*jnred)++;
      
      /* Test if the band part of the matrix and the corner element 
	 interfere, and in that case move elements from the band part
	 of the matrix to the corner element.                         */
      
      if (nfirst[kj] < incont)
      {
	 for (ki=nfirst[kj]; ki<incont; ki++)
	    ew[(kj-im+inlr)*incont+ki] += ea[kj*ik+ik-1-nlast[ki]+nfirst[ki]];
	 nfirst[kj] = incont;
      }
   }
   
   /* Equation system adusted. Copy factors to output array. */
   
   memcopy(efac, sr2, incont*incont, DOUBLE);
   *jstat = 0;
   goto out;
   
   /* Error in scratch allocation. */
   
   err101 : *jstat = -101;
   s6err("s1947", *jstat, 0);
   goto out;

   /* Error in knot vector. */
   
   err112 : *jstat = -112;
   s6err("s1947", *jstat, 0);
   goto out;
   
   out:
      /* Free scratch occupied by local arrays. */
      
      if (sr1 != SISL_NULL) freearray(sr1);
      
      return;
}
