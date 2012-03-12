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
 * $Id: s1774.c,v 1.2 2001-03-19 15:58:53 afr Exp $
 *
 */
#define S1774

#include "sislP.h"


#if defined(SISLNEEDPROTOTYPES)
void
s1774(SISLCurve *crv, double point[], int dim, double epsge,
	   double start, double end, double guess, double *clpar, int *stat)
#else
void s1774(crv, point, dim, epsge, start, end, guess, clpar, stat)
     SISLCurve  *crv;
     double point[];
     int dim;
     double epsge;
     double start;
     double end;
     double guess;
     double *clpar;
     int    *stat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Newton iteration on the distance function between
*              a curve and a point, to find a closest point or an
*              intersection point.
*              If a bad choice for the guess parameter is given in, the
*              iteration may end at a local, not global closest point.
*
*
* INPUT      : crv     - The curve in the closest point problem.
*              point   - The point in the closest point problem.
*              dim     - Dimension of the geometry.
*              epsge   - Geometrical resolution.
*              start   - Curve parameter giving the start of the search
*                        interval.
*              end     - Curve parameter giving the end of the search
*                        interval.
*              guess   - Curve guess parameter for the closest point
*                        iteration.
*
*
*
* OUTPUT     : clpar   - Resulting curve parameter from the iteration.
*              stat    - status messages
*                                = 2   : A minimum distanse found.
*                                = 1   : Intersection found.
*                                < 0   : error.
*
*
* METHOD     : Newton iteration.
*
*
* REFERENCES :
*
*
* WRITTEN BY : Johannes Kaasa, SINTEF, Aug 1995
*
*********************************************************************
*/
{
   int kpos = 0;             /* Error indicator. */
   SISLPoint* ppoint = SISL_NULL; /* SISL point.      */
   
   /* Generate a SISL point. */
   
   ppoint = newPoint(point, dim, 0);
   
   /* Call s1771. */
   
   s1771(ppoint, crv, epsge, start, end, guess, clpar, stat);
   if (*stat < 0) goto error;

   goto out;

   /* Error in lower level routine.  */

   error :
   s6err("s1774", *stat, kpos);
   goto out;

 out:    
    if (ppoint != SISL_NULL) freePoint(ppoint);
}

