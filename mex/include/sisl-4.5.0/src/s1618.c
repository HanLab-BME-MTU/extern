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
 * $Id: s1618.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S1618

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void s1618(double ematrix[], double eright[], double esol[], int in,
	   double *adiff)
#else
void s1618(ematrix, eright, esol, in, adiff)
     double ematrix[];
     double eright[];
     double esol[];
     int    in;
     double *adiff;
#endif
/*
*************************************************************************
*
* PURPOSE: To test the result of the solution of an equation system.
*
* INPUT:
*        Ematrix - The interpolation matrix. (length in*in.)
*        Eright  - The right hand side of the equation system.
*                  (length in.)
*        Esol    - The solutiuon of the equation system. (length in.)
*        In      - The dimension of the equation system.
*
* OUTPUT:
*        Adiff   - The maximal difference between left and right hand side.
*
* METHOD:
*        The product of the matrix and the solution is calculated
*        and subtracted from the right hand side.
*-
* Calls: No.
*
* Written by: A.M. Ytrehus, SI Oslo Oct.91.
* After FORTRAN (P1618), written by: T. Dokken  SI.
*****************************************************************
*/
{
  int ki, kj;
  int kdim;
  double tdiff, tdum;
  double tmax = (double) 0.0;
  double tmxel = (double) 0.0;


  /* Find greatest element of matrix. */

  kdim = in *in;

  for (ki = 0; ki < kdim; ki++)
    {
      tdum = fabs (ematrix[ki]);
      if (tdum > tmxel)
	tmxel = tdum;
    }
  if (tmxel == (double) 0.0)
    tmxel = (double) 1.0;


  /* Calculate product of matrix and solution, and
     subract from r.h.side. */

  for (ki = 0; ki < in; ki++)
    {
      tdum = (double) 0.0;

      for (kj = 0; kj < in; kj++)
	tdum += ematrix[ki * in +kj] *esol[kj];

      tdiff = tdum - eright[ki];

      tdum = fabs (tdiff) / tmxel;
      if (tdum > tmax)
	tmax = tdum;
    }

  *adiff = tmax;

  return;
}
