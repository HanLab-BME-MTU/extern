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
 * $Id: s6mulvec.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S6MULVEC

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s6mulvec (double ematrix[], double evect[], double eright[])
#else
void
s6mulvec (ematrix, evect, eright)
     double ematrix[];
     double evect[];
     double eright[];
#endif
/*
*************************************************************************
*
* Purpose: To multiply a 4*4 matrix by a 3-D vector.
*
* Input:
*        Ematrix - The matrix (length 16).
*        Evect   - The vector (length 3).
*
* Output:
*        Eright   - The resulting 3-D vector.
*
* Written by: A.M. Ytrehus, SI Oslo Nov.91.
* After FORTRAN (P6mvec), written by: S. Meen  SI.
*****************************************************************
*/
{
  double svect[4];
  double sright[4];
  double tdum;
  int ki;
  int kj;

  /* Store the 3D evect-array in the 4D svect-array. */

  for (ki = 0; ki < 3; ki++)
    svect[ki] = evect[ki];

  svect[3] = (double) 1.0;


  /* Multiply the matrix by the vector svect.
     The result is stored in the 4-D local vector sright.*/

  for (ki = 0; ki < 4; ki++)
    {
      tdum = (double) 0.0;

      for (kj = 0; kj < 4; kj++)
	tdum += ematrix[4 * ki + kj] * svect[kj];

      sright[ki] = tdum;
    }


  /* 
   * Check if the bottom line of the matrix is 0,0,0,1.
   * In that case, just store the first three elements of svect in evect.
   * Else, divide all 3 first elements of svect by svect[3]. 
   * --------------------------------------------------------------------
   */

  for (ki = 0; ki < 3; ki++)
    eright[ki] = sright[ki];

  if (!(ematrix[12] == (double) 0.0) && (ematrix[13] == (double) 0.0) &&
      (ematrix[14] == (double) 0.0) && (ematrix[15] == (double) 1.0))
    {
      tdum = sright[3];

      for (ki = 0; ki < 3; ki++)
	eright[ki] = sright[ki] / tdum;
    }

  goto out;

out:
  return;
}
