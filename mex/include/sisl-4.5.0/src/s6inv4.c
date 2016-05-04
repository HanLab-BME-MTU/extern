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
 * $Id: s6inv4.c,v 1.2 2005-02-28 09:04:50 afr Exp $
 *
 */



#define S6INV4

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s6inv4 (double em[], double einv[], int *jstat)
#else
void
s6inv4 (em, einv, jstat)
     double em[];
     double einv[];
     int *jstat;
#endif
/*
*************************************************************************
*
* Purpose: To invert a 4*4-matrix.
*
* Input:
*        Em      - Array (length 16) containing the matrix to be inverted.
*
* Output:
*        Einv    - The inverted matrix.
*
* Calls: s6err.
*
* Written by: A.M. Ytrehus, SI Oslo Feb.92.
*
*****************************************************************
*/
{
  int ki;
  double det;

  *jstat = 0;


  /* Calculate the determinant of the matrix em. */

  det = em[0] * em[5] * (em[10] * em[15] - em[14] * em[11])
    - em[0] * em[6] * (em[9] * em[15] - em[13] * em[11])
    + em[0] * em[7] * (em[9] * em[14] - em[13] * em[10])
    - em[1] * em[4] * (em[10] * em[15] - em[14] * em[11])
    + em[1] * em[6] * (em[8] * em[15] - em[12] * em[11])
    - em[1] * em[7] * (em[8] * em[14] - em[12] * em[10])
    + em[2] * em[4] * (em[9] * em[15] - em[13] * em[11])
    - em[2] * em[5] * (em[8] * em[15] - em[12] * em[11])
    + em[2] * em[7] * (em[8] * em[13] - em[12] * em[9])
    - em[3] * em[4] * (em[9] * em[14] - em[13] * em[10])
    + em[3] * em[5] * (em[8] * em[14] - em[12] * em[10])
    - em[3] * em[6] * (em[8] * em[13] - em[12] * em[9]);


  /* Calculate the inverse matrix. */

  einv[0] = em[5] * (em[10] * em[15] - em[14] * em[11])
    - em[6] * (em[9] * em[15] - em[13] * em[11])
    + em[7] * (em[9] * em[14] - em[13] * em[10]);

  einv[4] = -em[4] * (em[10] * em[15] - em[14] * em[11])
    + em[6] * (em[8] * em[15] - em[12] * em[11])
    - em[7] * (em[8] * em[14] - em[12] * em[10]);

  einv[8] = em[4] * (em[9] * em[15] - em[13] * em[11])
    - em[5] * (em[8] * em[15] - em[12] * em[11])
    + em[7] * (em[8] * em[13] - em[12] * em[9]);

  einv[12] = -em[4] * (em[9] * em[14] - em[13] * em[10])
    + em[5] * (em[8] * em[14] - em[12] * em[10])
    - em[6] * (em[8] * em[13] - em[12] * em[9]);


  einv[1] = -em[1] * (em[10] * em[15] - em[14] * em[11])
    + em[2] * (em[9] * em[15] - em[13] * em[11])
    - em[3] * (em[9] * em[14] - em[13] * em[10]);

  einv[5] = em[0] * (em[10] * em[15] - em[14] * em[11])
    - em[2] * (em[8] * em[15] - em[12] * em[11])
    + em[3] * (em[8] * em[14] - em[12] * em[10]);

  einv[9] = -em[0] * (em[9] * em[15] - em[13] * em[11])
    + em[1] * (em[8] * em[15] - em[12] * em[11])
    - em[3] * (em[8] * em[13] - em[12] * em[9]);

  einv[13] = em[0] * (em[9] * em[14] - em[13] * em[10])
    - em[1] * (em[8] * em[14] - em[12] * em[10])
    + em[2] * (em[8] * em[13] - em[12] * em[9]);


  einv[2] = em[1] * (em[6] * em[15] - em[14] * em[7])
    - em[2] * (em[5] * em[15] - em[13] * em[7])
    + em[3] * (em[5] * em[14] - em[13] * em[6]);

  einv[6] = -em[0] * (em[6] * em[15] - em[14] * em[7])
    + em[2] * (em[4] * em[15] - em[12] * em[7])
    - em[3] * (em[4] * em[14] - em[12] * em[6]);

  einv[10] = em[0] * (em[5] * em[15] - em[13] * em[7])
    - em[1] * (em[4] * em[15] - em[12] * em[7])
    + em[3] * (em[4] * em[13] - em[12] * em[5]);

  einv[14] = -em[0] * (em[5] * em[14] - em[13] * em[6])
    + em[1] * (em[4] * em[14] - em[12] * em[6])
    - em[2] * (em[4] * em[13] - em[12] * em[5]);


  einv[3] = -em[1] * (em[6] * em[11] - em[10] * em[7])
    + em[2] * (em[5] * em[11] - em[9] * em[7])
    - em[3] * (em[5] * em[10] - em[9] * em[6]);

  einv[7] = em[0] * (em[6] * em[11] - em[10] * em[7])
    - em[2] * (em[4] * em[11] - em[8] * em[7])
    + em[3] * (em[4] * em[10] - em[8] * em[6]);

  einv[11] = -em[0] * (em[5] * em[11] - em[9] * em[7])
    + em[1] * (em[4] * em[11] - em[8] * em[7])
    - em[3] * (em[4] * em[9] - em[8] * em[5]);

  einv[15] = em[0] * (em[5] * em[10] - em[9] * em[6])
    - em[1] * (em[4] * em[10] - em[8] * em[6])
    + em[2] * (em[4] * em[9] - em[8] * em[5]);


  for (ki = 0; ki < 16; ++ki)
    einv[ki] /= det;

  return;
}
