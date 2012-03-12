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
 * $Id: s6crss.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S6CRSS

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s6crss(double e1[],double e2[],double e3[])
#else
void s6crss(e1,e2,e3)
     double e1[];
     double e2[];
     double e3[];
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To make the cross product of two 3-D vectors
*
* INPUT      : e1      - First 3-D vector
*              e2      - Second 3-D vector
*
* OUTPUT     : e3      - The vector containing the cross product e1xe2
*
*
* METHOD     : The cross product is calculated by using its definition.
*
*-
* CALLS      :
*
* WRITTEN BY : Tor Dokken, SI, Oslo, Norway. 1988-may-03
*
*********************************************************************
*/
{
  e3[0] = e1[1]*e2[2] - e1[2]*e2[1];
  e3[1] = e1[2]*e2[0] - e1[0]*e2[2];
  e3[2] = e1[0]*e2[1] - e1[1]*e2[0];
}
