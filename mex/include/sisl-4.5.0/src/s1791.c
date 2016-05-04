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
 * $Id: s1791.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S1791

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
int 
s1791(double et[],int ik,int in)
#else
int s1791(et,ik,in)
     double et[];
     int    ik;
     int    in;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Test if it is possible to insert new internal nots
*              any further.
*
*
*
* INPUT      : et     - The knot vector.
*              ik     - The order of the curve/surface.
*              in     - The number of the basic functions.
*
*
*
* OUTPUT     : s1791  - Result of the test.
*                       = 0 : It is not possible to insert new knots
*                             any further.
*                       = 1 : It is possible to insert new knots.
*
*
* METHOD     :
*
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Arne Laksaa, SI, 89-04.
*
*********************************************************************
*/                                     
{
  register double tstart= et[ik - 1];
  register double tend  = et[in];
  register double tmid  = (tstart+tend)*(double)0.5;
  
  /* Check if it is possible to divide the parameter interval.  */
  
  if (DEQUAL(tmid,tstart) || DEQUAL(tmid,tend)) 
    return  0;
  else 
    return  1;
}





