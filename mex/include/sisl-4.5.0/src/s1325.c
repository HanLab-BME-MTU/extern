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
 * $Id: s1325.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S1325

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
double 
s1325(double aradiu,double angle)
#else
double s1325(aradiu,angle)
     double aradiu;
     double angle;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To create the tangent length for interpolating a
*              circular arc with an almost equi-oscillating Hermit qubic
*
* INPUT      : aradiu  - The radius of the circular arc
*              angle   - The opening angle of the circular arc
*
* OUTPUT     : s1325   - The proposed tangent length
*
* METHOD     : A second degree equation giving the tanget length is
*              solved
*
* WRITTEN BY : Tor Dokken, SI, Oslo, Norway. 30. June 1988
*                                  
*********************************************************************
*/
{
  double tcos,tsin;          /* Dummy variables                     */
  double ta,tb,tc,tl;        /* Dummy variables                     */
  double tconst = (double)1.85530139760811990992528773586425;
                             /* Constant used in the calculation    */
  
  
  
  tcos = cos(angle);
  tsin = sin(angle);
  
  /*  Calculate length of tangents
   *   tconst = (3-2sqrt(2))**1/3 + (3+2sqrt(2))**1/3 - 0.5 */
  
  ta     = (double)0.6*tconst - (double)0.9*tcos;
  tb     = ((double)0.4*tconst+(double)1.8)*tsin;
  tc     = ((double)0.4*tconst+(double)1.0)
           * tcos - (double)0.4*tconst - (double)1.0;
  tl     = aradiu*(-tb+sqrt(tb*tb-4*ta*tc))/((double)2.0*ta);
  
  return(tl);
}
