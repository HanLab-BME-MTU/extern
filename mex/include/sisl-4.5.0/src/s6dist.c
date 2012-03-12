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
 * $Id: s6dist.c,v 1.2 2001-03-19 15:59:01 afr Exp $
 *
 */


#define S6DIST

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
double 
s6dist(double epoint1[],double epoint2[],int idim)
#else
double s6dist(epoint1,epoint2,idim)
     double epoint1[];
     double epoint2[];
     int    idim;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Compute the distance between the points epoint1 and
*              epoint2.
*
*
*
* INPUT      : epoint1 - First point in distance calculation.
*              epoint2 - Second point in distance calculation.
*              idim    - Dimension of the space in which the points lie.
*
*
*
* OUTPUT     : s6dist  - Distance between the points.
*
*
* METHOD     : Compute lenght of the vector epoint1-epoint2.
*
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Vibeke Skytt, SI, 88-06.
*
*********************************************************************
*/                                     
{
  register double *s1,*s2,*s3; /* Pointers used to travers epoint1 and epoint2
				  arrays.                                      */
  register double tdist=DZERO; /* Distance between the points.                 */
  
  for (s1=epoint1,s2=epoint2,s3=epoint1+idim; s1<s3; s1++,s2++)
    tdist += (*s1 - *s2)*(*s1 - *s2);
  
  return(sqrt(tdist));
}
