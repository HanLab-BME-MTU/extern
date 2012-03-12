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
 * $Id: s6scpr.c,v 1.2 2001-03-19 15:59:02 afr Exp $
 *
 */


#define S6SCPR

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
double 
s6scpr(double e1[],double e2[],int idim)
#else
double s6scpr(e1,e2,idim)
     double e1[];
     double e2[];
     int    idim;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To make the scalar product of two vectors
*
* INPUT      : e1      - The first vector in the scalar product
*              e2      - The second vector in the scalar product
*              idim    - Number of dimensions in the space the vectors lie
*
* OUTPUT     : s6scpr  - The value of the scalar product
*
* METHOD     : The scalar product is calculated according to the definition
*-
* CALLS      :
*
* WRITTEN BY : Tor Dokken, SI, Oslo, Norway. 1988-may-03
*                                  
*********************************************************************
*/
{
  register int ki;
  register double tsum=DZERO; 

  for (ki=0;ki<idim;ki++)
    tsum += e1[ki]*e2[ki];

  return(tsum);
}
