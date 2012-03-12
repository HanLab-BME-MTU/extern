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
 * $Id: s6lprj.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */
#define S6LPRJ

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
double
s6lprj(double e1[],double e2[],int idim)
#else
double s6lprj(e1,e2,idim)
     double e1[];
     double e2[];
     int    idim;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To calculate the length of the projection of a vector on 
*              an arbitrary axis.
*
* INPUT      : e1      - The vector to be projected
*              e2      - The arbitrary axis vector
*              dim     - Number of dimensions in the space the vectors lie
*
* OUTPUT     : s6lprj  - The length of the projection vector
*
* METHOD     : The length of the projection vector is calculated as:
*                       __     __
*                       e1 dot e2
*              ||e3|| = ---------
*                       __     __  1/2
*                      (e2 dot e2)
*-
* CALLS      : s6scpr, s6length
*
* WRITTEN BY : Per Evensen, SI, Oslo, Norway. 1991-aug-16
*                                  
*********************************************************************
*/
{
  int kstat;
  double scpr1,scpr2,lproj; 

  scpr1 = s6scpr(e1,e2,idim);
  scpr2 = s6length(e2,idim,&kstat);
  if (kstat == 0) scpr2=0.000001;
  
  lproj = scpr1/scpr2;
  return(lproj);
}

