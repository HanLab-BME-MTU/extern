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
 * $Id: s6ang.c,v 1.2 2001-03-19 15:59:00 afr Exp $
 *
 */


#define S6ANG

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
double 
s6ang(double evec1[],double evec2[],int idim)
#else
double s6ang(evec1,evec2,idim)
     double evec1[];
     double evec2[];
     int    idim;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Compute the angle (in radians) between two vectors
*
*
*
* INPUT      : evec1   - First vector 
*              evec2   - Second vector 
*              idim    - Dimension of the space in which the vectors lie.
*
*
*
* OUTPUT     : s6ang   - Angle in radians between vectors
*
*
* METHOD     : Make cosine of the angle by computing the scalar product,
*              then divide by the length of the two vectors.
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Tor Dokken SI, 88-07.
*              Arne Laksaa SI, 89-07.
*
*********************************************************************
*/                                     
{
  double tscpr,tang,tlength1,tlength2,tcos;
  int    kstat1,kstat2;
  
  tscpr = s6scpr(evec1,evec2,idim);
  
  tlength1 = s6length(evec1,idim,&kstat1);
  tlength2 = s6length(evec2,idim,&kstat2);
  
  if (!kstat1 || !kstat2)
    tang = DZERO;
  else
    {
      tcos = fabs(tscpr/(tlength1*tlength2));
      tcos = MIN((double)1.0,tcos);
      tang = acos(tcos);
    }
  
  return(tang);
}
