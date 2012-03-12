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
 * $Id: s6affdist.c,v 1.2 2001-03-19 15:59:00 afr Exp $
 *
 */


#define S6AFFDIST

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
double
     s6affdist(double e1[],double e2[],double emat[],int idim)
#else
double s6affdist(e1,e2,emat,idim)
   double e1[];
   double e2[];
   double emat[];
   int    idim;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : Compute the distance between two points using an
*              affine metric described by the matrix emat.
*
*
*
* INPUT      : e1     - First point.
*              e2     - Second point.
*              emat   - Matrix of affine metric.
*              idim   - Dimension of geometry space.
*              
*
* OUTPUT     : s6affdist - Distance between two points.
*
*
* METHOD     : 
*
* REFERENCES :
*
*-
* CALLS      : 
*
* WRITTEN BY : Vibeke Skytt, SI, 91-03.
*
*********************************************************************
*/
{
   int ki,kj;              /* Counters.  */
   double tdist = DZERO;   /* Distance.  */
   
   for (ki=0; ki<idim; ki++)
      for (kj=0; kj<idim; kj++)
	 tdist += emat[ki*idim+kj]*(e1[ki]-e2[ki])*(e1[kj]-e2[kj]);
   
   tdist = sqrt(idim*tdist);
   
   return tdist;
}
