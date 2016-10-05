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
 * $Id: s6chpar.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S6CHPAR

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s6chpar(double ecoef1[],int in1,int in2,int idim,double ecoef2[])
#else
void s6chpar(ecoef1,in1,in2,idim,ecoef2)
     double ecoef1[];
     int    in1;
     int    in2;
     int    idim;
     double ecoef2[];
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : Change parameter directions of vertices of surface.
*
*
*
* INPUT      : ecoef1 - Vertices of original surface.
*              in1    - Number of vertices in first parameter direction.
*              in2    - Number of vertices in second parameter direction.
*              idim   - Dimension of the space in which the surfac lies.
*
*
* OUTPUT     : ecoef2 - Vertices after the changing of parameter directions.
*
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
* WRITTEN BY : Vibeke Skytt, SI, 88-11.
*
*********************************************************************
*/
{
  register int ki,kj,kk;  /* Counters.  */
  
  for (ki=0; ki<in1; ki++)
    for (kj=0; kj<in2; kj++)
      for (kk=0; kk<idim; kk++)
	ecoef2[(ki*in2+kj)*idim+kk] = ecoef1[(kj*in1+ki)*idim+kk];
}
