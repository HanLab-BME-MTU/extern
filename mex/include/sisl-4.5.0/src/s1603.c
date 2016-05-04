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
 * $Id: s1603.c,v 1.2 2005-02-28 09:04:48 afr Exp $
 *
 */


#define S1603

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1603(SISLSurf *psurf,double *cmin1,double *cmin2,double *cmax1,double *cmax2,int *jstat)
#else
void s1603(psurf,cmin1,cmin2,cmax1,cmax2,jstat)
     SISLSurf   *psurf;
     double *cmin1;
     double *cmin2;
     double *cmax1;
     double *cmax2;
     int    *jstat;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To pick the parametrization of a B-spline surface
*
* INPUT      : pc     - The B-spline surface.   
*
* OUTPUT     : cmin1  - Start parameter in first parameter directon. 
*              cmin2  - Start parameter in second parameter directon.
*              cmax1  - End   parameter in first parameter directon. 
*              cmax2    End   parameter in second parameter directon.
*              jstat  - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error  
*             
* METHOD     : 
*
*
* REFERENCES :
*
*-                                                 
* CALLS      : s6err
*              
*
* WRITTEN BY : Qyvind Hjelle SI, Oslo, Norway. Nov 1988
*
*********************************************************************
*/
{
  int kpos=0;              /* Position of error          */
  
  /* Check surf pointer */
  
  if (!psurf) goto err118;
  
  /* Pick parametrization */
  
  *cmin1 = psurf->et1[psurf->ik1-1];
  *cmax1 = psurf->et1[psurf->in1];
  *cmin2 = psurf->et2[psurf->ik2-1];
  *cmax2 = psurf->et2[psurf->in2];
  
  *jstat = 0;
  goto out;
  
  /* Error in input, no B-spline surface given */
  
 err118: 
  *jstat = -118;
  s6err("s1603",*jstat,kpos);
  goto out;
  
 out:
  
  return;
}          
