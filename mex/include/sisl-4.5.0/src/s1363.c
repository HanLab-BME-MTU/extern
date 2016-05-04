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
 * $Id: s1363.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S1363

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1363(SISLCurve *pc,double *cmin,double *cmax,int *jstat)
#else
void s1363(pc,cmin,cmax,jstat)
     SISLCurve  *pc;
     double *cmin;
     double *cmax;
     int    *jstat;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To pick the parametrization of a B-spline curve
*
* INPUT      : pc     - The B-spline curve.   
*
* OUTPUT     : 
*              jstat  - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*              cmin   - Start of parametrization of curve
*              cmax   - End of parametrization of curve
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
* WRITTEN BY : Tor Dokken, SI, Oslo, Norway. Nov 1988
*
*********************************************************************
*/
{
  int kstat;          /* Local status variable                           */
  int kpos=0;         /* Position of error                               */
  
  /* Check if curve is correct */
  
  s1707(pc,&kstat);
  if (kstat<0) goto error;
  
  /* Pick parametrization */
  
  *cmin = pc->et[pc->ik - 1];
  *cmax = pc->et[pc->in];
  
  *jstat = 0;
  goto out;
  
  /* Error in lower level function */
  
 error:  *jstat = kstat;
  s6err("s1363",*jstat,kpos);
  goto out;
 out:
  
  return;
}          
