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
 * $Id: s1333count.c,v 1.2 2001-03-19 15:58:45 afr Exp $
 *
 */


#define S1333_COUNT
#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1333_count(int inbcrv,SISLCurve *vpcurv[],int *jcont,int *jstat)
#else
void s1333_count(inbcrv,vpcurv,jcont,jstat)
     int    	inbcrv;
     SISLCurve  *vpcurv[];
     int        *jcont;
     int        *jstat;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To count the continuity at the start and end of the curves
*              based on the multiplicity of knots and return the continuity.
*
*
* INPUT      : inbcrv - Number of curves in the curve-set.
*              vpcurv  - Array (length inbcrv) of pointers to the
*                       curves in the curve-set.
*
* OUTPUT     : jstat  - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*              *jcount - The continuity at the start or end based on
*                        multiplicity of knots
*-
* CALLS      : 
*
* WRITTEN BY : Tor Dokken  SI  Oslo,Norway.  Feb 1992
*
*********************************************************************
*/
{
  int kmult1,kmult2,kmult;   /* Multiplicities */
  int kcont=0;               /* Continuity so far */
  int kpos=0;
  int kleft = 0;
  int kstat;
  int ki;
  SISLCurve *curve=SISL_NULL;     /* Pointer to curve being tested */

  *jcont = -1;
  

  for (ki=0 ; ki<inbcrv ; ki++)
    {
       curve = vpcurv[ki];
       kmult1 = s6knotmult(curve->et,curve->ik,curve->in,&kleft,
                           curve->et[curve->ik-1],&kstat);
       if (kstat<0)goto error;
       
       kmult2 = s6knotmult(curve->et,curve->ik,curve->in,&kleft,
                           curve->et[curve->in],&kstat);
       if (kstat<0)goto error;
       kmult = MAX(kmult1,kmult2);
       kmult = MIN(kmult,curve->ik);
       
       if (ki==0)
	 {
	   
	   kcont = curve->ik - kmult - 1;
	 }
       else
	 {   
           kcont = MIN(kcont, curve->ik - kmult - 1); 
	 }
     }
  

  /* Task done */
  
  *jcont = kcont;
  
  *jstat = 0;
  goto out; 
  
  /* Error in lower level routine.  */

  error : 
    *jstat = kstat;     
  s6err("s1333_count",*jstat,kpos);
  goto out;
 out:
   
  return;
}
