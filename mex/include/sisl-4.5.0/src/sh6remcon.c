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
 * $Id: sh6remcon.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define SH6REMCON

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
      sh6remcon(SISLIntdat **pintdat,SISLIntpt *pt1,SISLIntpt *pt2,int *jstat)
#else
void sh6remcon(pintdat,pt1,pt2,jstat)
   SISLIntdat **pintdat;
   SISLIntpt *pt1;
   SISLIntpt *pt2;
   int       *jstat;
#endif   
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Remove connection between pt1 and pt2.
*
* INPUT      : pintdat  - Pointer to pointer to the SISLIntdat data.
*              pt1      - First point.
*              pt2      - Second point.
*              jstat    - Error flag.
*                         jstat =  0  => successful
*                         jstat = -1  => pt1 and pt2 are not connected
*                         jstat = -2  => error in space allocation
*                         jstat <  0  => error in lower level routine
*
*
* METHOD     : 
*
* CALLS      : s6err      - Gives error message.
*
* REFERENCES :
*
* WRITTEN BY : Michael Floater, SI, Oslo, Norway. June 91.
*
*********************************************************************
*/
{
  int kstat;                 /* Local status variable.                     */
  
   *jstat = 0;
  

  sh6disconnect(pt1,pt2,&kstat);
  if (kstat < 0) goto error;

  /* Do we want to remove either point from pintdat? */

  
  goto out;  
  

  /* Error in lower level routine. */

  error : *jstat = kstat;
  s6err("sh6remcon",*jstat,0);
  goto out;                       
  
  
  out: ;
}
