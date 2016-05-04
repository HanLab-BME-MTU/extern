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
 * $Id: sh6inspnt.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define SH6INSERTPT

#include "sislP.h"


#if defined(SISLNEEDPROTOTYPES)
void
sh6insertpt (SISLIntpt * pt1, SISLIntpt * pt2, SISLIntpt * ptnew, int *jstat)
#else
void 
sh6insertpt (pt1, pt2, ptnew, jstat)
     SISLIntpt *pt1;
     SISLIntpt *pt2;
     SISLIntpt *ptnew;
     int *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Place a new intersection point between two (connected)
*              points.
*              pt1 and pt2 must be linked.
*
* INPUT      : pt1      - First point.
*              pt2      - Second point.
*              ptnew    - Point to be placed between pt1 and pt2.
*              jstat    - Error flag.
*                         jstat =  0  => successful
*                         jstat = -1  => pt1 and pt2 are not connected
*                         jstat <  0  => error in lower level routine
*
*
* METHOD     :
*
*
* REFERENCES :
*
* WRITTEN BY : Michael Floater, SI, Oslo, Norway. June 91.
*              UJK, Connect before disconnect, to
*              keep the index order
*********************************************************************
*/
{
  int kstat;			/* Local status variable.                  */
  int index1=0,index2=0;
  int crv_dir1=0,crv_dir2=0;
  
  *jstat = 0;

  sh6getlist (pt1, pt2, &index1, &index2, &kstat);
  if (kstat < 0)
    goto error;			/* Error. */
  if (kstat == 1)
    goto err1;			/* pt1 and pt2 are not linked. */

  /* Save info in curve_dir */
  crv_dir1 = pt1->curve_dir[index1];
  crv_dir2 = pt2->curve_dir[index2];


  /* Check pt1,pt2,ptnew. */

  sh6connect (pt1, ptnew, &kstat);
  if (kstat < 0)
    goto error;			/* Error. */

  /* Set values in curve_dir */
  sh6getlist (pt1, ptnew, &index1, &index2, &kstat);
  pt1->curve_dir[index1]   = crv_dir1;
  ptnew->curve_dir[index2] = crv_dir2;

  sh6connect (pt2, ptnew, &kstat);
  if (kstat < 0)
    goto error;			/* Error. */

  /* Set values in curve_dir */
  sh6getlist (pt2, ptnew, &index1, &index2, &kstat);
  pt2->curve_dir[index1] = crv_dir2;
  ptnew->curve_dir[index2] = crv_dir1;


  sh6disconnect (pt1, pt2, &kstat);
  if (kstat < 0)
    goto error;			/* Error. */
  if (kstat == 1)
    goto err1;			/* pt1 and pt2 are not linked. */


  goto out;


/* Error. pt1 and pt2 are not linked.  */

err1:*jstat = -1;
  s6err ("sh6insertpt", *jstat, 0);
  goto out;

/* Error in sub function.  */

error:*jstat = kstat;
  s6err ("sh6insertpt", *jstat, 0);
  goto out;

out:
  return;
}
