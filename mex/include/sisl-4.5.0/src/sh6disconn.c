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
 * $Id: sh6disconn.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define SH6DISCONNECT

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
      sh6disconnect(SISLIntpt *pt1,SISLIntpt *pt2,int *jstat)
#else
void sh6disconnect(pt1,pt2,jstat)
   SISLIntpt *pt1;
   SISLIntpt *pt2;
   int       *jstat;
#endif   
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Remove the connection between pt1 and pt2 (which
*              must be unique) if any.
*
* INPUT      : pt1      - First point.
*              pt2      - Second point.
*              jstat    - Error flag.
*                         jstat =  0  => Successful.
*                         jstat =  1  => pt1 and pt2 are not connected
*                         jstat = -1  => Error in data structure.
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
  int kstat;                 /* Local status variable.            */
  int index1,index2;         /* Indices for pt1 and pt2.          */
  
   *jstat = 0;
  

   /* Check if pt1 and pt2 are connected. */

   sh6getlist(pt1,pt2,&index1,&index2,&kstat);
   if(kstat < 0) goto err1;
   if(kstat == 1)
   {
       *jstat = 1;
       goto out;
   }


   /* Disconnect. */

   pt1->no_of_curves--;
   pt1->pnext[index1] = pt1->pnext[pt1->no_of_curves];
   pt1->curve_dir[index1] = pt1->curve_dir[pt1->no_of_curves];

   pt2->no_of_curves--;
   pt2->pnext[index2] = pt2->pnext[pt2->no_of_curves];
   pt2->curve_dir[index2] = pt2->curve_dir[pt2->no_of_curves];

  
  goto out;  
  

  
  /* No connection exists. */

  err1 : *jstat = -1;
  s6err("sh6disconnect",*jstat,0);
  goto out;                       
  
  out: ;
}
