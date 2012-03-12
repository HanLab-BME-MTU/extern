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
 * $Id: sh6getlist.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define SH6GETLIST

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
      sh6getlist(SISLIntpt *pt1,SISLIntpt *pt2,int *index1,int *index2,int *jstat)
#else
void sh6getlist(pt1,pt2,index1,index2,jstat)
   SISLIntpt *pt1;
   SISLIntpt *pt2;
   int       *index1;
   int       *index2;
   int       *jstat;
#endif   
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Given two Intpts, find the two indices of the
*              (unique) link which connects them.
*              If no such link exists, index1 and index2 are set to -1.
*              pt1 and pt2 can be the same point.
*              Error if data structure
*              is inconsistent.
*
*
* INPUT      : pt1       - Pointer to first Intpt.
*              pt2       - Pointer to second Intpt.
*
*
* OUTPUT     : index1    - Index of list at pt1.
*              index2    - Index of list at pt2.
*              jstat     - Error flag.
*                         jstat =  0  => OK.
*                         jstat =  1  => No connection, indices = -1.
*                         jstat = -1  => Data structure inconsistent.
*
*
* METHOD     : 
*
*
* REFERENCES :
*
* WRITTEN BY : Michael Floater, SI, Oslo, Norway. June 91.
*
*********************************************************************
*/
{
   *index1 = -1;
   *index2 = -1;

   *jstat=0;

   /* Find "next" link from pt1 to pt2. */

   *index1 = sh6getprev(pt1,pt2);
   *index2 = sh6getprev(pt2,pt1);

   if(*index1 >= 0 && *index2 < 0) goto err1;
   if(*index2 >= 0 && *index1 < 0) goto err1;

   if(*index1 < 0 && *index2 < 0) *jstat=1;


   goto out;

err1:
   /* Error --  bad data structure. */

   *jstat = -1;
   s6err("sh6getlist",*jstat,0);
   goto out;
   
   out :
      return;
}
