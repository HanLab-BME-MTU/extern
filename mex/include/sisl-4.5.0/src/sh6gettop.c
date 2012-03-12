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
 * $Id: sh6gettop.c,v 1.2 2001-03-19 15:59:08 afr Exp $
 *
 */


#define SH6GETTOP

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
      sh6gettop(SISLIntpt *pt,int ilist,int *left1,int *right1,int *left2,int *right2,int *jstat)
#else
void sh6gettop(pt,ilist,left1,right1,left2,right2,jstat)
   SISLIntpt *pt;
   int       ilist;
   int       *left1;
   int       *right1;
   int       *left2;
   int       *right2;
   int       *jstat;
#endif   
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Given an Intpt and a list which it lies in, get the
*              pre-topology information. If list does not
*              exist give an error message.
*
*
* INPUT      : pt       - Pointer to the Intpt.
*              ilist    - Index specifying a list.
*
*
* OUTPUT     : pt       - Pointer to the updated Intpt.
*              left1    - pre-topology data.
*              right1   - pre-topology data.
*              left2    - pre-topology data.
*              right2   - pre-topology data.
*              jstat    - Error flag.
*                         jstat =  0  => OK.
*                         jstat = -1  => ilist is out of range.
*                         jstat = -2  => Error.
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
   *jstat=0;

   /* Check pt. */

   if(pt == SISL_NULL) goto err2;

   if(ilist >= 0 && ilist < pt->no_of_curves)
   {
       *left1=pt->left_obj_1[ilist];
       *right1=pt->right_obj_1[ilist];
       *left2=pt->left_obj_2[ilist];
       *right2=pt->right_obj_2[ilist];
   }
   else if(pt->no_of_curves == 0 && ilist == 0)
   {
       *left1=pt->left_obj_1[0];
       *right1=pt->right_obj_1[0];
       *left2=pt->left_obj_2[0];
       *right2=pt->right_obj_2[0];
   }
   /* UJK */
   else if( ilist == -1)
   {
       *left1=pt->left_obj_1[0];
       *right1=pt->right_obj_1[0];
       *left2=pt->left_obj_2[0];
       *right2=pt->right_obj_2[0];
   }
   else goto err1;


   /* Data is set. */

   goto out;
   

err1:
   /* Error. ilist is out of range. */
   
   *jstat = -1;
   s6err("sh6gettop",*jstat,0);
   goto out;

err2:
   /* Error in input. pt is SISL_NULL. */
   
   *jstat = -2;
   s6err("sh6gettop",*jstat,0);
   goto out;
   
   
   out :
      return;
}
