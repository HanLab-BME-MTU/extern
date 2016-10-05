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
 * $Id: sh6tohelp.c,v 1.2 2001-03-19 16:06:04 afr Exp $
 *
 */


#define SH6TOHELP

#include "sislP.h"


#if defined(SISLNEEDPROTOTYPES)
void 
      sh6tohelp(SISLIntpt *pt,int *jstat)
#else
void sh6tohelp(pt,jstat)
   SISLIntpt *pt;
   int       *jstat;
#endif   
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Check if pt is a mai point. If it is, transform to
*              help point, if not give a message.
*
*
* INPUT      : pt       - Pointer to the Intpt.
*
*
* OUTPUT     : pt       - Pointer to the updated Intpt.
*              jstat    - Error flag.
*                         jstat =  0  => pt was a main point, now help.
*                         jstat =  1  => pt is not a main point.
*                         jstat = -1  => Error in pt.
*                         jstat = -2  => Illegal to convert status
*                         jstat = -3  => Error in data structure.
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
   int kstat; /* Local status */
   int num; 

   *jstat=0;

   if(pt == SISL_NULL) goto err1;

   if(sh6ismain(pt))  /* If pt is a help point. */
   {
      /* ??????????? */
      /* if(pt->no_of_curves > 2) goto err2; */
      
      num=sh6nmbmain(pt,&kstat);
      /* Problem in sh6edgred when starting reduction */
      /* if(num > 1) goto err2; */

       pt->iinter = -pt->iinter;  /* Convert status to main point. */
   }
   else
   {
       *jstat=1;
   }

   goto out;
   

err1:
   /* Error in input. pt is null. */
   
   *jstat = -1;
   s6err("sh6tohelp",*jstat,0);
   goto out;
   
   /* Error, Illegal to change status. */
   
   /* err2:
    *jstat = -2;
   s6err("sh6tohelp",*jstat,0);
   goto out; */
   
   
   out :
      return;
}
