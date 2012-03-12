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
 * $Id: sh6tomain.c,v 1.2 2001-03-19 16:06:04 afr Exp $
 *
 */


#define SH6TOMAIN

#include "sislP.h"


#if defined(SISLNEEDPROTOTYPES)
void 
      sh6tomain(SISLIntpt *pt,int *jstat)
#else
void sh6tomain(pt,jstat)
   SISLIntpt *pt;
   int       *jstat;
#endif   
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Check if pt is a help point. If it is, transform to
*              main point, if not give a message.
*
*
* INPUT      : pt       - Pointer to the Intpt.
*
*
* OUTPUT     : pt       - Pointer to the updated Intpt.
*              jstat    - Error flag.
*                         jstat =  0  => pt was a help point, now main.
*                         jstat =  1  => pt is not a help point.
*                         jstat = -1  => Error in pt.
*
*
* METHOD     : 
*
*
* REFERENCES :
*
* WRITTEN BY : Michael Floater, SI, Oslo, Norway. June 91.
* MODYFIED BY: UJK, SI, Oslo, Norway. September 91.
*********************************************************************
*/
{
   int ki; /* Loop variable. */
   int num; 
   int kstat;
   
   *jstat=0;

   if(pt == SISL_NULL) goto err1;

   if(sh6ishelp(pt))  /* If pt is a help point. */
   {
       pt->iinter = -pt->iinter;  /* Convert status to main point. */

       /* Go through all neighbours and keep invariant:
	  not more than one mainpoint connected to a help point. */
       for(ki=0; ki<pt->no_of_curves; ki++) 
       {
	   if(sh6ishelp(pt->pnext[ki]))
	   {
	      /* UJK, change all NON-terminators to main */
	      /* num=sh6nmbmain(pt->pnext[ki],&kstat); */
	       num = pt->pnext[ki]->no_of_curves;
	       if(num > 1) sh6tomain(pt->pnext[ki],&kstat);
	   }
       }

   }
   else
   {
       *jstat=1;
   }

   goto out;
   

err1:
   /* Error in input. pt is null. */
   
   *jstat = -1;
   s6err("sh6tomain",*jstat,0);
   goto out;
   
   
   out :
      return;
}
