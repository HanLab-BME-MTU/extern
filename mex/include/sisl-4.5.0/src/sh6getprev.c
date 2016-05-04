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
 * $Id: sh6getprev.c,v 1.2 2001-03-19 15:59:07 afr Exp $
 *
 */


#define SH6GETPREV

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
int 
      sh6getprev(SISLIntpt *pt1,SISLIntpt *pt2)
#else
int sh6getprev(pt1,pt2)
   SISLIntpt *pt1;
   SISLIntpt *pt2;
#endif   
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Given an Intpt pt1 and a pointer to another Intpt pt2,
*              fetch the index of the pt1 array corresponding
*              to pt2. If no such index exists return -1.
*
*
* INPUT      : pt1       - Pointer to the Intpt.
*              pt2     - Pointer to another Intpt.
*
*
*
* METHOD     : 
*
*
* REFERENCES :
*
* WRITTEN BY : Michael Floater, SI, Oslo, Norway. May 91.
*
*********************************************************************
*/
{
   int       ncurv;   /* number of curves pt1 is connected to       */
   int       index;   /* index number for pnext array              */

   index = -1;

   if(pt1 == SISL_NULL || pt2 == SISL_NULL) goto out;

   ncurv = pt1->no_of_curves;  /* note ncurv can be zero */

   index=0;
   while(index < ncurv && pt1->pnext[index] != pt2) index++;
   if(index == ncurv) index = -1;  /* no index found */

   goto out;

   out :
      return index;
}
