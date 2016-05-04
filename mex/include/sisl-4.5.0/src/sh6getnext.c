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
 * $Id: sh6getnext.c,v 1.2 2001-03-19 15:59:07 afr Exp $
 *
 */


#define SH6GETNEXT

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
SISLIntpt* 
      sh6getnext(SISLIntpt *pt,int index)
#else
SISLIntpt* sh6getnext(pt,index)
   SISLIntpt *pt;
   int       index;
#endif   
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Given an Intpt and an index, fetch the next point
*              given by index.
*              If error, return SISL_NULL.
*
*
* INPUT      : pt       - Pointer to the Intpt.
*              index    - Index of link at pt.
*
*
* OUTPUT     : 
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

   SISLIntpt *nextpt = SISL_NULL;

   /* check if index is within range */

   if(pt != SISL_NULL &&
      index >= 0 &&
      index < pt->no_of_curves) nextpt = pt->pnext[index];

   goto out;

   
   out :
      return nextpt;
}
