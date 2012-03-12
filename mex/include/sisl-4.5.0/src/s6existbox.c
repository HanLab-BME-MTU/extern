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
 * $Id: s6existbox.c,v 1.2 2001-03-19 15:59:01 afr Exp $
 *
 */


#define S6EXISTBOX

#include "sislP.h" 

#if defined(SISLNEEDPROTOTYPES)
int s6existbox(SISLbox *pbox,int itype,double aepsge)
#else
int s6existbox(pbox,itype,aepsge)
     SISLbox *pbox;
     int    itype;
     double aepsge;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Test if a particular box exist within an existing
*              box instance.
*
*
*
* INPUT      : pbox   - Box to test.
*              itype  - Kind of box to test existance of.
*                       = 0 : Do not expand box.
*                       = 1 : Make a totally expanded box.
*                       = 2 : Make a box expanded in the inner of the
*                             object, and reduced along the edges/endpoints.
*              aepsge - Geometry resolution.
*
* OUTPUT     : s6existbox -  Status.
*                            -1 : Kind of box exist, but is expanded
*                                 with another tolerance.
*                             0 : Kind of box do not exist.
*                             1 : Requested box exist.
*                                                                     
*
*
* METHOD     :
*
*
* REFERENCES :
*
*-
* CALLS      :  
*
* WRITTEN BY : Vibeke Skytt, SI, 91-01.
*
*********************************************************************
*/                                     
{
   if (pbox->e2min[itype] == SISL_NULL) return(0);  /* No box is made. */
   
   if (itype != 0 && DNEQUAL(pbox->etol[itype],aepsge))
      return(-1);  /* Box exist, but with another size of the expansion. */
   
   return(1);
}
   
