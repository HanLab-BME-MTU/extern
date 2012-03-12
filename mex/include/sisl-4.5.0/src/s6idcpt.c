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
 * $Id: s6idcpt.c,v 1.2 2001-03-19 15:59:01 afr Exp $
 *
 */


#define S6IDCPT

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s6idcpt(SISLIntdat *pintdat,SISLIntpt *pintpt,SISLIntpt **rintpt)
#else
void s6idcpt(pintdat,pintpt,rintpt)
     SISLIntdat *pintdat;
     SISLIntpt  *pintpt;
     SISLIntpt  **rintpt;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : To find the point which is closest to pintpt
*              in the parametric space. If pintpt is the only
*              point in pintdat *rintpt is SISL_NULL.
*
*
*
* INPUT       :pintpt   - Pointer to an intersection point.
*              pintdat  - Pointer to intersection data.
*
*
* OUTPUT     : rintpt   - Pointer to a pointer to a point closest
*                         to pintpt.
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
* WRITTEN BY : Arne Laksaa, 05.89.
*
*********************************************************************
*/                                     
{
  if (pintdat == SISL_NULL)
    *rintpt = SISL_NULL;
  else
    {
      int ki,knr;                /* Counters.          */
      double tdist,td;           /* To store distanse. */
      
      if (pintpt == pintdat->vpoint[0])
        tdist = HUGE;
      else
        tdist = s6dist(pintdat->vpoint[0]->epar,pintpt->epar,pintpt->ipar);
      
      for (knr=0,ki=1; ki<pintdat->ipoint; ki++)
        {
	  if (pintpt == pintdat->vpoint[ki])
	    td = HUGE;
	  else
	    td = s6dist(pintdat->vpoint[ki]->epar,pintpt->epar,pintpt->ipar);
	  
	  if (td < tdist)
	    {
	      knr = ki;
	      tdist = td;
	    }
        }
      
      if (tdist == HUGE)
        *rintpt = SISL_NULL;
      else
        *rintpt = pintdat->vpoint[knr];
    }
}

