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
 * $Id: s6idklist.c,v 1.2 2001-03-19 15:59:01 afr Exp $
 *
 */


#define S6IDKLIST

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s6idklist(SISLIntdat **pintdat,SISLIntlist *pintlist,int *jstat)
#else
void s6idklist(pintdat,pintlist,jstat)
     SISLIntdat  **pintdat;
     SISLIntlist *pintlist;
     int     *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : To remove an intersection list including all intersection points
*              in the list. The mother pintdat is updated.
*              If pintdat is empty, pintdat is killed and set to SISL_NULL.
*
*
*
* INPUT/OUTPUT:pintlist - Pointer to a list.
*              pintdat  - Pointer to a pointer to intersection data.
*
*
* OUTPUT     :jstat    - status messages  
*                               = 1      : Pintlist is not in pintdat.
*                               = 0      : OK!
*                               < 0      : error
*
*
* METHOD     : 
*
*
* REFERENCES :
*
*-
* CALLS      : s6err      - Gives error message.
*              s6idkpt    - Kills a point.
*              freeIntpt  - free instant of intpt structure.
*
* WRITTEN BY : Ulf J. Krystad, SI, 08.89.
*
*********************************************************************
*/                                     
{
  SISLIntpt *qkillpt,*qnext,*qdum1,*qdum2;
  
  int ki,knum,kstat;  
  
  *jstat = 0;
  
  /* We have to be sure that we have an intdat structure. */
  
  if ((*pintdat) == SISL_NULL)
    goto out;
  
  if (pintlist == SISL_NULL)
    {
      *jstat = 1;
      goto out;
    }
  
  /* Now we have to find the index in the vlist array in pintdat. */
  
  
  for (ki=0,knum = -1; ki < (*pintdat)->ilist; ki++)
    if ((*pintdat)->vlist[ki] == pintlist)
      {
	knum = ki;
	break;
      }
  
  if (knum == -1)
    /* Not in the pintdat list. */
    *jstat = 1;
  else
    {
      pintlist->plast->pcurve = SISL_NULL;
      
      /* Kill all points in the list. */
      for (ki=0,qkillpt=pintlist->pfirst,qnext=qkillpt->pcurve;
	   qnext!=SISL_NULL;
	   qkillpt=qnext,qnext=qnext->pcurve)
	{
	  s6idkpt(pintdat,&qkillpt,&qdum1,&qdum2,&kstat);
	  if (kstat < 0) goto error;
	}
      s6idkpt(pintdat,&qkillpt,&qdum1,&qdum2,&kstat);
      if (kstat < 0) goto error;
      
      /* Update pintdat. */
      if ((*pintdat) != SISL_NULL)
	{
	  (*pintdat)->vlist[knum] = (*pintdat)->vlist[(*pintdat)->ilist-1];
	  ((*pintdat)->ilist)--;
	  (*pintdat)->vlist[(*pintdat)->ilist] = SISL_NULL;
	}
      freeIntlist(pintlist);
    }
  
  goto out;  
  
  error : *jstat = kstat;
  s6err("s6idklist",*jstat,0);
  goto out;                       
  
  out: ;
}
