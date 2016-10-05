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
 * $Id: sh6trmlist.c,v 1.2 2001-03-19 16:06:04 afr Exp $
 *
 */


#define SH6TRIMLIST

#include "sislP.h"


#if defined(SISLNEEDPROTOTYPES)
void
sh6trimlist (SISLIntpt * pt, SISLIntpt *** ptlist, int *no_of_points,
	     int *no_alloc)
#else
void
sh6trimlist (pt, ptlist, no_of_points, no_alloc)
     SISLIntpt *pt;
     SISLIntpt ***ptlist;
     int *no_of_points;
     int *no_alloc;

#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : To find the maximal trim boundary (of coincidence)
*              containg the given point pt.
*
*
* INPUT      : pt            - Pointer to int point to be examined
*
*
* INPUT/OUTP:  ptlist        - Pointer to an array containing
*                              pointers to all intersection points
*                              that are contigous trim neighbours.
*               no_of_points - Number of points in ptlist array
*               no_alloc     - Allocation size of ptlist array
* OUTPUT     :
*              jstat     - Error flag.
*                         jstat =  0  => OK.
*                         jstat = -1  => Data structure inconsistent.
*
*
* METHOD     :
*
*
* REFERENCES :
*
* WRITTEN BY : Ulf J. Krystad, SI, Oslo, Norway. October 91.
*
*********************************************************************
*/
{
  int clean_up = FALSE;		/* Clean up on top level */
  int incr = 20;		/* Allocation size       */
  int ki;			/* Loop control          */
  /* --------------------------------------------------- */


  /* Check if point is a TRIM point */
  if (pt->iinter != SI_TRIM)
    goto out;

  /* Check if point is treated */
  if (pt->marker == -90)
    goto out;

  /* Mark point as treated */
  pt->marker = -90;


  if (*no_alloc <= *no_of_points)
    {
      if (*no_alloc == 0)
	{
	  clean_up = TRUE;
	  (*no_alloc) += incr;
	  *ptlist = newarray (*no_alloc, SISLIntpt *);
	  if (*ptlist == SISL_NULL)
	    goto out;
	}
      else
	{
	  clean_up = FALSE;
	  (*no_alloc) += incr;
	  *ptlist = increasearray (*ptlist, *no_alloc, SISLIntpt *);
	  if (*ptlist == SISL_NULL)
	    goto out;
	}
    }

  /* Fill in */
  (*ptlist)[*no_of_points] = pt;
  (*no_of_points)++;

  /* Treat all neighbours */
  for (ki = 0; ki < pt->no_of_curves; ki++)
    sh6trimlist (pt->pnext[ki], ptlist, no_of_points, no_alloc);


/* Must unmark the points in array if no_alloc == 0 */
  if (clean_up)
    for (ki = 0; ki < (*no_of_points); ki++)
      (*ptlist)[ki]->marker = 0;

  goto out;


out:
  return;
}
