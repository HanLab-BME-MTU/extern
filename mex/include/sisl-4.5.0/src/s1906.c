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
 * $Id: s1906.c,v 1.2 2001-03-19 15:58:55 afr Exp $
 *
 */


#define S1906

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1906 (double *epoint, int *etype, int icnsta, int icnend, int inbpnt,
       int idim, double **opoint, int **otype, int *knbpnt, int *jstat)
#else
void
s1906 (epoint, etype, icnsta, icnend, inbpnt, idim, opoint,
       otype, knbpnt, jstat)
     double *epoint;
     int    *etype;
     int    icnsta;
     int    icnend;
     int    inbpnt;
     int    idim;
     double **opoint;
     int    **otype;
     int    *knbpnt;
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE: To transform from one type of interpolation conditions to another.
*
* INPUT: epoint - Array (length idim*inbpnt) containing the 'old'
*                 interpolation conditions.
*        etype  - Array (length inbpnt) containing type indicator for
*                 points/derivatives/second-derivatives:
*                  1 - Ordinary point.
*                  2 - Knuckle point. (Is treated as an ordinary point.)
*                  3 - Derivative to next point.
*                  4 - Derivative to prior point.
*                  5 - Second derivative to next point.
*                  6 - Second derivative to prior point.
*                 13 - Start-point of tangent to next point.
*                 14 - End-point of tangent to prior  point.
*        icnsta - Additional condition at the start of the curve:
*                  0 : No additional condition.
*                  1 : Zero curvature at start.
*        icnend - Additional condition at the end of the curve:
*                  0 : No additional condition.
*                  1 : Zero curvature at end.
*        inbpnt - No. of points/derivatives in the epoint array.
*        idim   - The dimension of the space in which the points lie.
*
* OUTPUT:
*	 opoint - The new interpolation conditions. (length idim*inbpnt)
*        otype  - Array containing kind of condition. Dimension
*                 is inbpnt.
*                 =  0 : A point is given.
*                 =  d : The d'th derivatative condition to the
*                        previous point is given.
*                 = -d : The d'th derivatative condition to the
*                        next point is given.
*        knbpnt - no. of points.
*	 jstat	- Status variable.
*                  > 0     : warning
*                  = 0     : ok
*                  < 0     : error
*
* METHOD     :
*
* REFERENCES :
*
* CALLS      : s6err.
*
* WRITTEN BY :  Trond Vidar Stensby, SI, 1991-07
*
*********************************************************************
*/
{
  int kpos = 0;
  int count, count2;		/* Loop control variables. */
  int dummy;			/* Dummy variables used when adressing arrays. */
  int start, start2;

  *jstat = 0;


  /* Allocate output arrays. */

  if (icnsta && icnend)
    *knbpnt = inbpnt + 2;
  else if (icnsta || icnend)
    *knbpnt = inbpnt + 1;
  else
    *knbpnt = inbpnt;

  *opoint = newarray ((*knbpnt) * idim, DOUBLE);
  if (*opoint == SISL_NULL)
    goto err101;
  *otype = newarray (*knbpnt, INT);
  if (*otype == SISL_NULL)
    goto err101;


  /* Insert additional interpolation conditions. */

  if (icnsta != 0)
    {
      for (count = 0; count < idim; count++)
	(*opoint)[count] = (double) 0.0;

      (*otype)[0] = -2;
    }

  if (icnend != 0)
    {
      dummy = (*knbpnt) * idim;
      for (count = ((*knbpnt) - 1) * idim; count < dummy; count++)
	(*opoint)[count] = (double) 0.0;

      (*otype)[(*knbpnt) - 1] = 2;
    }

  /* Copy the rest of the points. */

  if (icnsta != 0)
    start = 1;
  else
    start = 0;

  for (count = 0; count < inbpnt; count++)
    {
      /* Transfor interpolation conditions. */

      if (etype[count] == 13)
	{
	  start2 = (count + start) * idim;
	  for (count2 = 0; count2 < idim; count2++)
	    {
	      (*opoint)[start2 + count2] = epoint[(count + 1) * idim + count2] -
		epoint[count * idim + count2];
	    }
	}
      else if (etype[count] == 14)
	{
	  start2 = (count + start) * idim;
	  for (count2 = 0; count2 < idim; count2++)
	    {
	      (*opoint)[start2 + count2] = epoint[count * idim + count2] -
		epoint[(count - 1) * idim + count2];
	    }
	}
      else
	{
	  start2 = (count + start) * idim;
	  for (count2 = 0; count2 < idim; count2++)
	    {
	      (*opoint)[start2 + count2] = epoint[count * idim + count2];
	    }
	}

      /* Transform derivative indicators. */

      if (etype[count] == 1 || etype[count] == 2)
	(*otype)[count + start] = 0;
      else if (etype[count] == 3)
	(*otype)[count + start] = -1;
      else if (etype[count] == 4)
	(*otype)[count + start] = 1;
      else if (etype[count] == 5)
	(*otype)[count + start] = -2;
      else if (etype[count] == 6)
	(*otype)[count + start] = 2;
      else if (etype[count] == 13)
	(*otype)[count + start] = -1;
      else if (etype[count] == 14)
	(*otype)[count + start] = 1;
    }

  /* Ok. */

  goto out;


  /* Error in scratch allocation. */

err101:
  *jstat = -101;
  s6err ("s1906", *jstat, kpos);
  goto out;

out:
  return;
}
