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
 * $Id: s1937.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S1937

#include "sislP.h"


#if defined(SISLNEEDPROTOTYPES)
void
s1937 (double et[], int iordr, int ref, int left, double alfa[], double etref[])
#else
void
s1937 (et, iordr, ref, left, alfa, etref)
     double et[];
     int iordr;
     int ref;
     int left;
     double alfa[];
     double etref[];

#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE: To calculate the discrete B-spline values number iref
*	   of the refinement etref of the array et.
*
*
* INPUT:   et	 - The original knot vector.
*	   iordr - The order of the discrete B-spline to be
*		   calculated.
*	   ref	 - The index of the discrete B-spline to be
*		   calculated.
*	   left	 - Arrayindex, satisfying:
*		   et[left-1] <= etref[ref-1] < et[left]
*	   etref - Refined knot vector.
*
*
* OUTPUT:  alfa	 - The values of the calculated discrete B-splines.
*		   alfa[0]    - Corresponds to number left-iordr+1.
*		   alfa[1]    - Corresponds to number left-iordr+2.
*		   alfa[left] - Corresponds to number left.
*
* METHOD: We use the Oslo-algorithm developped by Cohen, Lyche and
*         Riesenfeld.
*
* REFERENCES: Cohen, Lyche, Riesenfeld: Discrete B-splines and subdivision
*	      techniques in computer aided geometric design, computer
*	      graphics and image processing, vol 14, no.2 (1980)
*
* CALLS: No.
*
* WRITTEN BY :  Christophe R. Birkeland, SI, 1991-07
*
*********************************************************************
*/
{
  int ki, kl, kr, low;		/* Loop control parameters. 	*/
  int stop, start;		/* and array indicators.	*/
  double tj, td1, td2;		/* Parameters used to improve.	*/
  double beta1, beta;		/* algorithm.			*/


  /* We have et[left-1] <= etref[ref-1] < et[left]
     So the discrete B-splines can be calculated. */

  low = left - iordr;
  start = left - 1;
  stop = iordr - 1;
  alfa[stop] = 1;

  for (kr = 0; kr < stop; kr++)
    {
      beta1 = (double) 0.0;
      tj = etref[ref + kr];
      if (start < 0)
	start = 0;

      for (ki = start; ki < left; ki++)
	{
	  kl = ki - low;
	  td1 = tj - et[ki];
	  td2 = et[ki + kr + 1] - tj;
	  beta = alfa[kl] / (td1 + td2);
	  alfa[kl - 1] = td2 * beta + beta1;
	  beta1 = td1 * beta;
	}
      alfa[iordr - 1] = beta1;
      start--;
    }
  return;
}
