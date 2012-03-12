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
 * $Id: s1235.c,v 1.3 2001-03-19 15:58:42 afr Exp $
 *
 */


#define S1235

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1235(double et[],int in,int ik,int *jnbreak,double **gbreak,int *jstat)
#else
void s1235(et,in,ik,jnbreak,gbreak,jstat)
     double et[];
     int    in;
     int    ik;
     int    *jnbreak;
     double **gbreak;
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Find break points in a knot vector. The first and last
*              parameter values are break values.
*
*
*
* INPUT      : et     - Knot vector to find break points in.
*              in     - Number of vertices of the curve corresponding
*                       to the knot vector.
*              ik     - Order of the curve corresponding to et.
*
*
*
* OUTPUT     : jnbreak - Number of break points found.
*              gbreak  - Array containing parameter values of break points.
*              jstat   - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : The knot vector has a break point at a knot if the 
*              multiplicity of the knot is ik-1.
*
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Vibeke Skytt, SI, 88-11.
* REVISED BY : Vibeke Skytt, SINTEF, 9801. Correction in loop counter
*                                          for the periodic case.
*
*********************************************************************
*/
{
  int kpos = 0;   /* Position of error.                    */
  int kj;         /* Counter.                              */
  int kbreak;     /* Current number of break points found. */
  int kmult;      /* Multiplisity of current knot.         */
  double tprev;   /* Value of previous knot.               */
  double *sbreak; /* Pointer into break point array.       */
  double *st;     /* Pointer into knot vector.             */
  
  /* Allocate space for an array that is at least as great as the
     number of break points.                                       */
  
  *gbreak = SISL_NULL;
  if ((*gbreak = newarray(in+2,double)) == SISL_NULL) goto err101;
  
  /* Set local pointer to and counter of break points.  */
  
  sbreak = *gbreak;
  kbreak = 0;
  
  /* Find break point in start of parameter interval and internal breaks. */
  
  tprev = et[ik-1];
  kmult = ik - 1;
  for (st=et+ik,kj=ik; kj<in; st++,kj++)
    {
      
      if (*st == tprev) kmult++;
      else
	{
	  if (kmult >= ik-1)
	    {
	      
	      /* New break point found.  */
	      
	      *(sbreak++) = tprev;
	      kbreak++;
	    }
	  tprev = *st;
	  kmult = 1;
	}
    }
  
  /* Find break point in end of interval.  */
  
  if (et[in] != tprev && kmult >= ik-1)
    {
      
      /* Remember last internal break point.  */
      
      *(sbreak++) = tprev;
      kbreak++;
    }
  *(sbreak++) = et[in];
  kbreak++;
  
  /* Reduce break point array to correct size.  */
  
  if (kbreak < in+2)
    if ((*gbreak = increasearray(*gbreak,kbreak,double)) == SISL_NULL) goto err101;
  
  /* Break points found.  */
  
  *jnbreak = kbreak;
  *jstat = 0;
  goto out;
  
  /* Error in space allocation.  */
  
 err101: 
  *jstat = -101;
  s6err("s1235",*jstat,kpos);
  goto out;
  
 out: 
  return;
}
