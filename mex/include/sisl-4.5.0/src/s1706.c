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
 * $Id: s1706.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S1706

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1706(SISLCurve *pc)
#else
void s1706(pc)
     SISLCurve *pc;
#endif
/*
*******************************************************************
*
*********************************************************************
*
* PURPOSE    : Turn the direction of a curve.
*              The start of the new parameter is the same as the start
*              of the old parameter.
*              This rutine turn the direction of the orginal curve.
*              If you want a copy with a turned direction, just
*              make a copy and turn the direction of the copy.
*
*
*
* INPUT      : pc      -The curve.
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
* WRITTEN BY : Arne Laksaa, SI, 88-06.
* REVISED BY : Johannes Kaasa, SI, Sep 1991 (Introduced NURBS).
*
********************************************************************/
{
  int  kk=pc->ik;             /* Order of the input curve.             */
  int  kn=pc->in;             /* Number of vertices in the input curve.*/
  int  kdim=pc->idim;         /* Dimensjon of the space in whice curve
				 lies.                                 */
  register double *s1,*s2;
  register double *s3; 	       /* Pointers used in loop.               */
  register double t1,t2;       /* Help variables.                      */
  
  /* Now curve to turn. */
  
  if (!pc) goto out;
  
  /* Here we are turning the knot vector such that the first
     element have the same value as the old first element. */
  
  for (s1=pc->et,s2=s1+kk+kn-1,t1=(*s1)+(*s2); s1<=s2; s1++,s2--)
    {
      t2 = *s1;
      *s1 = t1 - *s2;
      *s2 = t1 - t2;
    }
  
  /* Here we just turn the vertices. */
  
  for (s1=pc->ecoef,s2=s1+kdim*(kn-1); s1<s2; s2-=2*kdim)
    for (s3=s1+kdim; s1<s3; s1++,s2++)
      {
	t1 = *s1;
	*s1 = *s2;
	*s2 = t1;
      }

  /* If necessary turn rational vertices. */

  if (pc->ikind == 2 || pc->ikind == 4)
    {
      kdim++;
      for (s1=pc->rcoef,s2=s1+kdim*(kn-1); s1<s2; s2-=2*kdim)
        for (s3=s1+kdim; s1<s3; s1++,s2++)
          {
	    t1 = *s1;
            *s1 = *s2;
	    *s2 = t1;
          }
    }
  
 out:
  return;
}

