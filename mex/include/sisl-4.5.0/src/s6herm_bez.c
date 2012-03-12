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
 * $Id: s6herm_bez.c,v 1.3 2001-03-19 15:59:01 afr Exp $
 *
 */


#define S6HERMITE_BEZIER

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
    s6hermite_bezier(SISLSurf* s,double a[],double b[],int idim, double c[],
		     int* jstat)
#else
void s6hermite_bezier(s,a,b,idim,c,jstat)
     SISLSurf *s;
     double a[],b[];
     int idim;
     double c[];
     int *jstat;
#endif
/*********************************************************************
*
*********************************************************************
*
* PURPOSE     : Returning the Hermite interpolant to the restriction of the
*               surface s to line segment [a,b].
*               The Hermite interpolant is a Bezier
*               curve of degree 3 parametrized over the intervall [0,1].
*
*
*
*
* INPUT      : s          - Pointer to surface object.
*              a[0:1]     - Start point of line segment.
*              b[0:1]     - End point of line segment.
*	       idim       - space dimension.
*              c[0:3*idim]- allocated space;
*
*
* OUTPUT     : c          - Bezier coeffs of the Hermite interpolant.
*              jstat      - status messages
*                                         < 0      : error
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
*
* WRITTEN BY : Kyrre Strom, SI, 93-01.
* Revised by : Paal Fugelli, SINTEF, Oslo, Norway, Nov.1994. Initialized
*              'jstat' to zero when no error.
*
**********************************************************************/
{
  int i,kstat,left1=0,left2=0;
  double dblocal[9];
  double *derive=SISL_NULL;


  if (DEQUAL(a[0],b[0]) && DEQUAL(a[1],b[1])) goto error;
  if (s->idim != idim) goto error;

  if ( idim > 3)
  {
    derive = newarray(3*idim,double);
    if (derive == SISL_NULL) goto err101;
  }
  else
    derive = dblocal;

  /* evaluate s and its derivative at a */

  s1424(s,1,1,a,&left1,&left2,derive,&kstat);
  if (kstat < 0) goto error;
  for (i=0; i < idim; i++)
  {
    c[i] = derive[i];
    c[idim+i] = c[i] + (derive[idim+i]*(b[0]-a[0])
			+ derive[2*idim+i]*(b[1]-a[1]))/3.0;
  }

  /* evaluate s and its derivative at b */

  s1424(s,1,1,b,&left1,&left2,derive,&kstat);
  if (kstat < 0) goto error;
  for (i=0; i < idim; i++)
  {
    c[3*idim+i] = derive[i];
    c[2*idim+i] = c[3*idim+i] - (derive[idim+i]*(b[0]-a[0])
				 + derive[2*idim+i]*(b[1]-a[1]))/3.0;
  }

  *jstat = 0;
  goto out;

  /* Error in space allocation.  */

  err101 :
    *jstat = -101;
    goto out;


  /* Error in lower level routine.  */

  error :
    *jstat = kstat;
    goto out;


  out :

    if (derive != SISL_NULL && derive != dblocal)
      freearray(derive);

  return;

}
