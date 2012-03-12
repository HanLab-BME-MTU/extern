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
 * $Id: s6identify.c,v 1.2 2001-03-19 15:59:01 afr Exp $
 *
 */


#define S6IDENTIFY

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
    s6identify(SISLSurf* s,double a[], double b[], double level_val,
	       double eps1,double eps2,int* jstat)
#else
void s6identify(s,a,b,level_val,eps1,eps2,jstat)
     SISLSurf* s;
     double a[],b[];
     double level_val;
     double eps1,eps2;
     int* jstat;
#endif
/***********************************************************************
*
*********************************************************************
*
* PURPOSE     : To check if the to points a,b in the parameter plane of
*               s should be identified according to certain criteria.
*               See METHOD below for the criteria.
*
*
*
*
* INPUT      : s          - Pointer to 1-dimensional surface object.
*              a[0:1]     - first point in parameter plane
*              b[0:1]     - second point in parameter plane
*	       idim       - space dimension.
*              level_val  - constant crucial for the identification criterium.
*                           ( represents usually the constant surface with
*                             wich s is intersected )
*              eps1	  - radius of ball used in the criteria for
*			    separate points in parameter plane.
*              eps2       - resolution in space.
*
*
* OUTPUT     : 		- jstat      - status messages
*					  = 1      : identify points a,b
*                                         = 0      : a,b are separate points
*                                         < 0      : error
*
*
* METHOD     : a and b are identified iff (i) and (ii) are satisfied:
*                (i)           |a-b| <= eps1.
*		 (ii)          |c(g)-level_val| <= eps2 . 
*		 	Here
*			    g is the degree 3 polynomial Hermite interpolant 
*                           to the restriction of
*                           the surface s to the line segment [a,b], i.e.
*                           g(0) = s(a), g'(0) = Ds(a)(b-a),
*                           g(1) = s(b), g'(1) = Ds(b)(b-a) where
*                           Ds is the gradient of s.
*                           Furthermore,
*                           c(g) is the spline control polygon of g for the 
*                           spline representation of g with k-tuple knots 
*                           at 0,1/2,1.
*                           
*
*
* REFERENCES :
*
*-
* CALLS      :
*              
*
* WRITTEN BY : Kyrre Strom, SI, 93-01.
* MODIFIED BY :
*
**********************************************************************/
{
  double c[4],cref[8];
  int i,kstat;

  if ( s == SISL_NULL ||
      (a[0] < s->et1[0] || a[0] > s->et1[s->in1]) ||
      (a[1] < s->et2[0] || a[1] > s->et2[s->in2]) ||
      (b[0] < s->et1[0] || b[0] > s->et1[s->in1]) ||
      (b[1] < s->et2[0] || b[1] > s->et2[s->in2])   )
    goto err109;

  if (DEQUAL(a[0],b[0]) && DEQUAL(a[1],b[1]))
    {
      kstat = 1;
      goto out;
    }
  if ( sqrt((a[0]-b[0])*(a[0]-b[0]) + (a[1]-b[1])*(a[1]-b[1])) > eps1 )
    kstat = 0;
  else
    {
      s6hermite_bezier(s,a,b,1,c,&kstat);
      if (kstat < 0) goto error;

      s6deCasteljau(c,0.0,1.0,0.5,4,cref,&kstat);
      if (kstat < 0) goto error;

      kstat = 1;
      for (i=0; i<8; i++)
	if (fabs(cref[i]-level_val) > eps2)
	  kstat = 0;
    }
  
  goto out;

 err109: kstat = -109;
  s6err("s6identify",kstat,0);
  goto out;


 error: 
  s6err("s6identify",kstat,0);
  goto out;

 out: 
    *jstat = kstat;
    return ;
} 
