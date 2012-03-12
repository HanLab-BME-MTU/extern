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
 * $Id: s6curvrad.c,v 1.2 2001-03-19 15:59:01 afr Exp $
 *
 */


#define S6CURVRAD

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
      s6curvrad(double epnt1[],double epnt2[],double etang[],int idim,
		double *crad,int *jstat)
#else	 
void s6curvrad(epnt1,epnt2,etang,idim,crad,jstat)
     int idim,*jstat;
     double epnt1[],epnt2[],etang[],*crad;
#endif     
/*
*********************************************************************
*                                                                   
* PURPOSE    : Given both endpoints of a curve segment and the tangent in
*              one of the endpoints, estimate the curvature radius of the 
*              curve segment in the endpoint where the tangent is given.
*
*
* INPUT      : epnt1   - First endpoint.
*              epnt2   - Second endpoint.
*              etang   - Given tangent.
*              idim    - Dimension of geometry space.
*                       
*
* OUTPUT     : crad    - Estimated curvature radius.
*              jstat   - status messages  
*                                         = 1      : Tangent vector zero.
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : Pass a circle through the points having the given tangent
*              in the given endpoint. Compute the radius of this circle.
*
* REFERENCES : 
*
* USE        : 
*
*-
* CALLS      : s6dist - Distance between two points.  
*              s6diff - Difference vector between two vectors. 
*
* WRITTEN BY : Vibeke Skytt, SI, 05.90.
*
*********************************************************************
*/
{
  int kstat = 0;        /* Status variable.                  */
  double tdist;         /* Distance between the endpoints.   */
  double tdot;          /* Scalar product between tangent and vector
                           between endpoints.                */
  double tlmid;         /* Length of tangent vector.         */
  double tcos;          /* Cosinus to the angle between the tangent
                           and the vector between the endpoints. */
  double tang;          /* Angle at the centre of the circle
                           between the vectors from origo to
                           the given endpoints.              */
  double tdum;          /* Denominator in expression to find 
                           the  radius of the circle.        */
  double trad;          /* The curvature radius, i.e. the radius
                           of the circle.                    */
  double sdiff[3];      /* Difference vector between the endpoints. */
  
  /* Test input.  */

  if (idim != 3) goto err104;
  
  /* Estimate curvature radius based on endpoints and tangent.  */

  tdist = s6dist(epnt1,epnt2,idim);
  s6diff(epnt2,epnt1,idim,sdiff);

  tdot = s6scpr(etang,sdiff,idim);
  tlmid = s6length(etang,idim,&kstat);
  
  tcos = (tlmid*tdist != DZERO) ? fabs(tdot/(tlmid*tdist)) : fabs(tdot);
  tcos = MIN((double)1.0,tcos);
  
  tang = 2*acos(tcos);
  tdum = sqrt(2-2*cos(tang));
  trad = (tdum > REL_COMP_RES) ? tdist/tdum : -1.0;
  
  /* Set curvature radius. */
  
  *crad = trad;
  *jstat = 0;
  goto out;

  /* Error in input.  */

  err104 :
    *jstat = -104;
  goto out;

  out :
    return;
}    
