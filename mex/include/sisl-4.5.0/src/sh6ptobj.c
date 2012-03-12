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
 * $Id: sh6ptobj.c,v 1.2 2001-03-19 16:06:03 afr Exp $
 *
 */


#define SH6PTOBJ

#include "sislP.h"


#if defined(SISLNEEDPROTOTYPES)
void 
sh6ptobj(double *point, SISLObject *obj, double aepsge,
	   double start[], double result[], int *jstat)
#else
void sh6ptobj(point, obj, aepsge, start, result, jstat)
     double *point;
     SISLObject *obj;
     double aepsge;
     double start[];
     double result[];
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Newton iteration on the distance function between
*              a point and an object to find a closest point
*              or an intersection point.
*
*
* INPUT      : point     - Point in the intersection.
*              obj       - Object in the intersection.
*              aepsge    - Geometry resolution.
*              start[]   - Start parameter value of the iteration on
*                          the object.
*
*
*
* OUTPUT     : result[]  - Parameter value of the object in intersection
*                          point.
*              jstat   - status messages  
*                                = 2   : A minimum distanse found.
*                                = 1   : Intersection found.
*                                < 0   : error.
*
*
* METHOD     : Newton iteration.
*
*
* REFERENCES :
*
*
* WRITTEN BY : Arne Laksaa, SI, Sep 1991
*              UJK, SI, Sep 1991, include point object
*********************************************************************
*/                       
{                        
  int kstat = 0;            /* Local status variable.                      */
  int kpos = 0;             /* Position of error.                          */
  double pstart[2];
  double pend[2];
  SISLPoint *sislpt = SISL_NULL;
  double loc_start[2];
  
  /* Test input.  */
  
  if (obj == SISL_NULL) goto err106;
		   
  if ( obj->iobj == SISLSURFACE)
  {
     if ((sislpt = newPoint(point, obj->s1->idim, 0)) == SISL_NULL)
        goto error;

     memcopy(loc_start,start,2,double);
     
     pstart[0] = obj->s1->et1[obj->s1->ik1 - 1];
     pstart[1] = obj->s1->et2[obj->s1->ik2 - 1];
     pend[0]   = obj->s1->et1[obj->s1->in1];
     pend[1]   = obj->s1->et2[obj->s1->in2];
     
     s1773(sislpt, obj->s1, aepsge,
	   pstart, pend, loc_start, result, &kstat);
     if (kstat < 0) goto error;
  }
  else if ( obj->iobj == SISLCURVE)
  {
     if ((sislpt = newPoint(point, obj->c1->idim, 0)) == SISL_NULL)
        goto error;
     
     pstart[0] = obj->c1->et[obj->c1->ik - 1];
     pend[0]   = obj->c1->et[obj->c1->in];
  
     loc_start[0] = start[0];
     s1771(sislpt, obj->c1, aepsge,
	   pstart[0], pend[0], loc_start[0], result, &kstat);
     if (kstat < 0) goto error;
  }
  else if ( obj->iobj == SISLPOINT)
  {
     if(s6dist(point,obj->p1->ecoef,obj->p1->idim) < aepsge)
	kstat = 1;
     else
        kstat = 2;
  }
  else goto err106;
  
  *jstat = kstat;
  goto out;
  
  /* Error in input. */
  
 err106: *jstat = -106;
  s6err("sh6ptobj",*jstat,kpos);
  goto out;                  
  
  /* Error in lower level routine.  */
  
  error : *jstat = kstat;
  s6err("sh6ptobj",*jstat,kpos);
  goto out;                  
	 
 out:    if (sislpt != SISL_NULL) freePoint(sislpt);
}
