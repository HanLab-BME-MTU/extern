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
 * $Id: s6dplane.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S6DPLANE

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
double
  s6dplane(double eq1[],double eq2[],double eq3[],double epoint[],
	   int idim,int *jstat)
#else
double s6dplane(eq1,eq2,eq3,epoint,idim,jstat)
   double eq1[];
   double eq2[];
   double eq3[];
   double epoint[];
   int    idim;
   int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Compute the distance between a plane given by three
*              points in the plane and a fourth point.
*
*
*
* INPUT      : eq1    - Point in the plane
*              eq2    - Point in the plane
*              eq3    - Point in the plane
*              epoint - Point. Dimension is idim.
*              idim   - Dimension of geometry space.
*
*
*
* OUTPUT     : s6dplane - Distance between plane and point.
*              jstat   - status messages  
*                        > 0      : warning
*                        = 0      : ok
*                        < 0      : error
*
*
* METHOD     : 
*              
*
*
* REFERENCES :
*
*-
* CALLS      : s6scpr   -  Scalar product between two vectors.
*              s6diff   -  Difference vector between two vectors.  
*              s6norm   -  Normalize vector.
*              s6crss   -  Cross product between two vectors.
*              s6dist   -  Distance between points.
*
* WRITTEN BY : Vibeke Skytt, SI, 91-02.
*
*********************************************************************
*/
{
   int kstat = 0;         /* Local status varaible.           */
   double tdist;          /* Distance between point and line. */
   double snorm[3];       /* Normal vector to the plane.      */
   double sdiff1[3];      /* Difference vector between points in the plane. */
   double sdiff2[3];      /* Difference vector between points in the plane. */
   double sdiff3[3];      /* Difference vector.               */
   
   /* Test dimension.     */
   
   if (idim != 3) goto err104;
   
   /* Compute difference vectors.  */
   
   s6diff(eq2,eq1,idim,sdiff1);
   s6diff(eq3,eq1,idim,sdiff2);
   s6diff(epoint,eq1,idim,sdiff3);
   
   /* Compute normalized plane normal.  */
   
   s6crss(sdiff1,sdiff2,snorm);
   (void)s6norm(snorm,idim,snorm,&kstat);
   
   /* Compute distance to closest point in plane. */
   
   if (kstat)
      tdist = fabs(s6scpr(sdiff3,snorm,idim));
   else 
      tdist = s6dist(eq1,epoint,idim);   /* Normal of zero length.  */

   /* Set status.  */
   
   *jstat = 0;
   goto out;
   
   /* Error in input, dimension not equal to 3.  */
   
   err104 : *jstat = -104;
   goto out;
   
   out :
      return tdist;
 }
