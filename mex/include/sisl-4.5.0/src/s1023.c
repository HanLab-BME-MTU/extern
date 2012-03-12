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
 * $Id: s1023.c,v 1.3 2001-03-19 15:58:41 afr Exp $
 *
 */


#define S1023

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1023(double center[], double axis[], double equator[], int latitude,
	  int longitude, SISLSurf **sphere, int *stat)
#else
void s1023(center, axis, equator, latitude, longitude, sphere, stat)
     double center[];
     double axis[];
     double equator[];
     int latitude;
     int longitude;
     SISLSurf  **sphere;
     int    *stat;
#endif
/*
*********************************************************************
*
* PURPOSE    : To describe octants of a sphere as a NURBS. This can also
*              be used to describe the complete sphere.
*
*
* INPUT      : center        - Center point of the sphere
*              axis          - Axis of the sphere (towards the north pole)
*              equator       - Vector from center to start point on the equator
*              latitude      - Flag indicating number of octants in north/south
*                              direction:
*                              = 1 : Octants in the northern hemisphere
*                              = 2 : Octants in both hemispheres
*              longitude     - Flag indicating number of octants along the
*                              equator:
*                              = 1 : Octants in 1. quadrant
*                              = 2 : Octants in 1. and 2. quadrant
*                              = 3 : Octants in 1., 2. and 3. quadrant
*                              = 4 : Octants in all quadrants
*                              This is counted counter-clockwise from equator
*
*
* OUTPUT     :
*              stat          - status messages
*                                             > 0      : warning
*                                             = 0      : ok
*                                             < 0      : error
*              spher         - Pointer to the sphere produced
*
* METHOD     :
*
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Johannes Kaasa, SI, Oslo, Norway, Jan. 93
* Revised by : Paal Fugelli and Johannes Kaasa, SINTEF, Oslo, Norway, 08-94.
*
*********************************************************************
*/
{
   int kstat;            /* Status variable.                            */
   int kpos=0;             /* Position of error.                          */
   int ki, kj, kl;       /* Indexes in for loops.                       */
   int in1;              /* Number of vertices along the longitude.     */
   int in2;              /* Number of vertices along the latitude.      */
   int ik1 = 3;          /* Order along the longitude.                  */
   int ik2 = 3;          /* Order along the latitude.                   */
   double *et1 = SISL_NULL;   /* Knot vector along the longitude.            */
   double *et2 = SISL_NULL;   /* Knot vector along the latitude.             */
   double *rcoef = SISL_NULL; /* Coefficients of the sphere.                 */
   int kind = 2;         /* Rational Bspline surface.                   */
   double weight;        /* Rational weight.                            */
   double radius;        /* Radius of the sphere.                       */
   double norm;          /* Length of vectors.                          */
   double x_axis[3];     /* axis normalized with radius length.         */
   double z_axis[3];     /* Radius vector normal to x_axis and equator. */
   double w1, w2;        /* Rational weights in both directions.        */
   double x_comp;        /* Component in local x direction.             */
   double y_comp;        /* Component in local y direction.             */
   double z_comp;        /* Component in local z direction.             */

   /* Do necessary initiation and allocation. */

   *sphere = SISL_NULL;
   weight = (double)1.0/sqrt(2.0);
   in1 = 1 + 2*latitude;
   in2 = 1 + 2*longitude;

   radius = s6length(equator, 3, &kstat);
   if (kstat < 0) goto error;
   norm = s6length(axis, 3, &kstat);
   if (kstat < 0) goto error;
   for (ki = 0; ki < 3; ki++)
      x_axis[ki] = radius*axis[ki]/norm;
   s6crss(x_axis, equator, z_axis);
   norm = s6length(z_axis, 3, &kstat);
   if (kstat < 0) goto error;
   for (ki = 0; ki < 3; ki++)
      z_axis[ki] = radius*z_axis[ki]/norm;

   if((et1 = newarray(in1 + ik1, DOUBLE)) == SISL_NULL) goto err101;
   if((et2 = newarray(in2 + ik2, DOUBLE)) == SISL_NULL) goto err101;
   if((rcoef = newarray(4*in1*in2, DOUBLE)) == SISL_NULL) goto err101;

   /* Initiate the knot vectors. */

   for (ki = 0; ki < ik1; ki++)
      et1[ki] = (double)0.;
   for (ki = 0; ki < latitude; ki++)
   {
      et1[ik1 + 2*ki]     = (ki + 1)*PIHALF;
      et1[ik1 + 2*ki + 1] = (ki + 1)*PIHALF;
   }
   et1[in1 + ik1 - 1] = latitude*PIHALF;

   for (ki = 0; ki < ik2; ki++)
      et2[ki] = (double)0.;
   for (ki = 0; ki < longitude; ki++)
   {
      et2[ik2 + 2*ki]     = (ki + 1)*PIHALF;
      et2[ik2 + 2*ki + 1] = (ki + 1)*PIHALF;
   }
   et2[in2 + ik2 - 1] = longitude*PIHALF;

   /* Initiate the coefficient vector. */

   for (ki = 0; ki < in2; ki++)
   {
      if (ki == 1 || ki == 3 || ki == 5 || ki == 7)
	 w2 = weight;
      else
         w2 = (double)1.;

      if (ki == 0 || ki == 1 || ki == 7 || ki == 8)
	 y_comp = (double)1.;
      else if (ki == 3 || ki == 4 || ki == 5)
	 y_comp = - (double)1.;
      else
         y_comp = (double)0.;

      if (ki == 1 || ki == 2 || ki == 3)
	 z_comp = (double)1.;
      else if (ki == 5 || ki == 6 || ki == 7)
	 z_comp = - (double)1.;
      else
         z_comp = (double)0.;

      for (kj = 0; kj < in1; kj++)
      {
	 if (kj == 1 || kj == 3)
	    w1 = weight;
	 else
	    w1 = (double)1.;

	 if (kj == 0 || kj == 1)
	    x_comp = (double)1.;
	 else if (kj == 3 || kj == 4)
	    x_comp = - (double)1.;
	 else
	    x_comp = (double)0.;

	 w1 *= w2;

         if (kj == 0 || kj == 4)
	 {
	    for (kl = 0; kl < 3; kl++)
	       rcoef[4*(ki*in1 + kj) + kl] = w1*(center[kl] + x_comp*x_axis[kl]);
         }
	 else
	 {
	    for (kl = 0; kl < 3; kl++)
	       rcoef[4*(ki*in1 + kj) + kl] =  w1*(center[kl] + x_comp*x_axis[kl]
	                        + y_comp*equator[kl] + z_comp*z_axis[kl]);
	 }
	 rcoef[4*(ki*in1 + kj) + 3] = w1;
      }
   }
   (*sphere) = newSurf(in1, in2, ik1, ik2, et1, et2, rcoef, kind, 3, 1);
   if ((*sphere) == SISL_NULL) goto err101;

   if (et1 != SISL_NULL) freearray(et1);
   if (et2 != SISL_NULL) freearray(et2);
   if (rcoef != SISL_NULL) freearray(rcoef);

   *stat = 0;
   goto out;

   /* Error in curve allocation.  */

   err101:
      *stat = -101;
      s6err("s1023",*stat,kpos);
      goto out;

   /* Error in lower level routine. */

   error:
      *stat = kstat;
      s6err("s1023", *stat, kpos);
      goto out;

   out:
      return;
}
