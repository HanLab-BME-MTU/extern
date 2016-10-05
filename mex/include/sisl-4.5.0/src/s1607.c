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
 * $Id: s1607.c,v 1.3 2001-03-19 15:58:51 afr Exp $
 *
 */


#define S1607

#include "sislP.h"


#if defined(SISLNEEDPROTOTYPES)
void
s1607 (SISLCurve * pc1, SISLCurve * pc2, double aepsge,
       double aend1, double afil1, double aend2, double afil2,
       int itype, int idim, int ik, SISLCurve ** rc, int *jstat)
#else
void
s1607 (pc1, pc2, aepsge, aend1, afil1, aend2, afil2, itype,
       idim, ik, rc, jstat)
     SISLCurve *pc1;
     SISLCurve *pc2;
     double aepsge;
     double aend1;
     double afil1;
     double aend2;
     double afil2;
     int itype;
     int idim;
     int ik;
     SISLCurve **rc;
     int *jstat;
#endif
/*
*********************************************************************
*
* PURPOSE    : To calculate a fillet curve between two curves.
*              Parameter values indicate between which points
*              the fillet is to be produced.
*
* INPUT      : pc1    - The first input curve.
*              pc2    - The second input curve.
*              aepsge - Geometry resolution.
*              aend1  - Parameter value on the first curve telling that
*                       the part of the curve lying on this side of afil1
*                       shall not be replaced by the fillet.
*              afil1  - Parameter value of the starting point of the fillet
*                       on the first curve.
*              aend2  - Parameter value on the second curve telling that the
*                       part of the curve lying on this side of afil2 shall
*                       not be replaced by the fillet.
*              afil2  - Parameter value of the starting point of the fillet on
*                       the second curve.
*              itype  - Indicator of type of fillet.
*                     = 1  - Circle, interpolating tangent on first curve,
*                            not on curve 2.
*                     = 2  - Conic if possible
*                     else - Polynomial segment
*              idim   - Dimension of space.
*              ik     - Order of fillet curve.
*
* OUTPUT     : rc     - Fillet curve produced
*              jstat  - status messages
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
* METHOD     : The points and tangents specified by the input parameter values
*              are calculated, then the interpolation curve of type given by
*              the  fillet indicator is calculated if possible. If a conic is
*              specified and not possible to make, a qubic spline segment is
*              made in stead.
*
* USE        :
*
* REFERENCES :
*
* NOTE	     : Routine s1611 has been changed to ensure
*	       that output status value may be positive (=9005).
*              The test line 201-219 in this routine might
*	       be changed. The limit is to be evaluated more
*	       precisely.
*
* CALLS      : s1707,s1221,s1227,s6norm,s6dist,s6crss,s1611,s1334,s6err.
*
*
* WRITTEN BY : Christophe Birkeland, SI, 1991-07
* REVISED BY : Christophe Birkeland, SI, July 92 (parameter TYPE changed
*              back to DOUBLE)
* REVISED BY : Johannes Kaasa, SI, Aug. 1992 (Changed status 9005 to 105
*              and -9006 to -104).
* REVISED BY : Johannes Kaasa, SI, Aug. 1992 (Made planar curve even if
*              start and end tangents are parallell. I also checked if we
*              in fact have circular conditions and itype != 1, to prevent
*              a singular matrix equation in s1616 I set itype = 1 in such
*              a case).
* Revised by : Christophe Rene Birkeland, SINTEF Oslo, May 1993.
* Revised by : Paal Fugelli, SINTEF, Oslo, Norway, Sept. 1994.  Fixed memory
*              leak from 'rc'.
*
*********************************************************************
*/
{
  int ki;		      /* Loop control parameter			*/
  int plane;		      /* Shape indicator for curve:
				 = 1  : Planar curve
				 else : Not planar			*/
  int kstat = 0;	      /* Status variable	 		*/
  int kpos = 0;		      /* Position of error  			*/
  int left;		      /* Array index used in call s1221 & s1227	*/
  int iopen;		      /* Used as a boolean variable:
			       * TRUE (=1) : open curve
			       * FALSE(=0) : closed curve		*/
  int knpoin;		      /* # of points to be used in interpolation */
  int knbpar;		      /* Number of different parameter values
			       * (Not used in this routine)		*/
  double sum1, sum2, sum3;    /* Used in normalization of vector algorithm  */
  double dum;
  double limit, dist;
  double tstpar;	      /* Parameter value to be used at the start
			       * of the curve (routine s1611)		*/
  double cndpar;	      /* Parameter value used at the end of the curve
			       * (Routine s1611)			*/
  double *ipar = SISL_NULL;	      /* Array containing the parameter values of the
			       * points in the curve. (not used)	*/
  double *kpoint1 = SISL_NULL;     /* Contains position and tangen of curve1	*/
  double *kpoint2 = SISL_NULL;     /* Contains position and tangen of curve2	*/
  double *type = SISL_NULL;	      /* Array (length inbpnt)containing type indicator
			       * for points/tangents :
			       *     1 - Ordinary point.
			       *     2 - Knuckle point. (Is treated as an
			       *	   ordinary point.)
			       *     3 - Tangent to next point.
			       *     4 - Tangent to prior point.	*/
  double *sdum = SISL_NULL;	      /* Is used to store cross product of tangent
			       * vectors				*/
  double tresol;	      /* Relative resulution of double precision
			       * numbers				*/
  double norder1[3];
  double norder2[3];
  double tcos, sang;
  int kstat1, kstat2;

  *jstat = 0;

  /* Check if curves are  correct. */

  s1707 (pc1, &kstat);
  if (kstat < 0) goto error;

  s1707 (pc2, &kstat);
  if (kstat < 0) goto error;

  /* Check if conflicting dimensions. */

  if (pc1->idim != pc2->idim) goto err106;
  if (pc1->idim != idim) goto err106;

  /* Perform filleting. */

  tresol = REL_COMP_RES;

  /* Calculate position and tangent of fillet point on curve1. */

  if((kpoint1 = newarray (2 * idim, DOUBLE)) == SISL_NULL) goto err101;

  left = 0;
  if (aend1 >= afil1)
    {
      /* Right hand position and tangent is to be calculated. */

      s1221 (pc1, 1, afil1, &left, kpoint1, &kstat);
      if (kstat < 0) goto error;

      /* Turn direction of tangent. */

      for (ki = idim; ki < 2 * idim; ki++)
	kpoint1[ki] = -kpoint1[ki];
    }
  else
    {
      /* Left hand position and tangent is to be calculated. */

      s1227 (pc1, 1, afil1, &left, kpoint1, &kstat);
      if (kstat < 0) goto error;
    }

  /* Calculate position and tangent of fillet point on curve2. */

  if((kpoint2 = newarray (2 * idim, DOUBLE)) == SISL_NULL) goto err101;

  left = 0;
  if (aend2 >= afil2)
    {
      /* Right hand position and tangent is to be calculated. */

      s1221 (pc2, 1, afil2, &left, kpoint2, &kstat);
      if (kstat < 0) goto error;
    }
  else
    {
      /* Left hand position and tangent is to be calculated. */

      s1227 (pc2, 1, afil2, &left, kpoint2, &kstat);
      if (kstat < 0) goto error;

      /* Turn direction of tangent. */

      for (ki = idim; ki < 2 * idim; ki++)
	kpoint2[ki] = -kpoint2[ki];
    }

  /* Find if conic arc choice will result in too refined B-spline. */

  if (itype == 2)
    {
      limit = (double) 0.0001;

      s6norm (&kpoint1[idim], idim, norder1, &kstat);
      s6norm (&kpoint2[idim], idim, norder2, &kstat);
      dist = s6dist (kpoint1, kpoint2, idim);

      for (ki = 0; ki < idim; ki++)
	{
	  norder1[ki] = kpoint1[ki] + dist * norder1[ki];
	  norder2[ki] = kpoint2[ki] - dist * norder2[ki];
	}
      if (aepsge / s6dist (norder1, norder2, idim) < limit)
	itype = 0;
    }

  /* Set type indicator. */

  type = newarray (4, DOUBLE);
  if (type == SISL_NULL) goto err101;

  type[0] = 1;
  type[1] = 4;
  type[2] = 1;
  type[3] = 4;

  /* Test if the points and tangents describe a planar curve
   * if we have 3-D curves. */

  plane = 0;
  if (idim == 3)
    {
      if((sdum = newarray (3, DOUBLE)) == SISL_NULL) goto err101;

      /* Compute cross product between tangent vectors. */

      s6crss (&kpoint1[idim], &kpoint2[idim], sdum);

      /* Normalize tangents and normal vector. */

      sum1 = (double) 0.0;
      sum2 = (double) 0.0;
      sum3 = (double) 0.0;

      for (ki = 0; ki < 3; ki++)
	{
	  dum = kpoint1[idim + ki];
	  sum1 += dum * dum;
	  dum = kpoint2[idim + ki];
	  sum2 += dum * dum;
	  sum3 += sdum[ki] * sdum[ki];
	}
      sum1 = sqrt (sum1);
      sum2 = sqrt (sum2);
      sum3 = sqrt (sum3);
      for (ki = 0; ki < 3; ki++)
	{
	  if (sum1 > (double) 0.0)
	    kpoint1[idim + ki] /= sum1;
	  else
	    kpoint1[idim + ki] = 0.0;
	  if (sum2 > (double) 0.0)
	    kpoint2[idim + ki] /= sum2;
	  else
	    kpoint2[idim + ki] = 0.0;
	  if (sum3 > (double) 0.0)
	    sdum[ki] /= sum3;
	  else
	    sdum[ki] = 0.0;
	}

      /* Find deviation from plane */

      dum = (double) 0.0;
      for (ki = 0; ki < 3; ki++)
	dum += (kpoint1[ki] - kpoint2[ki]) * sdum[ki];
      if (fabs (dum) < aepsge)
	plane = 1;
    }

  /* Test if a conic in fact is a circle. */

  if ((plane == 1 || idim == 2) && itype != 1)
     {
        dum = s6scpr(&kpoint1[idim], &kpoint2[idim], idim);

        sum1 = s6length(&kpoint1[idim], idim, &kstat1);
        sum2 = s6length(&kpoint2[idim], idim, &kstat2);

        if (!kstat1 || !kstat2)
          sang = DZERO;
        else
          {
            tcos = dum/(sum1*sum2);
            tcos = MIN((double)1.0,tcos);
            sang = acos(tcos);
          }

        if ((PI - sang) < ANGULAR_TOLERANCE)
	   {
	      /* Opposite parallell start and end tangents. */

	      for (ki = 0; ki < idim; ki++)
		 norder1[ki] = kpoint2[ki] - kpoint1[ki];
	      if (fabs(s6ang(norder1, &kpoint1[idim], idim) - PIHALF)
		  < ANGULAR_TOLERANCE)
		 itype = 1;
	   }
	else if (sang > ANGULAR_TOLERANCE)
	   {
	      /* Not parallell start and end tangents. */

	      if (idim == 2)
		 {
		    norder1[0] = kpoint1[idim + 1];
		    norder1[1] = -kpoint1[idim];
		    norder2[0] = kpoint2[idim + 1];
		    norder2[1] = -kpoint2[idim];
		    sum1 = (double) 0.0;
		    sum2 = (double) 0.0;
		    for (ki = 0; ki < idim; ki++)
		    {
		       sum1 += norder1[ki]*norder1[ki];
		       sum2 += norder2[ki]*norder2[ki];
		    }
		    sum1 = sqrt(sum1);
		    sum2 = sqrt(sum2);
		    for (ki = 0; ki < idim; ki++)
		    {
		       norder1[ki] /= sum1;
		       norder2[ki] /= sum2;
		    }
		 }
	      else
		 {
		    s6crss(sdum, &kpoint1[idim], norder1);
		    s6crss(sdum, &kpoint2[idim], norder2);
		 }

	      for (ki = 0; ki < idim; ki++)
	      {
		 norder1[ki] -= norder2[ki];
		 norder2[ki] = kpoint2[ki] - kpoint1[ki];
	      }
              if (s6ang(norder1, norder2, idim) < ANGULAR_TOLERANCE)
		 itype = 1;
	   }
     }

  /* Copy kpoint1 and kpoint2 into one array only kpoint1 */

  kpoint1 = increasearray (kpoint1, idim * 4, DOUBLE);
  if (kpoint1 == SISL_NULL)
    goto err101;

  memcopy (&kpoint1[2 * idim], kpoint2, 2 * idim, DOUBLE);

  if ((idim == 2 || (idim == 3 && plane == 1)) && (itype == 2 || itype == 1))
    {
      /* The curve is planar, interpolate with a conic. */

      if (itype == 1)
	knpoin = 3;
      else
	knpoin = 4;
      tstpar = (double) 0.0;
      iopen = SISL_CRV_OPEN;

      s1611 (kpoint1, knpoin, idim, type, iopen, ik, tstpar, aepsge, &cndpar,
	     rc, &kstat);
      if (kstat < 0 && kstat != -104)
	goto error;

      if ((*rc)->in >(*rc)->ik)
	goto out;

      if (kstat == 105)
	{
	  knpoin = 3;
	  tstpar = (double) 0.0;

	  if (*rc)  freeCurve(*rc);  /* Re-use of 'rc' (PFU 05/09-94) */
	  *rc = SISL_NULL;

	  s1611 (kpoint1, knpoin, idim, type, iopen, ik, tstpar, aepsge,
		 &cndpar, rc, &kstat);
	  if (kstat != -104 && (*rc)->in >(*rc)->ik)
	    {
	      if (kstat < 0)
		goto error;

	      /* Fillet produced */

	      goto out;
	    }
	}
    }

  /* Produce polynomial fillet */

  knpoin = 4;
  iopen = SISL_CRV_OPEN;
  tstpar = (double) 0.0;

  if (*rc)  freeCurve(*rc);  /* Re-use of 'rc' (PFU 05/09-94) */
  *rc = SISL_NULL;


  s1334 (kpoint1, knpoin, idim, type, 0, 0, iopen, ik, tstpar, &cndpar,
	 rc, &ipar, &knbpar, &kstat);
  if (ipar != SISL_NULL) freearray (ipar);
  if (kstat < 0) goto error;

  goto out;

  /* Error in input, conflicting dimensions */

  err106:
    *jstat = -106;
    s6err ("s1607", *jstat, kpos);
    goto out;

  /* Error in memory allocation */

  err101:
    *jstat = -101;
    s6err ("s1607", *jstat, kpos);
    goto out;

  /* Error in lower level function */

  error:
    *jstat = kstat;
    s6err ("s1607", *jstat, kpos);
    goto out;

  out:
    if (kpoint1 != SISL_NULL) freearray(kpoint1);
    if (kpoint2 != SISL_NULL) freearray(kpoint2);
    if (type != SISL_NULL) freearray(type);
    if (sdum != SISL_NULL) freearray(sdum);
    return;
}
