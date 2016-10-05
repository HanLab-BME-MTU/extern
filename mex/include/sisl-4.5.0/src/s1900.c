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
 * $Id: s1900.c,v 1.2 2001-03-19 15:58:55 afr Exp $
 *
 */


#define S1900

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
   s1900 (double param[], double knots[], double econd[], int ntype[],
       int inpt, int ik, int idim, int iopen,
       double *cendpar, SISLCurve ** rcurve,
       double **gpar, int *jnbpar, int *jstat)
#else
void
s1900 (param, knots, econd, ntype, inpt, ik, idim, iopen,
       cendpar, rcurve, gpar, jnbpar, jstat)
     double param[];
     double knots[];
     double econd[];
     int ntype[];
     int inpt;
     int ik;
     int idim;
     int iopen;
     double *cendpar;
     SISLCurve **rcurve;
     double **gpar;
     int *jnbpar;
     int *jstat;
#endif
/*
*********************************************************************
*
* PURPOSE    : Compute a B-spline curve interpolating a set of points.
*              The points can be assigned derivative conditions. The
*              curve can be open, closed, or closed and periodic.
*
* INPUT      : param - The parametrization of the point set.
*	       knots - The knot vector for the point set.
*              econd  - Array of interpolation conditions. Dimension
*                       is inpt*idim.
*              ntype  - Array containing kind of condition. Dimension
*                       is inpt.
*                       =  0 : A point is given.
*                       =  d : The d'th derivatative condition to the
*                              previous point is given.
*                       = -d : The d'th derivatative condition to the
*                              next point is given.
*              inpt   - Number of interpolation conditions.
*              ik     - Order of interpolating curve.
*              idim   - Dimension of geometry space.
*              iopen - Indicates if the curve is to be open, closed or
*                       periodic.
*
* OUTPUT     : cendpar - End parameter of parametrization.
*              rcurve  - Interpolating curve.
*	       gpar    - The distinct parameter values.
*	       jnbpar  - Number of distinct parameter values.
*              jstat   - status messages
*                        = 1      : Specified parametrization method
*                                   replaced by cord length parametrization.
*                                         = 0      : ok
*                                         < 0      : error
*
* METHOD     :
*
* REFERENCES :
*
* CALLS      : s1908, s1891, s1713, s1750
*
* WRITTEN BY : Vibeke Skytt, SI, 91-04.
* REVISED BY : Trond Vidar Stensby, SI, 91-07
*
*********************************************************************
*/
{
  int kstat = 0;		/* Status variable.                             */
  int kpos = 0;
  int ki;			/* Counter.                                     */
  int knpt;			/* Number of accepted interpolation conditions. */
  int kn;			/* Number of coefficients of B-spline curve.    */
  int kordr;			/* Local order of curve.                        */
  int kright = 1;		/* One equation system to solve in interpolation. */
  int knlr = 0;			/* Indicates shape of interpolation matrix.     */
  int knrc = 0;			/* Indicates shape of interpolation matrix.     */
  double *scoef = SISL_NULL;		/* Coefficients of curve.                          */
  int *ltype = SISL_NULL;		/* Type of accepted interpolation conditions.   */
  double *scond = SISL_NULL;		/* Array containing interpolation conditions.   */
  double *lpar = SISL_NULL;		/* Array containing new parameter valued. */
  int *sder = SISL_NULL;		/* Vector of derivative indicators.                */
  SISLCurve *qc = SISL_NULL;		/* Interpolation curve.                */
  SISLCurve *qc2 = SISL_NULL;	/* Interpolation curve. */

  *jstat = 0;

  /* Test interpolation conditions */

  s1908 (econd, ntype, param, inpt, ik, idim, iopen,
	 &scond, &ltype, &lpar, &knpt, &kstat);
  if (kstat < 0) goto error;

  /* Allocate scratch for derivative indicator. */

  if ((sder = newarray (knpt, INT)) == SISL_NULL)
    goto err101;

  for (ki = 0; ki < knpt; ki++)
    sder[ki] = (int) fabs ((double) ltype[ki]);

  /* Set local order.  */

  kordr = MIN (ik, knpt);

  if (iopen != SISL_CRV_OPEN)
    {
      knlr = kordr / 2;
      knrc = kordr - knlr - 1;
      knpt--;
    }
  /* Perform interpolation.  */

  s1891 (lpar, scond, idim, knpt, kright, sder, iopen, knots,
	 &scoef, &kn, kordr, knlr, knrc, &kstat);
  if (kstat < 0) goto error;

  /* Express the curve as a curve object.  */

  qc = newCurve (kn, kordr, knots, scoef, 1, idim, 1);
  if (qc == SISL_NULL) goto err101;

  if (!(iopen == SISL_CRV_OPEN))
    {
      /* A closed, non-periodic curve is expected. Pick the part of the
	 interpolation curve that has got a full basis.  */

      s1713 (qc, knots[kordr - 1], knots[kn], &qc2, &kstat);
      if (kstat < 0) goto error;

      if (qc != SISL_NULL) freeCurve (qc);
      qc = qc2;
    }

  if (kordr < ik)
    {
      /* The order of the curve is less than expected. Increase the order. */

      qc2 = SISL_NULL;
      s1750 (qc, ik, &qc2, &kstat);
      if (kstat < 0) goto error;
      if (qc != SISL_NULL) freeCurve (qc);
      qc = qc2;
    }

  /* Set open/closed parameter of curve.  */

  qc->cuopen = iopen;

  /* Interpolation performed. */

  /* Find distinct parameter values. */

  *gpar = lpar;

  *jnbpar = 1;
  for (ki = 1; lpar[ki] < *cendpar; ki++)
    {
      if (lpar[ki - 1] < lpar[ki])
	(*gpar)[(*jnbpar)++] = lpar[ki];
    }
  (*gpar)[(*jnbpar)++] = lpar[ki];

  *gpar = increasearray (*gpar, *jnbpar, DOUBLE);

  *rcurve = qc;
  goto out;

  /* Error in scratch allocation.  */

  err101:
    *jstat = -101;
    s6err ("s1900", *jstat, kpos);
    goto out;

  /* Error in lower level routine. */

  error:
    *jstat = kstat;
    s6err ("s1900", *jstat, kpos);
    goto out;

out:
  /* Free scratch occupied by local arrays. */

  if (scond != SISL_NULL)    freearray (scond);
  if (scoef != SISL_NULL)    freearray (scoef);
  if (sder != SISL_NULL)     freearray (sder);
  if (ltype != SISL_NULL)    freearray (ltype);

  return;
}
