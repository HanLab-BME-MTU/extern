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
 * $Id: s1370.c,v 1.3 2001-03-19 15:58:47 afr Exp $
 *
 */


#define S1370

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1370 (SISLCurve * pcurv, double earray[], int idim, int inarr,
       int ratflag, SISLCurve ** rcurv, int *jstat)
#else
void
s1370 (pcurv, earray, idim, inarr, ratflag, rcurv, jstat)
     SISLCurve *pcurv;
     double earray[];
     int idim;
     int inarr;
     int ratflag;
     SISLCurve ** rcurv;
     int *jstat;
#endif
/*
*********************************************************************
*
* PURPOSE    : To put a curve description into the implicit
*              second order surface described by the input array.
*
* INPUT      : pcurv  - Pointer to input curve
*              earray - The description of the input array
*                       dimension (idim+1)x(idim+1) (xinarr)
*              idim   - Put curve into implicit equation
*              inarr  - Number of parallel matrices in earray.
*                       inarr should be less or equal to 3.
*              ratflag - If pcurv is nonrational it is ignored.
*                        Otherwise:
*                        If ratflag = 0 rcurv is the nonrational numerator
*                        If ratflag = 1 rcurv is a full rational curve
*
* OUTPUT     : rcurv  - The resulting curve
*              jstat  - status messages
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
* METHOD     : Dependent on the type of object we make:
*
*        F(S,T) = (P(s,t),1)  EARRAY (P(s,t),1)
*
*     by sampling enough point to use interpolation for reproduction.
*
* REFERENCES :
*
* CALLS      : s1893,s6err.
*
* WRITTEN BY : Tor Dokken, SI, Oslo , Norway
* REVISED BY : Mike Floater, SI, Oslo 11/4/91 for rational curves
* REVISED BY : Mike Floater, SI, Oslo 11/9/91 -- ratflag.
* REVISED BY : Michael Floater, SI, June 92. The rational stuff
*              was completely messed up after translation of
*              the fortran part to c. But it works now.
* REVISED BY : Christophe Birkeland, SI, July 1992.
* REVISED BY : Christophe Rene Birkeland, SINTEF Oslo, May 1993.
*              jcurve removed, other minor changes
* Revised by : Paal Fugelli, SINTEF, Oslo, Norway, September 1994.
*              Didn't work for rationals - '(*rcurve)->rcoef' returned from
*              s1893() was SISL_NULL (must copy from ecoef).
*
*********************************************************************
*/
{
  int kpos = 0;
  int kstat = 0;
  SISLCurve *icurve = SISL_NULL;	/* Temporary SISLCurve. */
  int kn;			/* Number of vertices of pcurv            */
  int kk;			/* Order in  pcurv                        */
  int kdim;			/* Number of dimesions in pcurv           */
  int kdimp1;			/* Dimension of  earray should be kdim+1  */
  double *st = SISL_NULL;		/* First knot vector is pcurv             */
  double *scoef = SISL_NULL;		/* Vertices of pcurv                      */
  int ikind;			/* kind of surface pcurv is               */
  double *rscoef = SISL_NULL;	/* Scaled coefficients if pcurv is rational       */
  double wmin, wmax;		/* min and max values of the weights if rational  */
  double scale;			/* factor for scaling weights if rational         */
  int i;			/* loop variable                          */
  double *sarray = SISL_NULL;	/* Array for calculating denominator if used      */
  int knarr;			/* Number of parallel arrays to use.              */
  int nkind;			/* Kind of output curve (rcurf).                  */

  *jstat = 0;

  /* Make local pointers. */

  kn = pcurv->in;
  kk = pcurv->ik;
  kdim = pcurv->idim;
  st = pcurv->et;
  ikind = pcurv->ikind;

  kdimp1 = kdim + 1;

  /* Test input. */

  if (kdim != idim || (kdim != 2 && kdim != 3))
    goto err104;
  if (inarr < 1 || 3 < inarr) goto err172;

  /* rational surfaces are a special case. */
  if (ikind == 2 || ikind == 4)
    {
      kdim++;

      /* scale the coeffs so that min. weight * max. weight = 1. */

      rscoef = pcurv->rcoef;
      wmin = rscoef[kdim-1];
      wmax = rscoef[kdim-1];

      for (i = 2*kdim-1; i < kn * kdim; i += kdim)
	{
	  if (rscoef[i] < wmin)
	    wmin = rscoef[i];
	  if (rscoef[i] > wmax)
	    wmax = rscoef[i];
	}
      scale = (double) 1.0 / sqrt (wmin * wmax);
      scoef = newarray (kn * kdim, DOUBLE);
      if (scoef == SISL_NULL)
	goto err101;

      for (i = 0; i < kn * kdim; i++)
        scoef[i] = rscoef[i] * scale;
    }
  else
    scoef = pcurv->ecoef;

  icurve = newCurve (kn, kk, st, scoef, 1, kdim, 1);
  if (icurve == SISL_NULL)
    goto err171;

  icurve->cuopen = pcurv->cuopen;

  if ((ikind == 2 || ikind == 4) && ratflag == 1)
    {
      /* Output curve will also be rational. */

      nkind = 2;

      /* Add an extra parallel array to pick up the weights
	 of the subsequent homogeneous vertices of rcurv. */

      knarr = inarr + 1;
      sarray = new0array (kdimp1 * kdimp1 * knarr, DOUBLE);
      if (sarray == SISL_NULL) goto err101;

      memcopy (sarray, earray, kdimp1 * kdimp1 * inarr, DOUBLE);
      sarray[kdimp1 * kdimp1 * knarr - 1] = (DOUBLE) 1.0;
    }
  else
    {
      nkind = 1;
      knarr = inarr;
      sarray = earray;
    }

  /* Put curve into implicit surface. */

  s1893 (icurve, sarray, kdimp1, knarr, 0, 0, rcurv, &kstat);
  if (kstat < 0) goto error;

  if (*rcurv == SISL_NULL) goto err171;

  if ( ikind == 2 || ikind == 4 )
  {
    /* Free arrays. */

    if (scoef) freearray (scoef);
    if (ratflag && sarray) freearray (sarray);

    if ( ratflag == 1 )
    {
      /* Output from s1893 is a dim+1 non-rational curve. */
      /* Convert homogeneous curve to rational form (rcoef is SISL_NULL here). */

      (*rcurv)->rcoef = newarray((*rcurv)->in * (*rcurv)->idim, DOUBLE);
      memcopy((*rcurv)->rcoef, (*rcurv)->ecoef,
	      (*rcurv)->in * (*rcurv)->idim, DOUBLE);

      (*rcurv)->idim --;    /* Adjust from the homogeneus coordinates. */
      (*rcurv)->ikind = 2;  /* i.e. rational */

    }
  }


  /* Ok ! */

  goto out;

  /* Error in lower level function. */

  error:
    *jstat = kstat;
    s6err ("s1370", *jstat, kpos);
    goto out;

  /* Allocation problems.    */

  err101:
    *jstat = -101;
    s6err ("s1370", *jstat, kpos);
    goto out;

  /* Dimension not equal to 3.    */

  err104:
    *jstat = -104;
    s6err ("s1370", *jstat, kpos);
    goto out;

  /* Could not create curve */

  err171:
    *jstat = -171;
    s6err ("s1370", *jstat, kpos);
    goto out;

  /* Dimension inarr not equal to 1,2 or 3. */

  err172:
    *jstat = -172;
    s6err ("s1370", *jstat, kpos);
    goto out;

  out:
  if (icurve != SISL_NULL) freeCurve (icurve);
  return;
}
