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
 * $Id: s1520.c,v 1.3 2001-03-19 15:58:50 afr Exp $
 *
 */


#define S1520

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1520 (SISLCurve * pc, double angle, double ep[], double eaxis[],
       SISLSurf ** rs, int *jstat)
#else
void
s1520 (pc, angle, ep, eaxis, rs, jstat)
     SISLCurve *pc;
     double angle;
     double ep[];
     double eaxis[];
     SISLSurf **rs;
     int *jstat;
#endif
/*
*********************************************************************
*
* PURPOSE    : To create a NURBS rotational surface by rotating
*              the curve *pc around the axis defined by ep[] and eaxis[]
*              the given angle. This will be an exact representation.
*
*
* INPUT      : pc     - Pointer to curve to be rotated (NURBS or B-spline).
*              angle  - The rotational angle. Counter clockwise around axis.
*                       If the absolute value of the angle is greater than
*                       2 PI then a rotational surface closed in the
*                       rotation direction is made.
*              ep     - SISLPoint on rotational axis
*              eaxis  - Direction of rotational axis
*
*
* OUTPUT     :
*              jstat  - status messages
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*              rs     - Pointer to the surface produced (NURBS or B-spline).
*
* METHOD     : First a normalized circle segment spanning the actual angle
*              is generated. This circle is then translated to generate
*              the actual rows of control vertices of the surface
*
* REFERENCES :
*
*-
* CALLS      : s1713, s6rotax, s6mvec, newCurve, newSurf, s6err
*
* WRITTEN BY : Johannes Kaasa, SI, Oslo, Norway. 09. Aug. 1991
*              Based on s1302 with use of NURBS instead of B-splines.
* Revised by : Christophe Rene Birkeland, SINTEF Oslo, May 1993.
*              SISL_NULL tests included
* Revised by : Paal Fugelli, SINTEF, Oslo, Norway, Dec. 1994.  Added
*              initialization and check for allocation error of 'rs'.
*              The 'cuopen_1' flag is now set to CLOSED if angle >= 2PI, it
*              should really have been made periodic - and might at a later
*              date be updated.
*
*********************************************************************
*/
{
  double *st1;			/* Pointer to knot vector of circle segment   */
  double *scoef1;		/* Pointer to vertices of circle segment      */
  double *rcoef1;		/* Pointer to rational vertices of circle seg */
  int kn1;			/* Number of vertice of circle segment        */
  int kk1;			/* Order of circle segment                    */
  double *st2;			/* Pointer to knot vector of input curve      */
  double *scoef2;		/* Pointer to vertices of input curve         */
  double *rcoef2;		/* Pointer to rational vertices of input curve */
  int kn2;			/* Number of vertice of input curve           */
  int kk2;			/* Order of input curve                       */
  int kdim;			/* Dimension of space in which curve lies     */
  int ksurfdim;			/* Dimension of space in which surface lies   */
  int ki;			/* Control variable in loop                   */
  int kj;			/* Control variable in loop                   */
  int kl;			/* Control variable in loop                   */
  int kn;			/* Number of vertices in total circle         */
  int kk;			/* Order of total circle                      */
  double *st;			/* Pointer to knots of total circle           */
  double *scoef;		/* Pointer to vertices of total circle        */
  int kind;			/* Kind of total circle curve                 */
  int kcopy;			/* Copy flag for total circle curve           */
  double tangle;		/* Local positive rotational angle            */
  int quadrant;			/* Quadrant of the rotational angle           */
  double resang;		/* Residue angle in the actual quadrant       */
  double spar;			/* Start parameter of circle segment          */
  double epar;			/* End parameter of circle segment            */
  SISLCurve *totcurve=SISL_NULL;	/* Pointer to total normalized circle         */
  SISLCurve *pnorm;		/* Pointer to normalized circle segment       */
  double weight = (double) 1. / sqrt ((double) 2.);	/* Rational weight                         */
  double tfac;			/* Weights along the profile curve.           */
  double *sucof = SISL_NULL;		/* Pointer to vertex array for surface        */
  double smat[16];		/* Transformation matrix                      */
  int kstat;			/* Status variable                            */
  double *srow;			/* Pointer to row of vertices in surface      */
  double *scirc;		/* Pointer to vertices in circular arc        */
  int nbvec = 1;		/* Number of vectors                          */

  int kpos = 1;			/* Position of error                          */



  /* Ensure a valid output surface. */

  *rs = SISL_NULL;


  /* Make local pointers to description of curve */

  st2      = pc->et;
  kn2      = pc->in;
  kk2      = pc->ik;
  scoef2   = pc->ecoef;
  rcoef2   = pc->rcoef;
  kdim     = pc->idim;
  ksurfdim = kdim + 1;

  /* The routine is only working for dimension=3 */

  if (kdim != 3)
    goto err104;

  /* Calculate normalized NURBS circle */

  kn = 9;
  kk = 3;
  st = newarray (kn + kk, DOUBLE);
  scoef = newarray (kn * ksurfdim, DOUBLE);
  kind = 4;
  kcopy = 2;
  st[0] = (double) 0.;
  for (ki = 1; ki < kk; ki++)
    {
      st[ki]     = (double) 0.;
      st[2 + ki] = PIHALF;
      st[4 + ki] = PI;
      st[6 + ki] = THREEPIHALF;
      st[8 + ki] = TWOPI;
    }
  st[11] = TWOPI;
  for (ki = 0; ki < 36; ki++)
    {
      switch (ki)
	{
	case 1:
	case 2:
	case 6:
	case 8:
	case 10:
	case 14:
	case 17:
	case 18:
	case 22:
	case 24:
	case 26:
	case 30:
	case 33:
	case 34:
	  scoef[ki] = (double) 0.0;
	  break;
	case 0:
	case 3:
	case 9:
	case 11:
	case 19:
	case 27:
	case 32:
	case 35:
	  scoef[ki] = (double) 1.0;
	  break;
	case 16:
	case 25:
	  scoef[ki] = (double) -1;
	  break;
	case 4:
	case 5:
	case 7:
	case 13:
	case 15:
	case 23:
	case 28:
	case 31:
	  scoef[ki] = weight;
	  break;
	case 12:
	case 20:
	case 21:
	case 29:
	  scoef[ki] = -weight;
	  break;
	}
    }

  totcurve = newCurve (kn, kk, st, scoef, kind, kdim, kcopy);
  if(totcurve == SISL_NULL) goto err101;

  /* Pick out a part of the total curve */

  tangle = fabs (angle);
  if (tangle > TWOPI)
    tangle = TWOPI;
  quadrant = (int) floor (tangle / PIHALF);
  resang = tangle - quadrant * PIHALF;
  spar = (double) 0.;
  epar = PIHALF * (quadrant + ((double) 1. + (sqrt ((double) 2.) + (double) 1.) * tan ((resang - PI / (double) 4.) / (double) 2.)) / (double) 2.);
  s1713 (totcurve, spar, epar, &pnorm, &kstat);
  if (kstat < 0)
    goto error;

  /* Make local variables for curve description */

  st1    = pnorm->et;
  scoef1 = pnorm->ecoef;
  rcoef1 = pnorm->rcoef;
  kn1    = pnorm->in;
  kk1    = pnorm->ik;

  if (angle < DZERO)
    for (ki = 0; ki < kn1; ki++)
      scoef1[kdim * ki + 1] = -scoef1[kdim * ki + 1];

  /* Allocate vertex array for surface */

  sucof = newarray (kn1 * kn2 * ksurfdim, DOUBLE);
  if (sucof == SISL_NULL)
    goto err101;

  /* Make the surface vertices circle segment by circle segment */

  for (ki = 0; ki < kn2; ki++)
    {

      /*  Check for profile weights */

      if (pc->ikind == 2 || pc->ikind == 4)
	tfac = rcoef2[(ki + 1) * ksurfdim - 1];
      else
	tfac = (double) 1.;

      /*  Make transformation matrix for first vertex on curve to be rotated */

      s6rotax (ep, eaxis, &scoef2[ki * kdim], smat, &kstat);
      if (kstat < 0)
	goto error;

      /*  Transform the vertices of this row into right position */

      for (kj = 0; kj < kn1; kj++)
	{
	  srow  = sucof + ki * (kn1 * ksurfdim) + kj * ksurfdim;
	  scirc = scoef1 + kj * kdim;
	  s6mvec (smat, scirc, nbvec, srow);
	  weight = rcoef1[(kj + 1) * ksurfdim - 1] * tfac;
	  for (kl = 0; kl < kdim; kl++)
	    srow[kl] *= weight;
	  srow[kdim] = weight;
	}
    }

  /* Create the surface */

  *rs = newSurf (kn1, kn2, kk1, kk2, st1, st2, sucof, 2, kdim, 1);
  if ( *rs == SISL_NULL )  goto err101;

  if ( tangle >= TWOPI )
  {
    /* Set the flag indicating that the surface is closed in the first
       parameter direction.  It should really have been made cyclic/periodic. */

    (*rs)->cuopen_1 = SISL_SURF_CLOSED;
  }

  *jstat = 0;
  goto out;

  /* Error in space allocation.  */

  err101:*jstat = -101;
    s6err ("s1520", *jstat, kpos);
    goto out;

  /* Error in input, dimension not equal to 3 */

  err104:*jstat = -104;
    s6err ("s1520", *jstat, kpos);
    goto out;

  /* Error in lower level routine.  */

  error:*jstat = kstat;
    s6err ("s1520", *jstat, kpos);
    goto out;

  out:
    /* Free allocated arrays */

    if (sucof != SISL_NULL) freearray (sucof);
    if (totcurve != SISL_NULL) freeCurve (totcurve);
    if (pnorm != SISL_NULL) freeCurve (pnorm);
    return;
}
