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
 * $Id: shdivsurf.c,v 1.3 2005-02-28 09:04:50 afr Exp $
 *
 */


#define SH_DIV_SURF

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
    sh_div_surf (SISLSurf * ps, int which_end_1, int which_end_2,
		 double aepsge, SISLSurf ** rsnew, int *jstat)
#else
void 
   sh_div_surf (ps, which_end_1, which_end_2, aepsge, rsnew, jstat)
   SISLSurf *ps;
   int which_end_1;
   int which_end_2;
   double aepsge;
   SISLSurf **rsnew;
   int *jstat;
#endif
/*
********************************************************************
*
*********************************************************************
*
* PURPOSE     :To factorize a bezier surface S(u,v)  over the interval 
*              [a,b]x[c,d] into
*              
*              oldsurf = (u-a)/(b-a) *newsurf 
*              when which_end_1 eq 0 
*                 and
*              oldcurve = (b-u)/(b-a)*newsurf 
*              when which_end_1 eq 1
*                 and
*              oldsurf = (v-c)/(d-c) *newsurf 
*              when which_end_2 eq 0 
*                 and
*              oldcurve = (d-v)/(d-c)*newsurf 
*              when which_end_2 eq 1
*                  NB ! Both which_end_1 and which_end_2 can be
*                       used at the same time. To achieve no factorisation
*                       in one dir, set the proper which_end to -1.
*
*                ------------------------
*                |          1            |
*                |                       |
*             ^  |                       |
*             |  | 0                   1 |
*             |  |                       |
*             v  |          0            |
*                ------------------------
*                    u --->
*              The edge in question must be zero.
*
*
* INPUT      : ps           - Oldsurf to factorize.
*              which_end_1   - (-1,0,1) Branch parameter for zero edge 1. dir.
*              which_end_2   - (-1,0,1) Branch parameter for zero edge 2. dir.
*              aepsge       - Geometry tolerance.
*
*
*
* OUTPUT     : rsnew      -The new surf.
*              jstat     - status messages
*                                         > 0      : warning
*                                         = 0      : ok
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
* WRITTEN BY : Ulf J. Krystad, SI, 92-12.
* MODIFIED BY :
*
**********************************************************************/
{
  int kstat;			/* Local status variable.		*/
  int kdim = ps->idim;		/* Dimension of geometry space.        */
  int kkind = ps->ikind;	/* Kind of surface.                    */
  int kk1;			/* Order in 1. par. dir.               */
  int kk2;			/* Order in 2. par. dir.               */
  int kn1;			/* Number of vertices in 1. par. dir.  */
  int kn2;			/* Number of vertices in 2. par. dir.  */
  double *st1;			/* Knot vector in 1. par. dir.         */
  double *st2;			/* Knot vector in 2. par. dir.         */
  double *scoef1 = SISL_NULL;	/* Coefficients of input curve to
			           factorize in 1. par. dir.           */
  double *scoef2 = SISL_NULL;	/* Coefficients of factorized surface. */
  double *scoef;		/* Coefficients of factorized surface. */
  SISLCurve *qc1 = SISL_NULL;	/* Input curve to sh_div_crv in 1. par. dir.    */
  SISLCurve *qc2 = SISL_NULL;	/* Output curve from sh_div_crv in 1. par. dir. */
  SISLCurve *qc3 = SISL_NULL;	/* Output curve from sh_div_crv in 2. par. dir. */
  /* __________________________________________________________________ */

  if (which_end_1 > -1)
    {
      /* Factorize 1. dir,
	 first express the surface as a curve.  */

      if ((scoef1 = newarray (kdim * ps->in1 * ps->in2, double)) == SISL_NULL)
	goto err101;

      /* Change parameter directions of surface.  */

      s6chpar (ps->ecoef, ps->in1, ps->in2, kdim, scoef1);

      /* Create curve.  */

      qc1 = newCurve (ps->in1, ps->ik1, ps->et1, scoef1, kkind, kdim * ps->in2, 0);
      if (qc1 == SISL_NULL)
	goto err101;

      /* Factorize the curve.  */
      sh_div_crv (qc1, which_end_1, aepsge, &qc2, &kstat);
      if (kstat < 0)
	goto error;

      /* Change parameter directions of the coefficient array of
         the resulting curve.  */

      if ((scoef2 = newarray (qc2->in *ps->in2 * kdim, DOUBLE)) == SISL_NULL)
	goto err101;
      s6chpar (qc2->ecoef, ps->in2, qc2->in, kdim, scoef2);

      /* Set local parameters of factorized surface. */

      kk1 = qc2->ik;
      kn1 = qc2->in;
      kk2 = ps->ik2;
      kn2 = ps->in2;
      st1 = qc2->et;
      st2 = ps->et2;

      /* Free curve used as input to the sh_div_crv. */

      if (qc1 != SISL_NULL)
	freeCurve (qc1);
      qc1 = SISL_NULL;
    }
  else
    {
      /* Set local parameters of input surface. */

      kk1 = ps->ik1;
      kk2 = ps->ik2;
      kn1 = ps->in1;
      kn2 = ps->in2;
      st1 = ps->et1;
      st2 = ps->et2;
      scoef2 = ps->ecoef;
    }

  if (which_end_2 > -1)
    {
      /* Factorize in second parameter direction of the
	 surface. First express the surface as a curve.           */

      if ((qc1 = newCurve (kn2, ps->ik2, st2, scoef2, kkind, kn1 * kdim, 0))
	  == SISL_NULL)
	goto err101;



      
      /* Factorize the curve.  */
      sh_div_crv(qc1, which_end_2, aepsge, &qc3, &kstat);
      if (kstat < 0)
	goto error;

      /*	Set local parameters of the surface. */

      kk2 = qc3->ik;
      kn2 = qc3->in;
      st2 = qc3->et;
      scoef = qc3->ecoef;
    }
  else
    scoef = scoef2;

  /* Express result as a surface.  */

  if ((*rsnew = newSurf (kn1, kn2, kk1, kk2, st1, st2,
			 scoef, kkind, kdim, 1)) == SISL_NULL)
    goto err101;

  /* Exit.  */

  *jstat = 0;
  goto out;

  /* Error in scratch allocation.  */

err101:*jstat = -101;
  goto out;

  /* Error in lower level routine.  */

error:*jstat = kstat;
  goto out;

out:
  /* Free scratch occupied by local arrays and objects.  */

  if (which_end_1 > -1 && scoef1 != SISL_NULL)
    freearray (scoef1);
  if (which_end_1 > -1&& scoef2 != SISL_NULL)
    freearray (scoef2);
  if (qc1 != SISL_NULL)
    freeCurve (qc1);
  if (qc2 != SISL_NULL)
    freeCurve (qc2);
  if (qc3 != SISL_NULL)
    freeCurve (qc3);

  return;
}
