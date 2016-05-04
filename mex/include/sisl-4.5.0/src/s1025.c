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
 * $Id: s1025.c,v 1.3 2005-02-28 09:04:48 afr Exp $
 *
 */


#define S1025

#include "sislP.h"


#if defined(SISLNEEDPROTOTYPES)
void
s1025 (SISLSurf * ps, double epar1[], int inpar1, double epar2[],
	    int inpar2, SISLSurf ** rsnew, int *jstat)
#else
void 
s1025 (ps, epar1, inpar1, epar2, inpar2, rsnew, jstat)
     SISLSurf *ps;
     double epar1[];
     int inpar1;
     double epar2[];
     int inpar2;
     SISLSurf **rsnew;
     int *jstat;
#endif
/*
********************************************************************
*
*********************************************************************
*
* PURPOSE    : Insert a given set of knots in each parameter direction
*	       into the description of a B-spline surface.
* NOTE       : When the surf is periodic in one direction, the input
*              parameter values in this direction
*              must lie in the HALFOPEN [et[kk-1], et[kn), the function
*              will automatically update the extra knots and
*              coeffisients.
*
*
*
*
* INPUT      : ps	- Surface to be refined.
*	       epar1    - Knots to insert in 1. parameter direction.
*   	       inpar1   - Number of new knots in 1. par. dir.
*	       epar2    - Knots to insert in 2. parameter direction.
*   	       inpar2   - Number of new knots in 2. par. dir.
*
*
*
* OUTPUT     : rsnew	- The new, refined surface.
*              jstat	- status messages
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
* CALLS      : newSurf     - Allocate space for a new surface-object.
*	       newCurve    - Allocate space for a new curve-object.
*	       freeCurve   - Free space occupied by a curve-object.
*              s1018       - Insert knots in a B-spline curve.
*
* WRITTEN BY : Vibeke Skytt, SI, 01.91.
* CHANGED BY : Ulf J. Krystad, SI, 92-01
*              Treatment of periodic curves.
*
**********************************************************************/
{
  int kstat;			/* Local status variable.		*/
  int kdim = ps->idim;		/* Dimension of geometry space.        */
  int kkind = ps->ikind;	/* Kind of surface.                    */
  int kn1;			/* Number of vertices in 1. par. dir.  */
  int kn2;			/* Number of vertices in 2. par. dir.  */
  double *st1;			/* Knot vector in 1. par. dir.         */
  double *st2;			/* Knot vector in 2. par. dir.         */
  double *scoef1 = SISL_NULL;	/* Coefficients of input curve to
			           refinement in 1. par. dir.          */
  double *scoef2 = SISL_NULL;	/* Coefficients of refined surface.    */
  double *scoef;		/* Coefficients of refined surface.    */
  SISLCurve *qc1 = SISL_NULL;	/* Input curve to refinement in 1. par. dir.    */
  SISLCurve *qc2 = SISL_NULL;	/* Output curve from refinement in 1. par. dir. */
  SISLCurve *qc3 = SISL_NULL;	/* Output curve from refinement in 2. par. dir. */
  double *oldcoef;           	/* Pointer to vertices of old surf.                */


  
  if(kkind == 2 || kkind == 4)
  {
     oldcoef = ps->rcoef;
     kdim++;
  }
  else
  {
     oldcoef = ps->ecoef;
  }
  
  if (inpar1 > 0)
  {
     /* Insert knots in the first parameter direction of the
	 surface. First express the surface as a curve.  */

      if ((scoef1 = newarray (kdim * ps->in1 * ps->in2, double)) == SISL_NULL)
	goto err101;

      /* Change parameter directions of surface.  */

      s6chpar (oldcoef, ps->in1, ps->in2, kdim, scoef1);

      /* Create curve.  */

      qc1 = newCurve (ps->in1, ps->ik1, ps->et1, scoef1,1, kdim * ps->in2, 0);
      if (qc1 == SISL_NULL)
	goto err101;

      /* Set periodicity flag of curve.  */
	
      qc1->cuopen = ps->cuopen_1;
      
      /* Insert knots into this curve.  */

      s1018 (qc1, epar1, inpar1, &qc2, &kstat);
      if (kstat < 0)
	goto error;

      /* Change parameter directions of the coefficient array of
         the refined curve.     */

      if ((scoef2 = newarray (qc2->in *ps->in2 * kdim, DOUBLE)) == SISL_NULL)
	goto err101;
      s6chpar (qc2->ecoef, ps->in2, qc2->in, kdim, scoef2);

      /* Set local parameters of refined surface. */

      kn1 = qc2->in;
      kn2 = ps->in2;
      st1 = qc2->et;
      st2 = ps->et2;

      /* Free curve used as input to the knotinsertion for curves. */

      if (qc1 != SISL_NULL)
	freeCurve (qc1);
      qc1 = SISL_NULL;
    }
  else
    {
      /* Set local parameters of input surface. */

      kn1 = ps->in1;
      kn2 = ps->in2;
      st1 = ps->et1;
      st2 = ps->et2;
      scoef2 = oldcoef;
    }

  if (inpar2 > 0)
    {
      /* Insert knots into the second parameter direction of the
	 surface. First express the surface as a curve.           */

      if ((qc1 = newCurve (kn2, ps->ik2, st2, scoef2, 1, kn1 * kdim, 0))
	  == SISL_NULL)
	goto err101;

      /* Set periodicity flag of curve.  */
      
      qc1->cuopen = ps->cuopen_2;
      
      /* Insert knots into this curve.  */

      s1018 (qc1, epar2, inpar2, &qc3, &kstat);
      if (kstat < 0)
	goto error;

      /*	Set local parameters of the refined surface. */

      kn2 = qc3->in;
      st2 = qc3->et;
      scoef = qc3->ecoef;
    }
  else
    scoef = scoef2;

  /* Express result as a surface.  */

  if ((*rsnew = newSurf (kn1, kn2, ps->ik1, ps->ik2, st1, st2,
			 scoef, kkind, ps->idim, 1)) == SISL_NULL)
    goto err101;

  /* Copy periodicity flag from input surface.  */
  
  (*rsnew)->cuopen_1 = ps->cuopen_1;
  (*rsnew)->cuopen_2 = ps->cuopen_2;
  
  /* Refinement performed.  */

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

  if (inpar1 > 0 && scoef1 != SISL_NULL)
    freearray (scoef1);
  if (inpar1 > 0 && scoef2 != SISL_NULL)
    freearray (scoef2);
  if (qc1 != SISL_NULL)
    freeCurve (qc1);
  if (qc2 != SISL_NULL)
    freeCurve (qc2);
  if (qc3 != SISL_NULL)
    freeCurve (qc3);

  return;
}
