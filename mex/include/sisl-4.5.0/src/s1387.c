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
 * $Id: s1387.c,v 1.2 2001-03-19 15:58:48 afr Exp $
 *
 */


#define S1387

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1387(SISLSurf *ps,int ik1,int ik2,SISLSurf **rsnew,int *jstat)
#else
void s1387(ps,ik1,ik2,rsnew,jstat)
     SISLSurf *ps;
     int  ik1;
     int  ik2;
     SISLSurf **rsnew;
     int  *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : To express the B-spline surface as a B-spline surface
*              of higher orders.
*
*
*
* INPUT      : ps       - Surface 
*              ik1      - New order in first direction 
*              ik2      - New order in second direction
*
*
* OUTPUT     : rsnew    - The resulting order elevated surface
*              jstat    - status messages
*                                         = 1      : Input orders equal
*                                                    to surface orders.
*                                                    Pointer set to input
*                                                    surface.
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : We treat  the surface as a curve in both parameter directions
*              and utilize the curve order elevation method twice
*
*
* REFERENCES : Larry L. Schumaker, Spline Functions:Basic Theory. 
*              Carl de Boor, A Practial Guide to Spline.           
*
*-
* CALLS      : newCurve  - Allocate space for a new curve-object.
*              freeCurve - Free space occupied by given curve-object.
*
* WRITTEN BY : Tor Dokken, SI 1988-11
* REVISED BY : Johannes Kaasa, SI, May 1992 (Introduced NURBS)
* REVISED BY : Christophe Birkeland, SI, July 1992 (Parameters in call s1750)
*
**********************************************************************/
{
  SISLCurve *qc1 = SISL_NULL;	/* Temporary curve                          */
  SISLCurve *qc2 = SISL_NULL;	/* Temporary curve                          */
  SISLCurve *qc3 = SISL_NULL;	/* Temporary curve                          */
  SISLCurve *qc4 = SISL_NULL;	/* Temporary curve                          */
  int kk1;			/* Order in first parameter direction       */
  int kk2;			/* Order in second parameter direction      */
  int kn1;			/* NumberOrder in first parameter direction */
  int kn2;			/* Order in second parameter direction      */
  int kdim;			/* Dimension used in temporary calc         */
  int kstat = 0;		/* Local parameter value                    */
  int kpos = 0;			/* Position of error                        */
  double *st1 = SISL_NULL, *st2 = SISL_NULL;	/* Pointers to knot vectors                 */
  double *scoef = SISL_NULL;		/* Pointer to coefficients                  */
  int rdim;                     /* Potential rational dimension.            */
  double *rcoef;                /* Potential rational vertices.             */

  *jstat = 0;


  kk1 = ps->ik1;
  kk2 = ps->ik2;
  kn1 = ps->in1;
  kn2 = ps->in2;
  
  if (ps->ikind == 2 || ps->ikind == 4)
    {
      rdim = ps->idim + 1;
      rcoef = ps->rcoef;
    }
  else
    {
      rdim = ps->idim;
      rcoef = ps->ecoef;
    }
    
  if (ik1 < kk1 || ik2 < kk2)
    goto err183;

  if (ik1 == kk1 && ik2 == kk2)
    goto war01;


  /* Create curve representing the surface a a curve in the second parameter
     direction, copy input arrays. */

  kdim = (ps->in1) * rdim;

  qc1 = newCurve (ps->in2, ps->ik2, ps->et2, rcoef, 1, kdim, 1);
  if (qc1 == SISL_NULL)
    goto err171;


  /* Make the order elevation in second parameter direction. */

  s1750 (qc1, ik2, &qc2, &kstat);
  if (kstat < 0)
    goto error;


  /* Remember new knot vector in second parameter direction. */

  kk2 = qc2->ik;
  kn2 = qc2->in;
  st2 = newarray (kk2 + kn2, DOUBLE);
  if (st2 == SISL_NULL)
    goto err101;

  memcopy (st2, qc2->et, kk2 + kn2, DOUBLE);

  /* Allocate space for turned parameter directions. */

  scoef = newarray (kn1 * kn2 * rdim, DOUBLE);
  if (scoef == SISL_NULL)
    goto err101;


  /* Turn parameter directions. */

  s6chpar (qc2->ecoef, kn1, kn2, rdim, scoef);


  /* Represent the surface a curve using the first knot vector. */

  kdim = kn2 * rdim;

  qc3 = newCurve (ps->in1, ps->ik1, ps->et1, scoef, 1, kdim, 1);
  if (qc3 == SISL_NULL)
    goto err101;


  /* Make the order elevation in the first parameter direction. */

  s1750 (qc3, ik1, &qc4, &kstat);
  if (kstat < 0)
    goto error;


  /* Remember new knot vector in first parameter direction. */

  kk1 = qc4->ik;
  kn1 = qc4->in;
  st1 = newarray (kk1 + kn1, DOUBLE);
  if (st1 == SISL_NULL)
    goto err101;

  memcopy (st1, qc4->et, kk1 + kn1, DOUBLE);

  /* Turn parameter directions of coefficients to match surface. */


  /* Allocate space for turned parameter directions. */

  scoef = increasearray (scoef, kn1 * kn2 * rdim, DOUBLE);
  if (scoef == SISL_NULL)
    goto err101;

  s6chpar (qc4->ecoef, kn2, kn1, rdim, scoef);


  /* Create surface object containing the order elevated surface. */

  *rsnew = newSurf (kn1, kn2, kk1, kk2, st1, st2, scoef, (ps->ikind), (ps->idim), 1);
  if (*rsnew == SISL_NULL)
    goto err171;
  
  /* Set periodicity flag according to that of the input surface. */
		      
  (*rsnew)->cuopen_1 = ps->cuopen_1;	      
  (*rsnew)->cuopen_2 = ps->cuopen_2;	      
    
  goto out;


  /* Input orders equal to surface orders */

war01:
  *jstat = 1;
  *rsnew = ps;
  goto out;

  /* Error in space allocation */

err101:
  *jstat = -101;
  s6err ("s1387", *jstat, kpos);
  goto out;

  /* Could not create curve or surface. */

err171:
  *jstat = -171;
  s6err ("s1387", *jstat, kpos);
  goto out;

  /* Order(s) specified too low */

err183:
  *jstat = -183;
  s6err ("s1387", *jstat, kpos);
  goto out;

  /* Error in lower level function. */

error:
  *jstat = kstat;
  s6err ("s1387", *jstat, kpos);
  goto out;

  /* Free local used memory. */

out:
  if (qc1 != SISL_NULL)
    freeCurve (qc1);
  if (qc2 != SISL_NULL)
    freeCurve (qc2);
  if (qc3 != SISL_NULL)
    freeCurve (qc3);
  if (qc4 != SISL_NULL)
    freeCurve (qc4);
  if (st1 != SISL_NULL)
    freearray (st1);
  if (st2 != SISL_NULL)
    freearray (st2);
  if (scoef != SISL_NULL)
    freearray (scoef);

  return;
}
