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
 * $Id: sh1779_at.c,v 1.2 2005-02-28 09:04:50 afr Exp $
 *
 */


#define SH1779_AT

#include "sislP.h"


#if defined(SISLNEEDPROTOTYPES)
void
sh1779_at (SISLObject * po1, SISLObject * po2, SISLIntpt * pintpt,
	   int *jstat)
#else
void
sh1779_at (po1, po2, pintpt, jstat)
     SISLObject *po1;
     SISLObject *po2;
     SISLIntpt *pintpt;
     int *jstat;
#endif


/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Set pre-topology AT in 3-dimensional curve-surface
*              intersection.
*
*
* INPUT      : po1      - Pointer to the first object in the intersection.
*              po2      - Pointer to the second object in the intersection.
*              pintpt   - Current intersection point.
*
*
* OUTPUT     : pintpt   - Intersection point after updating pre-topology.
*              jstat    - status messages
*                                > 0   : Warning.
*                                = 0   : Ok.
*                                < 0   : Error.
*
*
* METHOD     :
*
* CALLS      : sh6gettop  -  Get topology of intersection point.
*              sh6settop  -  Set topology of intersection point.
*
* REFERENCES :
*
* WRITTEN BY : Ulf J. Krystad, SI, 06.91.
*********************************************************************
*/
{
  int kstat = 0;		/* Status variable.                        */
  int kpar1, kpar2;		/* Index of parameter value of object.     */
  int kn;			/* Number of vertices of curve.            */
  int kk;			/* Order of curve.                         */
  int lleft[2];			/* Array storing pre-topology information. */
  int lright[2];		/* Array storing pre-topology information. */
  int *ll1, *ll2, *lr1, *lr2;	/* Pointers into pre-topology arrays.   */
  double tref;			/* Referance value in equality test.       */
  double *st;			/* Knot vector of curve.                   */
  double *sptpar = pintpt->epar;/* Pointer to parameter values of int.pt.  */
  SISLCurve *qc;		/* Pointer to the curve.                   */
  SISLSurf *qs;			/* Pointer to the surface.                 */
  double sf_low_lim[2];
  double sf_high_lim[2];
  /* ---------------------------------------------------------------------- */
  /* Don't make pretop for help points ! */
  if (sh6ishelp (pintpt))
    {
      *jstat = 0;
      goto out;
    }

  /* Set pointers into the arrays storing pre-topology information. */
  if (po1->iobj == SISLCURVE)
    {
      qc = po1->c1;
      qs = po2->s1;
      kpar1 = 0;
      kpar2 = 1;
      ll1 = lleft;
      lr1 = lright;
      ll2 = lleft + 1;
      lr2 = lright + 1;
    }
  else
    {
      qc = po2->c1;
      qs = po1->s1;

      kpar1 = 2;
      kpar2 = 0;
      ll1 = lleft + 1;
      lr1 = lright + 1;
      ll2 = lleft;
      lr2 = lright;
    }

  kk = qc->ik;
  kn = qc->in;
  st = qc->et;
  tref = st[kn] - st[kk - 1];

  sf_low_lim[0] = qs->et1[qs->ik1 - 1] + REL_COMP_RES;
  sf_low_lim[1] = qs->et2[qs->ik2 - 1] + REL_COMP_RES;
  sf_high_lim[0] = qs->et1[qs->in1] - REL_COMP_RES;
  sf_high_lim[1] = qs->et2[qs->in2] - REL_COMP_RES;

  sh6gettop (pintpt, -1, lleft, lright, lleft + 1, lright + 1, &kstat);
  if (kstat < 0)
    goto error;
  /* Check endpoint of curve. */

  if (DEQUAL (sptpar[kpar1] + tref, st[kk - 1] + tref))
    *ll1 = SI_AT;
  if (DEQUAL (sptpar[kpar1] + tref, st[kn] + tref))
    *lr1 = SI_AT;

  /* Update pre-topology of intersection point.  */
  sh6settop (pintpt, -1, lleft[0], lright[0], lleft[1], lright[1], &kstat);
  if (kstat < 0)
    goto error;

  *jstat = 0;
  goto out;

  /* Error lower level routine.  */
error:*jstat = kstat;
  goto out;

out:
  return;
}
