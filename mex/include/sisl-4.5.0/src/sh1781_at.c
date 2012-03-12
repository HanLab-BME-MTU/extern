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
 * $Id: sh1781_at.c,v 1.2 2005-02-28 09:04:50 afr Exp $
 *
 */


#define SH1781_AT

#include "sislP.h"



#if defined(SISLNEEDPROTOTYPES)
void
sh1781_at (SISLObject * po1, SISLObject * po2,
	   SISLIntpt * pintpt,
	   int *jstat)
#else
void
sh1781_at (po1, po2, pintpt, jstat)
     SISLObject *po1;
     SISLObject *po2;
     SISLIntpt *pintpt;
     int *jstat;
#endif
/*
******************************************************************
*
*********************************************************************
*
* PURPOSE    : Set pre-topology data AT in 1-dimensional curve-point
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
* CALLS      : sh6gettop - Get topological info
*              sh6settop - Set topological info
* REFERENCES :
*
* WRITTEN BY : Ulf J. Krystad, SI, 06.91.
*********************************************************************
*/
{
  int kstat = 0;		/* Status variable.                        */
  int kn;			/* Number of vertices of curve.            */
  int kk;			/* Order of curve.                         */
  int lleft[2];			/* Array storing pre-topology information. */
  int lright[2];		/* Array storing pre-topology information. */
  int *ll1, *ll2, *lr1, *lr2;	/* Pointers into pre-topology arrays.      */
  double *st;			/* Pointer to knot vector of curve.        */
  double *sptpar = pintpt->epar;/* Pointer to parameter array of int.pt.   */
  double tref;			/* Referance value in equality test.       */
  SISLCurve *qc;		/* Pointer to current curve.               */
  /* --------------------------------------------------------------------- */


  /* Don't make pretop for help points ! */
  if (sh6ishelp (pintpt))
    {
      *jstat = 0;
      goto out;
    }

  /* Set pointers into the arrays storing pre-topology information. */

  if (po1->iobj == SISLCURVE)
    {
      ll1 = lleft;
      lr1 = lright;
      ll2 = lleft + 1;
      lr2 = lright + 1;
    }
  else
    {
      ll1 = lleft + 1;
      lr1 = lright + 1;
      ll2 = lleft;
      lr2 = lright;
    }

  /* Get pre-topology information. */
  sh6gettop (pintpt, -1, lleft, lright, lleft + 1, lright + 1, &kstat);
  if (kstat < 0)
    goto error;

  /* Test dimension of geometry space. */
  if (po1->iobj == SISLCURVE)
    qc = po1->c1;
  else
    qc = po2->c1;

  /* Store curve information in local parameters. */
  kn = qc->in;
  kk = qc->ik;
  st = qc->et;
  tref = st[kn] - st[kk - 1];

  /* Test if the intersection point lies at an endpoint of
     the curve. */

  if (DEQUAL (sptpar[0] + tref, st[kn] + tref))
    *lr1 = SI_AT;
  if (DEQUAL (sptpar[0] + tref, st[kk - 1] + tref))
    *ll1 = SI_AT;

  /* Update pretopology of intersection point.  */

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
