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
 * $Id: sh1780_at.c,v 1.2 2005-02-28 09:04:50 afr Exp $
 *
 */


#define SH1780_AT

#include "sislP.h"


#if defined(SISLNEEDPROTOTYPES)
void
sh1780_at (SISLObject * po1, SISLObject * po2,
	   SISLIntpt * pintpt, int *jstat)
#else
void
sh1780_at (po1, po2, pintpt, jstat)
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
* PURPOSE    : Set pre-topology data AT in curve-curve
*              intersection.
*
*
* INPUT      : po1      - Pointer to the first object in the intersection.
*              po2      - Pointer to the second object in the intersection.
*              pintpt   - Current intersection point.
*
*
* OUTPUT     : pintpt   - Intersection point with updated pre-topology info.
*              jstat    - status messages
*                                > 0   : Warning.
*                                = 0   : Ok.
*                                < 0   : Error.
*
*
* METHOD     :
*
* CALLS      : sh6gettop  -  Set topology of intersection point.
*              sh6settop  -  Set topology of intersection point.
*
* REFERENCES :
*
* WRITTEN BY : Ulf J. Krystad, SI, 09.91
*********************************************************************
*/
{
  int kstat = 0;		/* Status variable.                        */
  int kk1, kk2;			/* Orders of the two curves.               */
  int kn1, kn2;			/* Number of vertices in the curves.       */
  int lleft[2];			/* Array storing pre-topology information. */
  int lright[2];		/* Array storing pre-topology information. */
  double tref;			/* Reference value in equality test.       */
  double *st1, *st2;		/* Pointers to knot vectors of curves.     */
  double *sptpar = pintpt->epar;/* Parameter array of int.pt.              */
  /* --------------------------------------------------------------------- */

  /* Don't make pretop for help points ! */
  if (sh6ishelp (pintpt))
    {
      *jstat = 0;
      goto out;
    }


  /* Express the curve by local parameters.  */

  kn1 = po1->c1->in;
  kk1 = po1->c1->ik;
  st1 = po1->c1->et;
  kn2 = po2->c1->in;
  kk2 = po2->c1->ik;
  st2 = po2->c1->et;
  tref = MAX (st1[kn1] - st1[kk1 - 1], st2[kn2] - st2[kk2 - 1]);

  /* Update pre-topology of intersection point.  */
  sh6gettop (pintpt, -1, lleft, lright, lleft + 1, lright + 1, &kstat);

  /* Change the pre-topology information if the intersection point
	 lies at an endpoint of the curves.    */
  if (DEQUAL (sptpar[0] + tref, st1[kn1] + tref))
    {
      lright[0] = SI_AT;
    }
  if (DEQUAL (sptpar[0] + tref, st1[kk1 - 1] + tref))
    {
      lleft[0] = SI_AT;
    }
  if (DEQUAL (sptpar[1] + tref, st2[kn2] + tref))
    {
      lright[1] = SI_AT;
    }
  if (DEQUAL (sptpar[1] + tref, st2[kk2 - 1] + tref))
    {
      lleft[1] = SI_AT;
    }

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
