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
 * $Id: sh1786.c,v 1.2 2001-03-19 15:59:05 afr Exp $
 *
 */


#define SH1786

#include "sislP.h"


#if defined(SISLNEEDPROTOTYPES)
void
sh1786 (SISLObject * po1, SISLObject * po2, double aepsge,
	SISLIntdat ** rintdat, SISLIntpt * pintpt, int *jnewpt,
	int *jstat)
#else
void
sh1786 (po1, po2, aepsge, rintdat, pintpt, jnewpt, jstat)
     SISLObject *po1;
     SISLObject *po2;
     double aepsge;
     SISLIntdat **rintdat;
     SISLIntpt *pintpt;
     int *jnewpt;
     int *jstat;
#endif
/*
******************************************************************
*
*********************************************************************
*
* PURPOSE    : Set pre-topology data in n-dimensional curve-point
*              intersection. Also find help points if necessary.
*
*
* INPUT      : po1      - Pointer to the first object in the intersection.
*              po2      - Pointer to the second object in the intersection.
*              aepsge   - Geometry resolution.
*              rintdat  - Intersection data structure.
*              pintpt   - Current intersection point.
*
*
* OUTPUT     : pintpt   - Intersection point after updating pre-topology.
*              jnewpt   - Number of new int.pt. created.
*              jstat    - status messages
*                                > 0   : Warning.
*                                = 0   : Ok.
*                                < 0   : Error.
*
*
* METHOD     :
*
* CALLS      : shevalc  -  Evaluate curve using filtered coefficients.
*              s6ang    -  Angle between two vectors.
*              s6idcon  -  Connect two intersection points.
*              hp_newIntpt -  Create new intersection point.
*
* REFERENCES :
*
* WRITTEN BY : Michael Floater, SI, Oslo, Norway. 09.91.
* CHANGED BY : UJK, ALA, SI, 10.10.91:
*              Purpose of function is only to generate helppoints in 2D
* REWISED BY : VSK, 11-92.  Make sure that a help point always is an
*                           intersection point.
*********************************************************************
*/
{
  int kstat = 0;		/* Status variable.                        */
  int ki, kj;			/* Counters.                               */
  int kleft = 0;		/* Parameter to evaluator.                 */
  int korgleft = 0;		/* Knot index.                 		   */
  int kdim;			/* Dimension of geometry space.            */
  int kn;			/* Number of vertices of curve.            */
  int kk;			/* Order of curve.                         */
  int kpos = 0;			/* Current position in int.pt. array.      */
  double tpoint[3];		/* Value of point to intersect.            */
  double tpar0, tpar;		/* Parameter value of point on curve.      */
  double spar[1];		/* Parameter value of endpoint of curve.   */
  double sder[6];		/* Result of curve evaluation.             */
  double stang1[2];		/* Tangent vector of curve.                */
  double stang2[2];		/* Tangent vector of level value.          */
  double *st;			/* Pointer to knot vector of curve.        */
  double *sptpar = pintpt->epar;/* Pointer to parameter array of int.pt.   */
  double tref;			/* Referance value in equality test.       */
  SISLCurve *qc;		/* Pointer to current curve.               */
  SISLIntpt *uintpt[2];		/* Array storing new intersection points.  */
  double *ret_val;		/* Pointer to geo data from sh6getgeom     */
  double *ret_norm;		/* Pointer to geo data from sh6getgeom     */
  double *nullp = SISL_NULL;
  double dist;                  /* Distance from curve to point.           */
  double tot_ang;               /* Sum of angles between curve deriv. and 1*/
  int i;                        /* Loop variable.                          */
  int make_hp;                  /* Flag, make/not make help pt.            */
  /* --------------------------------------------------------------------- */

  /* Don't make help point for help points ! */
  if (sh6ishelp (pintpt))
    {
      *jstat = 0;
      goto out;
    }


  /* Test dimension of geometry space. */
  if (po1->iobj == SISLCURVE)
    {
      qc = po1->c1;
    }
  else
    {
      qc = po2->c1;
    }

  kdim = qc->idim;
  if (kdim != 2) goto err106;
	       
  /* Store curve information in local parameters. */

  kn = qc->in;
  kk = qc->ik;
  st = qc->et;
  tref = st[kn] - st[kk - 1];

  /* Fetch geometry information, point.  */
  sh6getgeom ((po1->iobj == SISLPOINT) ? po1 : po2,
	      (po1->iobj == SISLPOINT) ? 1 : 2,
	      pintpt, &ret_val, &ret_norm, aepsge, &kstat);
  if (kstat < 0)
    goto error;

  for(i=0; i<kdim; i++)
      tpoint[i] = ret_val[i];

  /* Fetch geometry information, curve.  */
  sh6getgeom ((po1->iobj == SISLCURVE) ? po1 : po2,
	      (po1->iobj == SISLCURVE) ? 1 : 2,
	      pintpt, &ret_val, &ret_norm, aepsge, &kstat);
  if (kstat < 0)
    goto error;

  s1219(st,kk,kn,&korgleft,sptpar[0],&kstat);
  if (kstat < 0) goto error;
  
  for(i=0; i<kdim*2; i++)
      sder[i] = ret_val[i];

/* Set tangent vectors. */

  stang1[0] = (double)1.0;
  stang2[0] = (double)1.0;
  stang2[1] = DZERO;
  tot_ang = (double)0.0;

  for(i=0; i<kdim; i++)
  {
      stang1[1] = sder[kdim+i];
      tot_ang += s6ang(stang1,stang2,2);
  }

  /*  tune */
  if (fabs(tot_ang) <= ANGULAR_TOLERANCE)
    {
      /* Test if the intersection point lies at the endpoint of
         the curve. */

      if (DEQUAL (sptpar[0] + tref, st[kn] + tref))
	{

	}
      else
	{
	  /* Find endpoint of coincidence interval in the positive
             direction of the curve. */
	  ki = 0;
	  tpar = sptpar[0] + (double) 2.0 *sqrt (aepsge);
	  tpar = min (tpar, st[kn]);
	  tpar0 = tpar = min (tpar, st[korgleft+1]);
	  shevalc (qc, 0, tpar, aepsge, &kleft, sder, &kstat);
          dist=s6dist(sder,tpoint,kdim);
          if (dist <= aepsge)
	    {
	       make_hp = TRUE;
	      for (ki = kleft - kk + 1; ki < kn; ki++)
		{
		  for (tpar = DZERO, kj = ki + 1; kj < ki + kk; kj++)
		    tpar += st[kj];
		  tpar /= (double) (kk - 1);

		  if (tpar > sptpar[0])
		    {
		      shevalc (qc, 0, tpar, aepsge, &kleft, sder, &kstat);
		      dist=s6dist(sder,tpoint,kdim);
		      if (dist > aepsge) break;
		      
		      tpar0 = tpar;  /* Store parameter value of 
					intersection point.          */
		    }
		}
	    }
	  else make_hp = FALSE;

	  /* Test if there is coincidence along the entire curve part. */

	  if (ki != kn && make_hp)
	    {
	      /* Create help point.  */
	       
	      spar[0] = tpar0;
	      uintpt[kpos] = SISL_NULL;
	      if ((uintpt[kpos] = hp_newIntpt (1, spar, DZERO, -SI_ORD,
					       SI_UNDEF, SI_UNDEF, SI_UNDEF,
				    SI_UNDEF, 0, 0, nullp, nullp)) == SISL_NULL)
		goto err101;

	      /* Insert the point into the data structure.  */

	      sh6idnpt (rintdat, &uintpt[kpos], 1, &kstat);
	      if (kstat < 0)
		goto error;

	      kpos++;
	    }
	}

      /* Test if the intersection point lies at the startpoint
         of the curve. */

      if (DEQUAL (sptpar[0] + tref, st[kk - 1] + tref))
	{
	}
      else
	{
	  /* Find endpoint of coincidence interval in the negative
             direction of the curve. */

	  ki = kn;
	  while (sptpar[0] == st[korgleft]) korgleft--;
	  tpar = sptpar[0] - (double) 2.0 *sqrt (aepsge);
	  tpar = max (tpar, st[kk - 1]);
	  tpar0 = tpar = max (tpar, st[korgleft]);	  
	  shevalc (qc, 0, tpar, aepsge, &kleft, sder, &kstat);
          dist=s6dist(sder,tpoint,kdim);
          if (dist <= aepsge)
	    {
	       make_hp = TRUE;
	       
	      for (ki = kleft; ki >= 0; ki--)
		{
		  for (tpar = DZERO, kj = ki + 1; kj < ki + kk; kj++)
		    tpar += st[kj];
		  tpar /= (double) (kk - 1);

		  if (tpar < sptpar[0])
		    {
		      shevalc (qc, 0, tpar, aepsge, &kleft, sder, &kstat);
		      dist=s6dist(sder,tpoint,kdim);
		      if (dist > aepsge) break;
		      
		      tpar0 = tpar;  /* Store parameter value of last
					found intersection point.  */
		    }
		}
	    }
	  else make_hp = FALSE;
	  
	  /* Test if there is coincidence along the entire curve part. */
	  if (ki >= 0 && make_hp)
	    {

	      /* Create intersection point.  */
	      spar[0] = tpar0;
	      uintpt[kpos] = SISL_NULL;
	      if ((uintpt[kpos] = hp_newIntpt (1, spar, DZERO, -SI_ORD,
					       SI_UNDEF,SI_UNDEF,SI_UNDEF,
					SI_UNDEF, 0, 0, nullp, nullp)) == SISL_NULL)
		goto err101;

	      /* Insert the point into the data structure.  */

	      sh6idnpt (rintdat, &uintpt[kpos], 1, &kstat);
	      if (kstat < 0)
		goto error;


	      kpos++;

	    }

	}
    }

  /* Join intersection points.  (kpos=0,1,2)*/
  for (ki = 0; ki < kpos; ki++)
    {
      sh6idnpt (rintdat, &uintpt[ki], 1, &kstat);
      if (kstat < 0)
	goto error;
      /* Mark that an intersection interval is found.  */
      if (sh6ishelp (uintpt[ki]) && uintpt[ki]->no_of_curves == 0)
	{
	  sh6idcon (rintdat, &uintpt[ki], &pintpt, &kstat);
	  if (kstat < 0)
	    goto error;
	}
    }

  /* Pre-topology information computed. */

  *jnewpt = kpos;
  *jstat = 0;
  goto out;

  /* Error in scratch allocation.  */

err101:*jstat = -101;
  goto out;

  /* Error in input. Incorrect dimension.  */

err106:*jstat = -106;
  goto out;

  /* Error lower level routine.  */

error:*jstat = kstat;
  goto out;

out:
  return;
}
