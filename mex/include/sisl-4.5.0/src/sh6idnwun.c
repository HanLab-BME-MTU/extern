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
 * $Id: sh6idnwun.c,v 1.2 2001-03-19 15:59:08 afr Exp $
 *
 */


#define SH6IDNEWUNITE

#include "sislP.h"


#if defined(SISLNEEDPROTOTYPES)
void
sh6idnewunite (SISLObject *po1, SISLObject *po2, SISLIntdat ** intdat, 
	       SISLIntpt ** pt1, SISLIntpt ** pt2, double weight, 
	       double aepsge, int *jstat)
#else
void 
sh6idnewunite (po1, po2, intdat, pt1, pt2, weight, aepsge, jstat)
     SISLObject *po1;
     SISLObject *po2;
     SISLIntdat **intdat;
     SISLIntpt **pt1;
     SISLIntpt **pt2;
     double weight;
     double aepsge;
     int *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Given two Intpt, unite them to one by finding the
*              weighted medium in the parameter values corresponding
*              to one object, and iterate to the closest point in
*              the other. If one of the objects is a point, no
*              iteration is necessary.
*
*
* INPUT      : po1      - First object in the intersection.
*              po2      - Second object in the intersection.
*              intdat	- Pointer to intersection data.
*	       pt1	- Pointer to the first Intpt.
*	       pt2	- Pointer to the second Intpt.
*	       weight	- Weight to middle parameter values. ((1-W)*p1+W*p2).
*              aepsge   - Geometry tolerance. 
*
*
* OUTPUT     : pt1   	- Pointer to joined Intpt.
*	       jstat   	- Status value
*
*
* METHOD     :
*
* CALLS      : sh6ptobj - Iterate to closest point in second object.
*              s1221    - Curve evaluation.
*              s1421    - Surface evaluation.
*
*
* REFERENCES :
*
* WRITTEN BY : Arne Laksaa, SI, Oslo, Norway. May 91.
* CORRECTED BY: UJK, SI, Oslo, Norway. October 91.
* RENAMED AND MODIFIED BY : Vibeke Skytt, SI, 92-11. Iterate in second
*                                                    object.
*********************************************************************
*/
{
   int ki, kstat;
   int kpar;           /* Number of parameter directions in 1. object. */
   int kiterate;       /* Indicates if iteration is necessary.    */
   int kleft1=0,kleft2=0; /* Parameters used in evaluation.       */
   double spar[4];     /* Parameter values of intersection point. */
   double start[2];    /* Start parameter value to iteration.     */
   double spoint[3];   /* Position in curve or surface.           */
   double snorm[3];    /* Dummy vector. Surface normal.           */
   SISLIntpt *lpt;
   SISLIntpt *lpt1;
   SISLIntpt *lpt2;
   
   /* Test if one object is a point.  */
   
   if (po1->iobj == SISLPOINT || po2->iobj == SISLPOINT)
   {
      kpar = po1->iobj + po2->iobj;
      kiterate = 0;
   }
   else
   {
      kpar = po1->iobj;
      kiterate = 1;
   }

  sh6idnpt (intdat, pt1, 0, &kstat);
  if (kstat < 0)
    goto error;
  sh6idnpt (intdat, pt2, 0, &kstat);
  if (kstat < 0)
    goto error;

  if (sh6ismain (*pt1))
    {
      lpt1 = (*pt1);
      lpt2 = (*pt2);
    }
  else
    {
      lpt1 = (*pt2);
      lpt2 = (*pt1);
      weight = 1.0 - weight;
    }

  sh6disconnect (lpt1, lpt2, &kstat);
  if (kstat < 0)
    goto error;

  /* UJK, Oct. 91 */
  /* for (ki=0;;ki++) */
  for (ki = 0;;)
    {
      if ((lpt = sh6getnext (lpt2, ki)) == SISL_NULL)
	break;

      sh6disconnect (lpt2, lpt, &kstat);
      if (kstat < 0)
	goto error;


      sh6connect (lpt1, lpt, &kstat);
      if (kstat < 0)
	goto error;
    }

  for (ki = 0; ki < kpar; ki++)
    spar[ki] = lpt1->epar[ki] * (1.0 - weight) + lpt2->epar[ki] * weight;
  
  if (kiterate)
  {
     /* Compute start parameter to iteration.  */
     
     for (; ki < lpt1->ipar; ki++)
	start[ki-kpar] = lpt1->epar[ki] * (1.0 - weight) + lpt2->epar[ki] * weight;
	
     /* Iterate to closest point in second object. First evaluate
	value of intersection point in first object.  */
     
     if (po1->iobj == SISLCURVE)
     {
	s1221(po1->c1,0,spar[0],&kleft1,spoint,&kstat);
	if (kstat < 0) goto error;
     }
     else
     {
	s1421(po1->s1,0,spar,&kleft1,&kleft2,spoint,snorm,&kstat);
	if (kstat < 0) goto error;
     }
     
     /* Iterate. */
     
     sh6ptobj(spoint,po2,aepsge,start,spar+kpar,&kstat);
     if (kstat < 0) goto error;
  }
  
  /* Copy new parameter values into intersection point. */
  
  memcopy(lpt1->epar,spar,lpt1->ipar,DOUBLE);
     

  sh6idkpt (intdat, &lpt2, 0, &kstat);
  if (kstat < 0)
    goto error;

  (*pt1) = lpt1;
  (*pt2) = lpt2;

  goto out;

error:
  *jstat = kstat;
  s6err ("sh6idunite", kstat, 0);
  goto out;
out:
  ;
}

