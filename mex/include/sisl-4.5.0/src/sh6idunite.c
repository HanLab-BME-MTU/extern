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
 * $Id: sh6idunite.c,v 1.2 2001-03-19 16:06:03 afr Exp $
 *
 */


#define SH6IDUNITE

#include "sislP.h"


#if defined(SISLNEEDPROTOTYPES)
void
sh6idunite (SISLIntdat ** intdat, SISLIntpt ** pt1, SISLIntpt ** pt2,
	    double weight, int *jstat)
#else
void 
sh6idunite (intdat, pt1, pt2, weight, jstat)
     SISLIntdat **intdat;
     SISLIntpt **pt1;
     SISLIntpt **pt2;
     double weight;
     int *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Given two Intpt, unite them to one.
*
*
* INPUT      : intdat	- Pointer to intersection data.
*	       pt1	- Pointer to the first Intpt.
*	       pt2	- Pointer to the second Intpt.
*	       weight	- Weight to middle parameter values. ((1-W)*p1+W*p2).
*
*
* OUTPUT     : pt1   	- Pointer to joined Intpt.
*	       jstat   	- Status value
*
*
* METHOD     :
*
*
* REFERENCES :
*
* WRITTEN BY : Arne Laksaa, SI, Oslo, Norway. May 91.
* CORRECTED BY: UJK, SI, Oslo, Norway. October 91.
*********************************************************************
*/
{
  int ki, kstat;
  SISLIntpt *lpt;
  SISLIntpt *lpt1;
  SISLIntpt *lpt2;

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

  for (ki = 0; ki < lpt1->ipar; ki++)
    lpt1->epar[ki] = lpt1->epar[ki] * (1.0 - weight) + lpt2->epar[ki] * weight;

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

