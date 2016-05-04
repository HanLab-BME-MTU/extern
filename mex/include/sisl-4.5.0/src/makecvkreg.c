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
 * $Id: makecvkreg.c,v 1.7 1994-11-30 14:37:16 pfu Exp $
 *
 */


#define MAKE_CV_KREG

#include "sislP.h"


#if defined(SISLNEEDPROTOTYPES)
void
    make_cv_kreg (SISLCurve * pc, SISLCurve ** rcnew, int *jstat)
#else
void
   make_cv_kreg (pc, rcnew, jstat)
     SISLCurve *pc;
     SISLCurve **rcnew;
     int *jstat;
#endif
/*
********************************************************************
*
*********************************************************************
*
* PURPOSE    : To convert a curve to a k-regular basis.
*
*
*
* INPUT      : pc	- Curve to be made k-regular.
*
*
*
* OUTPUT     : rcnew	- The new curve on a k-regular basis.
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
*
* WRITTEN BY : Ulf J. Krystad, SI, 04.92.
* Revised by : Paal Fugelli, SINTEF, Oslo, Norway, 94-08. Added error propagation.
* Revised by : Paal Fugelli, SINTEF, Oslo, Norway, Nov.1994. Set cuopen flag to
*              closed when changed from periodic.
**********************************************************************/
{
   int kn=pc->in;	/* Number of vertices in 1. par. dir.  */
   int kk=pc->ik;	/* Order in 1. par. dir.               */
   /* --------------------------------------------------------- */
   /* Pick part of curve */
   s1712 (pc, pc->et[kk-1], pc->et[kn], rcnew, jstat);
  if (*jstat < 0)  goto error;

   if (pc->cuopen == SISL_CRV_PERIODIC )
     (*rcnew)->cuopen = SISL_CRV_CLOSED;

  goto out;

  /* Error in lower level routine */
error:
  s6err ("make_cv_kreg", *jstat, 0);

out:;

}
