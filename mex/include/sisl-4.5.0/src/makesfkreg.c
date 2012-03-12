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
 * $Id: makesfkreg.c,v 1.6 1994-11-30 12:53:02 pfu Exp $
 *
 */


#define MAKE_SF_KREG

#include "sislP.h"


#if defined(SISLNEEDPROTOTYPES)
void make_sf_kreg (SISLSurf * ps, SISLSurf ** rsnew, int *jstat)
#else
void
   make_sf_kreg (ps, rsnew, jstat)
     SISLSurf *ps;
     SISLSurf **rsnew;
     int *jstat;
#endif
/*
********************************************************************
*
*********************************************************************
*
* PURPOSE    : To convert a surface to a k-regular basis.
*
*
*
* INPUT      : ps	- Surface to be made k-regular.
*
*
*
* OUTPUT     : rsnew	- The new surface on a k-regular basis.
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
*
**********************************************************************/
{
  int kn1=ps->in1;	/* Number of vertices in 1. par. dir.  */
  int kn2=ps->in2;	/* Number of vertices in 2. par. dir.  */
  int kk1=ps->ik1;	/* Order in 1. par. dir.               */
  int kk2=ps->ik2;	/* Order in 2. par. dir.               */
  /* --------------------------------------------------------- */

  s1001 (ps, ps->et1[kk1-1], ps->et2[kk2-1],
		ps->et1[kn1], ps->et2[kn2], rsnew, jstat);
  if (*jstat < 0)  goto error;

  goto out;

  /* Error in lower level routine */
error:
  s6err ("make_sf_kreg", *jstat, 0);

out:;

}
