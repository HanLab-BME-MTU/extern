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
 * $Id: crvlintang.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define CRV_LIN_TANG

#include "sislP.h"


#if defined(SISLNEEDPROTOTYPES)
void 
   crv_lin_tang(SISLCurve *pc1, double point[], double normal[],
	   double ang_tol, double guess_par, double *iter_par,
	   int *jstat)
#else
     void crv_lin_tang(pc1, point, normal, ang_tol, guess_par,
		       iter_par, jstat)
     SISLCurve   *pc1;
     double point[];
     double normal[];
     double ang_tol;
     double guess_par;
     double *iter_par;
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Newton iteration to find a tangent between the curve
*              pc1 and a line.
*
*
* INPUT      : pc1       - Pointer to the first curve.
*              point     - Original point on the line.
*              normal    - Normal to the line.
*              ang_tol   - Angular tolerance (in radians).
*              guess_par - Guess parameter values in pc1.
*
*
*
* OUTPUT     : iter_par  - Tangential parameter values in pc1.
*              jstat   - status messages  
*                                < 0   : error.
*
*
* METHOD     : Use po_crv_tang.c to find the tangent from point to pc1,
*              and check the tangent direction against normal.
*              Guess_par and iter_par must not be separated by a tangential
*              discontinuity.
*
*
* REFERENCES :
*
*
* WRITTEN BY : Johannes Kaasa, SI, March 1992.
*
*********************************************************************
*/                       
{                        
  int kstat = 0;            /* Local status variable.           */
  int kpos = 0;             /* Position of error.               */
  
  int kder = 0;             /* Evaluates only position.         */
  int kleft = 0;            /* Knot interval pointer.           */
  double iter_pnt[2];       /* Resulting point on the curve.    */ 
  double diffvec[2];        /* Vector between point and curve.  */
  double tangent[2];        /* Tangent along the line.          */
  int kdim = 2;             /* 2 dimensional.                   */
  double iter_ang;          /* Angular deviation.               */
  
  /* Test input.  */
  
  if (pc1->idim != 2) goto err106;
  
  /* Find the tangent from point to pc1. */
  
  po_crv_tang(pc1, point, ang_tol, guess_par, iter_par, &kstat);
  if (kstat < 0) goto error;
  
  /* Check the result. */

  s1221(pc1, kder, *iter_par, &kleft, iter_pnt, &kstat);
  if (kstat < 0) goto error;
  diffvec[0] = iter_pnt[0] - point[0];
  diffvec[1] = iter_pnt[1] - point[1];
  tangent[0] = -normal[1];
  tangent[1] = normal[0];
  iter_ang = s6ang(diffvec, tangent, kdim);
  if (iter_ang < ang_tol)
    *jstat = 1;
  else
    *jstat = 2;

  goto out;
  
  /* Error in input. Conflicting dimensions.  */
  
 err106: *jstat = -106;
  s6err("crv_lin_tang",*jstat,kpos);
  goto out;                  
  
  /* Error in lower level routine.  */
  
  error : *jstat = kstat;
  s6err("crv_lin_tang",*jstat,kpos);
  goto out;                  
  
 out: return;
}


