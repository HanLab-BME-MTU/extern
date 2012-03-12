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
 * $Id: s1530.c,v 1.2 2001-03-19 15:58:50 afr Exp $
 *
 */


#define S1530

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void s1530(double ep[],double eder10[],double eder01[],double eder11[],
	    double epar1[],double epar2[],
	    int im1,int im2,int idim,
	   SISLSurf **rsurf,int *jstat)
#else
void s1530(ep,eder10,eder01,eder11,epar1,epar2,
	    im1,im2,idim,rsurf,jstat)
     double ep[];
     double eder10[];
     double eder01[];
     double eder11[];
     double epar1[];
     double epar2[];
     int    im1;
     int    im2;
     int    idim;
     SISLSurf  **rsurf;
     int    *jstat;
#endif
/*
************************************************************************
*
* PURPOSE: To compute the cubic Hermite interpolant to the data given.
*          More specifically, given positions, 10, 01, and 11
*          derivatives at points of a rectangular grid, the routine
*          computes a cubic tensor-product B-spline interpolant to
*          the given data with double knots at each data (the first
*          knot vector will have double knots at all interior points
*          in epar1, quadruple knots at the first and last points,
*          and similarly for the second knot vector).
*
* INPUT:
*          ep     - Array of dimension idim*im1*im2 containing
*                   the positions of the nodes (using the same ordering
*                   as ecoef in the SISLSurf structure).
*
*          eder10 - Array of dimension idim*im1*im2 containing the
*                   first derivative in the first parameter direction.
*
*          eder01 - Array of dimension idim*im1*im2 containing the
*                   first derivative in the second parameter direction.
*
*          eder11 - Array of dimension idim*im1*im2 containing
*                   the cross derivative (twist vector).
*
*          epar1  - Array of dimension im1 containing the
*                   parametrization in the first direction.
*
*          epar2  - Array of dimension im2 containing the
*                   parametrization in the first direction.
*
*          im1    - The number of interpolation points in the
*                   first parameter direction.
*
*          im2    - The number of interpolation points in the
*                   second parameter direction.
*
*          idim   - Dimension of the space we are working in.
*
* Output:
*          rsurf - Pointer to the surf produced
*          jstat  - Status variable
*                    < 0 - Error.
*
* Method:
*     The interpolation is accomplished by using a one dimensional
*     routine for cubic Hermite spline interpolation. First, the data
*     is considered to be im2 positional and derivative vectors on
*     two curves in idim * im1 dimensional space sampled at the
*     points of epar2.
*     The first of these has position vectors given by ep and
*     derivative vectors given by eder01, the second position vectors
*     given by eder10 and derivative vectors given by eder11.
*     Running these curves through the one dimensional cubic Hermite
*     spline interpolation routine then produces two cubic splines rpos
*     and rder with coefficients of dimension (idim * im1) * (2 * im2)
*     on the knot vector et2 which is just the points of epar2 with
*     multiplicity 2 for the interior points and 4 for the
*     end points.
*     These coefficients are then considered to be im1 position vectors
*     and derivative vectors on a curve in 2*idim*im2 dimensional
*     space (after an appropriate tranposition) sampled at the
*     points of epar1. Running this data through the one dimensional
*     cubic Hermite spline routine results in a cubic spline
*     with coefficients of dimension (2 * idim * im2) * (2 * im1)
*     with knot vector et1 similar to et2.
*     A transposition of these coefficients yields the B-spline
*     coefficients of the bicubic Hermite tensor-product spline
*     interpolant.

* REFERENCES :
*
* CALLS      : 
*
* WRITTEN BY : Michael Floater, SI, June 1992.
*
*********************************************************************
*/                                                               
{
  SISLCurve *rpos;    /* Curve of positions of dimension idim*im1    */
  SISLCurve *rder;    /* Curve of derivatives of dimension idim*im1  */
  SISLCurve *rcurve;  /* Curve of dimension idim*im2                 */
  int kstat=0;        /* Status variable                             */
  int kpos=0;         /* Position of error                           */
  double *ph=SISL_NULL;    /* Transposed positions (in rpos)              */
  double *dh=SISL_NULL;    /* Transposed derivatives (in rder)            */
  double *scoef=SISL_NULL; /* Transposed positions in rcurve              */
  
  
  
  /* Check input */        
  
  if (im1 < 2 || im2 < 2 || idim < 1)   goto err102;
  

  /* Interpolate position and derivative
     in the second parameter direction. */

  s1379(ep,eder01,epar2,im2,idim * im1,&rpos,&kstat);
  if(kstat < 0) goto error;

  s1379(eder10,eder11,epar2,im2,idim * im1,&rder,&kstat);
  if(kstat < 0) goto error;



  /* Transpose the matrices (ecoef arrays) of the curves. */

  s1531(rpos->ecoef, idim, im1, rpos->in, &ph, &kstat);
  if(kstat < 0) goto error;

  s1531(rder->ecoef, idim, im1, rder->in, &dh, &kstat);
  if(kstat < 0) goto error;



  /* Interpolate in the first parameter direction. */

  s1379(ph,dh,epar1,im1,idim * rpos->in,&rcurve,&kstat);
  if(kstat < 0) goto error;

  /* Transpose coefficient matrix of rcurve. */

  s1531(rcurve->ecoef, idim, rpos->in, rcurve->in, &scoef, &kstat);
  if(kstat < 0) goto error;

  /* Create instance of surface. */

  (*rsurf) = newSurf(rcurve->in,rpos->in,rcurve->ik,rpos->ik,
			rcurve->et,rpos->et,scoef,1,idim,1);
  if((*rsurf) == SISL_NULL) goto err101;
  
  /* Set periodicity flag. */
  
  (*rsurf)->cuopen_1 = rcurve->cuopen;
  (*rsurf)->cuopen_2 = rpos->cuopen;

  /* Bicubic surface made. */
  
  *jstat = 0;
  goto out;
  
  /* Error in space allocation. */

 err101: *jstat = -101;
  s6err("s1530",*jstat,kpos);
  goto out;
  
  
  /* Error in input data. */

 err102: *jstat = -102;
  s6err("s1530",*jstat,kpos);
  goto out;
  
  
  /* Error in lower level routine. */

 error:  *jstat =kstat;
  s6err("s1530",*jstat,kpos);
  goto out;
  
 out:
  if (rpos != SISL_NULL) freeCurve(rpos);
  if (rder != SISL_NULL) freeCurve(rder);
  if (rcurve != SISL_NULL) freeCurve(rcurve);
  if (scoef != SISL_NULL) freearray(scoef);
  if (ph != SISL_NULL) freearray(ph);
  if (dh != SISL_NULL) freearray(dh);
  
  return;
}
