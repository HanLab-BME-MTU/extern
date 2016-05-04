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
 *
 *
 */


#define S1517

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void s1517(double ep[],double ev[],double epar[],int im, double mu,
           double **evnew,int *jstat)
#else
void s1517(ep,ev,epar,im,mu,evnew,jstat)
     double ep[];
     double ev[];
     double epar[];
     int    im;
     double mu;
     double **evnew;
     int    *jstat;
#endif
/*
************************************************************************
*
* Purpose:   Given positive values and arbitrary derivatives, adjust
*            the derivatives so that the cubic Hermite interpolant will
*            be positive (in fact so that the coefficients of each
*            cubic Bezier segment will be positive).
*
* Input:
*          ep     - Array containing the point in sequence
*                   (x,y,..,x,y,..), length im.
*          ev     - Pointer to array containing the derivatives in sequence
*                   (x,y,..,x,y,..), length im.
*          epar   - Parametrization array. The array should be increasing
*                   in value, length im.
*          im     - Number of point and derivatives
*          mu     - Factor used for adjustment (0 <= mu < 1):
*                    mu = 0 => all derivatives will be zero,
*                    mu = 1 => all Bezier coeffs will be >=0, some
*                               can be exactly 0.
*
* Output:
*          evnew  - Pointer to array containing the adjusted derivatives
*                   in sequence (x,y,..,x,y,..), length im.
*          jstat  - Status variable
*                    < 0 - Error.
* Method:
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Micahel Floater, SI 1993-10
*
*********************************************************************
*/
{
  int ki;             /* Loop variables                              */
  int kpos=0;         /* Position of error                           */
  double *evtemp;
  double mu3;
  double knotint1;
  double knotint2;




  /* Check input */

  if (im < 2) goto err102;
  if(mu < 0.0 || mu >= 1.0) goto err102;

  mu3 = 3.0 * mu;

  /* Allocate array for derivatives */

  evtemp    = newarray(im,DOUBLE);
  if (evtemp == SISL_NULL) goto err101;

  /* Adjust the derivatives. */


  knotint2 = epar[1] - epar[0];

  if(ev[0] < - mu3 * ep[0] / knotint2)
  {
      evtemp[0] = - mu3 * ep[0] / knotint2;
  }
  /* @JKA: This is copied from the uncommented version. */
  else if (mu == DZERO)
  {
     evtemp[0] = mu;
  }
  /* @JKA */
  else
  {
      evtemp[0] = ev[0];
  }


  for(ki=1; ki<im-1; ki++)
  {
      knotint1 = epar[ki] - epar[ki-1];
      knotint2 = epar[ki+1] - epar[ki];

      if(ev[ki] >  mu3 * ep[ki] / knotint1)
      {
          evtemp[ki] = mu3 * ep[ki] / knotint1;
      }
      else if(ev[ki] < - mu3 * ep[ki] / knotint2)
      {
          evtemp[ki] = - mu3 * ep[ki] / knotint2;
      }
      else
      {
          evtemp[ki] = ev[ki];
      }
  }

  knotint1 = epar[im-1] - epar[im-2];

  if(ev[im-1] > mu3 * ep[im-1] / knotint1)
  {
      evtemp[im-1] = mu3 * ep[im-1] / knotint1;
  }
  /* @JKA: This is copied from the uncommented version. */
  else if (mu == DZERO)
  {
     evtemp[im-1] = mu;
  }
  /* @JKA */
  else
  {
      evtemp[im-1] = ev[im-1];
  }


  /* Calculation completed */

  /* Set result. */

  (*evnew) = evtemp;

  *jstat = 0;
  goto out;



  /* Error in space allocation */

 err101: *jstat = -101;
  s6err("s1517",*jstat,kpos);
  goto out;


  /* Error in input. */

 err102: *jstat = -102;
  s6err("s1517",*jstat,kpos);
  goto out;


 out:

  return;
}
