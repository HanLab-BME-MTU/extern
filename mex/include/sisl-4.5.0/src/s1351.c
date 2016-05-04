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
 * $Id: s1351.c,v 1.2 2001-03-19 15:58:46 afr Exp $
 *
 */

#define S1351

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void s1351(double ep[],int ipar,
	   int im,int idim,int ik,
	   SISLCurve **rc, int *jstat)
#else
void s1351(ep,ipar,im,idim,ik,rc,jstat)
     double ep[];
     int ipar;
     int im;
     int idim;
     int ik;	   
     SISLCurve **rc;
     int    *jstat;
#endif
/*
************************************************************
*
* Purpose: To compute the piecewise linear interpolant to a set 
*          of datapoints and express it as a linear combination 
*          of B-splines of order ik using the parametrization
*          determined by ipar.
*
* Input :
*        Ep     - Array [idim,im] containing the points to
*                 be approximated.
*        Ipar   - Integer indicating type of parametrization:
*                 1: Chord length parametrization.
*                 2: Uniform parametrization.
*        Im     - The no. of data points.
*        Idim   - The dimension of the euclidean space in which the data
*                 points lie, i.e. the number of components of each data point.
*        Ipar   - Flag indicating the type of parameterization to be used:
*                  = 1 : Paramterize by accumulated cord length.
*                        (Arc length parametrization for the piecewise
*                        linear interpolant.)
*                  = 2 : Uniform parameterization.
*                  = 3 : Parametrization given by epar.
*                 If ipar<1 or ipar>3, it will be set to 1.
*        Ik     - The polynomial order of the approximation.
*
* Output:
*        Jstat  - Output staus:
*                  < 0 : Error.
*                  = 0 : Ok.
*                  > o : Warning.
*        Rc     - Pointer to curve.
*
* Method: 
*
* The fortran version was written by Knut M|rken,  Si.
* Written by: C.R.Birkeland  Si  Oslo,Norway April 1993.
*
********************************************************************
*/
{
  int stat = 0;                  /* Error-Status parameters       */ 
  int kpos = 0;
  int par=0;                     /* Type of parametrization       */
  double *epar = SISL_NULL;           /* Parametrization array         */ 
  int i;                         /* Loop index                    */ 
  int kpek1, kpek2;              /* Used in gen. of cord len. par.*/

  /* Check Input */
  
  if (im < 2 || idim < 1 || ik < 2) goto err103;
  if (ipar < 1 || ipar > 2) goto err103;

  /* Allocate array for parametrization */

  if( (epar = newarray(im, DOUBLE)) == SISL_NULL ) goto err101;
  epar[0] = (double) 0.0;

  /* Compute parametrization */

  if (ipar == 1)
    {
      /* Chord length parametrization */

      kpek1 = 0;
      for (i=1; i<im; i++)
	{
	  kpek2 = kpek1 + idim;
	  epar[i] = epar[i-1] + s6dist(&ep[kpek2], &ep[kpek1], idim);
	  kpek1 = kpek2;
	}
      if (epar[im-1] == 0.) par=2;
    }

  if (ipar == 2 || par == 2)
    /* Uniform parametrization */

    for (i=1; i<im; i++)
      epar[i] = i;
      
  /* Compute cubic spline hermite interpolant */

  s1350(ep, epar, im, idim, ik, rc, &stat);
  if (stat<0) goto error;
  
  /* Success */

  *jstat = 0;
  goto out;


  /* Error in scratch allocation.  */

  err101 :
    *jstat = -101;
    s6err("s1351",*jstat,kpos);
    goto out;

  /* Error in input */

  err103: 
    *jstat = -103;
    s6err("s1351",*jstat,kpos);
    goto out;
  
  /* Error in lower level routine. */

  error:
    *jstat = stat;
    s6err("s1351",*jstat,kpos);
    goto out;

  /* Exit */

  out:
    if(epar != SISL_NULL) freearray(epar);
    return;
}
