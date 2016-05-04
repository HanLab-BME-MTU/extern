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
 * $Id: s1380.c,v 1.2 2001-03-19 15:58:48 afr Exp $
 *
 */


#define S1380

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1380(double ep[],double ev[],int im,int idim,int ipar,
	   SISLCurve **rcurve,int *jstat)
#else
void s1380(ep,ev,im,idim,ipar,rcurve,jstat)
     double ep[];
     double ev[];
     int    im;
     int    idim;
     int    ipar;
     SISLCurve  **rcurve;
     int    *jstat;
#endif
/*
************************************************************************
*
* Purpose: To compute the cubic Hermit interpolant to the data given
*          by the points ep and the derivatives ev. 
*          The curve is represented as a B-spline curve.
*
* Input:
*          ep     - Array containing the point in sequence
*                   (x,y,..,x,y,..)
*          ev     - Array containing the derivatives in sequence
*                   (x,y,..,x,y,..)
*          im     - Number of point and derivatives
*          idim   - The dimension of the space the points and derivatives
*                   lie in
*          ipar   - Type of parametrization
*                    1 - Parametrization using cordlength between point
*                  !=1 - Uniform parametrization
* Output:
*          rcurve - Pointer to the curve produced
*          jstat  - Status variable
*                    < 0 - Error.
* Method:
*     First the parmaterization is calculated and then the interpolation
*     is performed using s1379.
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Tor Dokken, SI 1988-11
*
*********************************************************************
*/                                                               
{
  int ki;             /* Loop variables                              */
  int kpek1,kpek2;    /* Pointers into point array                   */
  int kstat;          /* Status variable                             */
  int kpos=0;         /* Position of error                           */
  double *spar=SISL_NULL;  /* Pointer to parametrization array            */
  
  
  
  /* Check input */        
  
  if (im < 2)   goto err181;
  if (idim < 1) goto err102;
  
  /* Allocate array for parametrization */
  
  spar = newarray(im,DOUBLE);
  if (spar == SISL_NULL) goto err101;
  
  spar[0] = (double)0.0;
  
  if (ipar == 1)                 
    {
      /*  Cord length parametrization */
      
      kpek1 = 0;
      for (ki=1 ; ki<im ; ki++)
	{
	  kpek2 = kpek1 + idim;
	  spar[ki] = spar[ki-1] + s6dist(&ep[kpek2],&ep[kpek1],idim);
	  kpek1 = kpek2;
	}
    }
  else
    {
      /*  Uniform parametrization */
      for (ki=0;ki<im;ki++)
	spar[ki] = ki;
    }
  
  /* Calculate Hermite interpolant */
  
  s1379(ep,ev,spar,im,idim,rcurve,&kstat);
  if (kstat<0) goto error;
  
  /* Calculation completed */
  
  *jstat = 0;
  goto out;
  
  
  /* Error in space allocation. Return zero. */
  
  
  /* Error in space allocation */
 err101: *jstat = -101;
  s6err("s1380",*jstat,kpos);
  goto out;
  
  
  /* Dimension less than 1*/
 err102: *jstat = -102;
  s6err("s1380",*jstat,kpos);
  goto out;
  
  /* Too few interpolation conditions */
  
 err181: *jstat = -181;
  s6err("s1380",*jstat,kpos);
  goto out;
  
  /* Error in lower level routine */
  
 error:  *jstat = kstat;
  s6err("s1380",*jstat,kpos);
  goto out;
  
  
 out:
  if (spar != SISL_NULL) freearray(spar);
  
  return;
}
