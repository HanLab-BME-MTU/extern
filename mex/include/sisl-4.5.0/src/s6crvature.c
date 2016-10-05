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
 * $Id: s6crvature.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S6CURVATURE

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
      s6curvature(double eder[],int idim,double ecurv[],int *jstat)
#else	 
void s6curvature(eder,idim,ecurv,jstat)
     int idim,*jstat;
     double eder[],ecurv[];
#endif     
/*
*********************************************************************
*                                                                   
* PURPOSE    : Given position, first and second derivative of a curve at
*              a point, compute curvature vector.
*
*
*
* INPUT      : eder    - Array containing position, 1. and 2. derivative of
*                        curve. Dimension is 3*idim.
*              idim    - Dimension of geometry space.
*                       
*
* OUTPUT     : ecurv   - Curvature vector.
*              jstat   - status messages  
*                                         = 1      : Tangent vector zero.
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : Express the curve using cord length parametrisation. Then the
*              curvature is equal to the 2. derivative of the curve.
*
* REFERENCES : 
*
* USE        : 
*
*-
* CALLS      : s6length - Length of vector.   
*              s6scpr   - Scalar product between two vectors.  
*
* WRITTEN BY : Vibeke Skytt, SI, 04.90.
*
*********************************************************************
*/
{
  int kstat = 0;    /* Status variable.  */
  int ki;           /* Counter.   */
  double tleng;     /* Length of first derivative.  */
  double tleng2;    /* Square of length.  */
  double tdot;      /* Scalar product between 1. and 2. derivative. */
  
  /* Compute length of 1. derivative. */

  tleng = s6length(eder+idim,idim,&kstat);
  
  if (kstat == 0)
    {
      /* The first derivative is zero. */

      for (ki=0; ki<idim; ki++) ecurv[ki] = (double)0.0;
      goto warn1;
    }
  else
    {
      tleng2 = tleng*tleng;
      tdot = s6scpr(eder+idim,eder+2*idim,idim);
      
      for (ki=0; ki<idim; ki++)
	ecurv[ki] = (eder[2*idim+ki] - eder[idim*ki]*tdot/tleng2)/tleng2;
    }
  
  *jstat = 0;
  goto out;
  
  warn1 :
    *jstat = 1;
  goto out;
  
  out :
    return;
}

  
