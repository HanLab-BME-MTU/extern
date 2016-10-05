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
 * $Id: s6findfac.c,v 1.4 2005-02-28 09:04:49 afr Exp $
 *
 */


#define S6FINDFAC

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s6findfac(double evecu[],double evecv[],double evecw[],double etang[],
               int idim,int isign,double *coef1,double *coef2,double *coef3,int *jstat)

#else
void s6findfac(evecu,evecv,evecw,etang,idim,isign,coef1,coef2,coef3,jstat)
     double evecu[];
     double evecv[];
     double evecw[];
     double etang[];
     int    idim;
     int    isign;
     double *coef1;
     double *coef2;
     double *coef3;
     int    *jstat;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : Given four vectors, evecu, evecv, evecw and etang, find 
*              the factors, coef1, coef2 and coef3, such that the vector
*              coef1*evecu + coef2*evecv + coef3*evecw = isign*etang.
*              
*
* INPUT      : evecu      - First vector.
*              evecv      - Second vector.
*              evecw      - Third vector.
*              etang      - Vector to approximate.
*              idim       - Dimension of geometry space.
*              isign      - Sign with wich etang is to be multiplied.
*
*
* OUTPUT     : coef1      - First factor.
*              coef2      - Second factor.
*              coef3      - Third factor.
*              jstat      - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : Minimize the square of the expression 
*                   dist(coef1*evecu+coef2*evecv,isign*etang)
*              over coef1 and coef2.
*              The expression is differentiated and set equal to
*              zero. Then this equation system of 2 equations
*              with two unknowns is solved. 
*              If the three vectors evecu, evecv and evecw span
*              3D, |coef3| = |d|/|evecw| where 
*              d = isign*etang - coef1*evecu - coef2*evecv
*
*********************************************************************
*/
{

  int kstat = 0;           /* Status variable.                    */
  int ki;                  /* Counter.                            */
  double tdotuu;           /* Scalar product of evecu and evecu.  */
  double tdotuv;           /* Scalar product of evecu and evecv.  */
  double tdotutang;        /* Scalar product of evecu and etang.  */
  double tdotvv;           /* Scalar product of evecv and evecv.  */
  double tdotvtang;        /* Scalar product of evecv and etang.  */
  double tdiv;             /* Determinant of equation system.     */
  double sdum[3];          /* Help vector.     */

  *jstat = 0;
  
  /* Test input.  */

  /* if (idim != 3) goto err104; */
  
  /* Set output to zero. */

  *coef1 = (double)0.0;
  *coef2 = (double)0.0;
  
  /* Compute coefficients of equation system.  */

  tdotuu = s6scpr(evecu,evecu,idim);
  tdotuv = s6scpr(evecu,evecv,idim);
  tdotutang = (double)isign*s6scpr(evecu,etang,idim);
  tdotvv = s6scpr(evecv,evecv,idim);
  tdotvtang = (double)isign*s6scpr(evecv,etang,idim);

  tdiv = tdotuv*tdotuv - tdotuu*tdotvv;
  if (DEQUAL(tdiv,DZERO))
    {
      if (DEQUAL(tdotuu,DZERO) && DEQUAL(tdotvv,DZERO));
      else if (DEQUAL(tdotuu,DZERO))
	  *coef2 = s6length(etang,idim,&kstat)/sqrt(tdotvv);
      else
	*coef1 = s6length(etang,idim,&kstat)/sqrt(tdotuu);
      goto out;
    }
  
  /* Compute the first two output factors.  */

  *coef1 = (tdotvtang*tdotuv - tdotutang*tdotvv)/tdiv;
  *coef2 = (tdotutang*tdotuv - tdotvtang*tdotuu)/tdiv;

  /* Find third output factor.  */

  for (ki=0; ki<idim; ki++) 
    sdum[ki] = (double)isign*etang[ki] - *coef1*evecu[ki] - *coef2*evecv[ki];
  *coef3 = s6length(sdum,idim,&kstat)/s6length(evecw,idim,&kstat);
  
  if (s6scpr(sdum,evecw,idim) < DZERO) (*coef3) *= -(double)1.0;

  goto out;


  out :
    return;
}
