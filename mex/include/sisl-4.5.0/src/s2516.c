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
 * $Id: s2516.c,v 1.3 2001-06-12 11:07:34 jbt Exp $
 *
 */


#define S2516

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
   s2516(SISLSurf *surf, int ider, double derive[], double normal[],
      double mehlum[], int *stat)
#else
 void s2516(surf, ider, derive, normal, mehlum, stat)
      SISLSurf *surf;
      int ider;
      double derive[];
      double normal[];
      double mehlum[];
      int *stat;
#endif
/*
***************************************************************************
*
***************************************************************************
*  PURPOSE      :  To compute the numerator and denominator of the Mehlum 
*                  curvature M(u,v) of a Surface for given
*                  values (u,v). This is a lower level routine, used
*                  for evaluation of many M(u,v)'s.
*  INPUT        :
*          surf     - Pointer to the surface to evaluate.
*          ider     - Only implemented for ider=0 (derivative order).
*       derive      - Array containing derivatives from routine s1421().
*                     Size = idim*6.
*       normal      - Array containing the normal from routine s1421().
*                     Size = 3.
*
*  OUTPUT       :
*     mehlum        - The nominator and denominator of the Mehlum curvature 
*                     for the surface in (parvalue[0],parvalue[1]).
*                     Size = 2.
*        stat       - Status messages
*
*                         = 0 : Ok.
*                         < 0 : Error.
*
*  METHOD       :  The Mehlum curvature is given by
*
*                      M(u,v) = (3((eG-2fF+gE)^2)/8 - (eg-f*f)(EG-F*F)/2)/
*                               (EG-F*F)^3,
*
*                  The variables E,F,G,e,f and g
*                  are the coefficients of the first and second fundamental form.
*                  They are given by: e = <N,Xuu>, f = <N,Xuv>, g = <N,Xvv>,
*                  E = <Xu,Xu>, F = <Xu,Xv> and G = <Xv,Xv>. The routine will
*                  test if the surface is degenerate (not regular) or close to
*                  degenerate. 
*
*  REFERENCES   :  Differential Geometry of Curves and Surfaces,
*                    (Manfredo P. Do Carmo, Prentice Hall,
*                      ISBN: 0-13-212589-7).
*-
*  CALLS        :
*
*  LIMITATIONS  :
*                (i) If the surface is degenerated (not regular) at the point
*                    (u,v), it makes no sense to speak about the Mehlum
*                    curvature M(u,v).
*               (ii) If the surface is closed to degenerate, the Mehlum
*                    curvature M(u,v) can be numerical unstable.
*              (iii) If the surface is Cr the curvature calculated is C(r-2).
*               (iv) The dimension of the space in which the surface lies must
*                    be 1,2 or 3, if not, stat = -105 is returned.
*
*
* WRITTEN BY   :  Johannes Kaasa, SINTEF, Oslo, Norway.            Date: 1995-8
*****************************************************************************
*/
{
   double fundform[6];  /* The coefficients of the fundamental forms.
			   The sequence is: E, F, G, e, f, g.         */
   double coef1, coef2; /* Utility coefficients.                      */

   if (ider != 0) goto err178;

   if (surf->idim == 1 || surf->idim == 3) /* 1D and 3D surface */
   {
      s2513(surf, ider, 2, 0, derive, normal, fundform, stat);
      if (*stat < 0) goto error;
      
      coef1 = fundform[3]*fundform[2] - 2*fundform[4]*fundform[1]
		     + fundform[5]*fundform[0];
      coef2 = fundform[0]*fundform[2] - fundform[1]*fundform[1];
      
      mehlum[0] = 3.*coef1*coef1/8. - (fundform[3]*fundform[5] 
				       - fundform[4]*fundform[4])*coef2/2.;
      mehlum[1] = coef2*coef2*coef2;
  }

  else if (surf->idim == 2) /* 2D surface */
  {
    /* The surface lies in a plane => K(u,v) = 0 */

     mehlum[0] = 0.0;
     mehlum[1] = 1.0;
  }
  else /* When surf->idim != 1,2 or 3 */
  {
    goto err105;
  }

  /* Successful computations  */

  *stat = 0;
  goto out;


   /* Error in input, surf->idim != 1,2 or 3 */
err105:
  *stat = -105;
  s6err("s2516",*stat,0);
  goto out;

  /* Illegal derivative requested. */
err178:
  *stat = -178;
  s6err("s2516",*stat,0);
  goto out;
  
error:
  s6err("s2516",*stat,0);
  goto out;

out:

  return;

}
