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
 * $Id: s2506.c,v 1.2 2001-03-19 15:58:59 afr Exp $
 *
 */


#define S2506

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
   s2506(SISLSurf *surf, int ider, int iside1, int iside2, double parvalue[],
      int *leftknot1,int *leftknot2, double *totalCurvature, int *jstat)
#else
 void s2506(surf, ider, iside1, iside2, parvalue, leftknot1, leftknot2,
	    totalCurvature, jstat)
      SISLSurf *surf;
      int    ider;
      int    iside1;
      int    iside2;
      double parvalue[];
      int *leftknot1;
      int *leftknot2;
      double *totalCurvature;
      int *jstat;
#endif
/*
***************************************************************************
*
***************************************************************************
*  PURPOSE      :  To compute the total curvature T(u,v) of a surface for
*                  given values (u,v) = (parvalue[0],parvalue[1]), where:
*
*                          et1[leftknot1] <= parvalue[0] < et1[leftknot1+1],
*                          et2[leftknot2] <= parvalue[1] < et2[leftknot2+1].
*
*  INPUT        :
*          surf     - Pointer to the surface to evaluate.
*          ider     - Number of derivatives to calculate.
*                     Only implemented for ider=0.
*                       < 0 : No derivative calculated.
*                       = 0 : Position calculated.
*                       = 1 : Position and first derivative calculated.
*                       etc.
*          iside1   - Indicator telling if the derivatives in the first
*                     parameter direction is to be calculated from the
*                     left or from the right:
*                        <  0 calculate derivative from the left hand side
*                        >= 0 calculate derivative from the right hand side.
*          iside2   - Indicator telling if the derivatives in the second
*                     parameter direction is to be calculated from the
*                     left or from the right:
*                        <  0 calculate derivative from the left hand side
*                        >= 0 calculate derivative from the right hand side.
*      parvalue     - Parameter-value at which to evaluate. Dimension of
*                     parvalue is 2.
*
*  INPUT/OUTPUT :
*     leftknot1     - Pointer to the interval in the knot vector in the
*                     first parameter direction where parvalue[0] is found,
*                     that is:
*                          et1[leftknot1] <= parvalue[0] < et1[leftknot1+1].
*                     leftknot1 should be set equal to zero at the first call
*                     to the routine.
*
*     leftknot1     - Pointer to the interval in the knot vector in the
*                     second parameter direction where parvalue[1] is found,
*                     that is:
*                          et2[leftknot2] <= parvalue[1] < et2[leftknot2+1].
*                     leftknot2 should be set equal to zero at the first call
*                     to the routine.
*
*  OUTPUT       :
*    totalCurvature - Total curvature of the surface in (u,v) =
*                     (parvalue[0],parvalue[1]).
*        jstat      - Status messages
*                         = 2 : Surface is degenerate at the point, that is,
*                               the surface is not regular at this point.
*                         = 1 : Surface is close to degenerate at the point.
*                               Angle between tangents is less than the angular
*                               tolerance.
*                         = 0 : Ok.
*                         < 0 : Error.
*
*
*  METHOD        :  The total curvature is given by
*
*                      T(x,y) = k1^2 + k2^2,
*
*                  if the surface (h(x,y)) is 1D, and
*
*                      T(u,v) = k1^2 + k2^2,
*
*                  if the surface (X(u,v)) is 3D. The variables E,F,G,e,f and g
*                  are the coefficients of the first and second fundamental form.
*                  They are given by: e = <N,Xuu>, f = <N,Xuv>, g = <N,Xvv>,
*                  E = <Xu,Xu>, F = <Xu,Xv> and G = <Xv,Xv>. The rutine will
*                  test if the surface is degenerate (not regular) or close to
*                  degenerate.
*
*  REFERENCES   :  Differential Geometry of Curves and Surfaces,
*                    (Manfredo P. Do Carmo, Prentice Hall,
*                      ISBN: 0-13-212589-7).
*-
*  CALLS        :  s1422() and s2507().
*
*  LIMITATIONS  :
*                (i) If the surface is degenerated (not regular) at the point
*                    (u,v), it makes now sense to speak about the absolute c.
*                    T(u,v). The routine return jstat == 2.
*               (ii) If the surface is closed to degenerate, the absolute c.
*                    T(u,v) can be numerical unstable. The routine return
*                    jstat == 1.
*              (iii) The surface should be C2, since the absolute c. is calculated
*                    from the second derivatives. But since the routine is using
*                    right derivatives, the absolute c. will be correct (provided
*                    that the surface is not degenerate).
*               (iv) The dimension of the space in which the surface lies must
*                    be 1,2 or 3.  The routine return jstat < 0.
*
*
* WRITTEN BY    :  Geir Westgaard,  SINTEF, Oslo, Norway.         Date: 1995-1
******************************************************************************
*/
{
  int kwarn = 0;         /* Local staus variable(warning).                  */
  int kistat = 0;        /* Local staus variable.                           */
  int der = 0;           /*  (dummy) */
  double derive[18];     /* Array containing the computed derivatives.      */
  double normal[3];      /* Array containing the computed normalvektor.     */


  if (ider != 0) goto err178;


  if (surf == SISL_NULL)  goto err150;
  else
  {
    /* Compute derivates and normal. */

    s1422(surf,2,iside1,iside2,parvalue,leftknot1,leftknot2,derive,normal,&kistat);
    if (kistat > 0) kwarn=kistat;

    if (kistat < 0) /* Error in lower level routine. */
    {
      goto error;
    }
    else if (kistat != 2) /* The surface is not degenerate */
    {
      s2507(surf, der, derive, normal, totalCurvature, &kistat);

      if (kistat < 0)
	goto error;
    }
    else if (kistat == 2) /* The surface is degenerated. */
    {
      *totalCurvature = 0.0;
      goto war002;
    }

  }



  /* Successful computations  */

  *jstat = kwarn;
  goto out;


   /* The surface is degenerated at (u,v) */
war002:
  *jstat = 2;
  goto out;

  /* Error. Input (surface) pointer is SISL_NULL. */
err150:
  *jstat = -150;
  s6err("s2506", *jstat, 0);
  goto out;

  /* Illegal derivative requested. */
err178:
  *jstat = -178;
  s6err("s2506",*jstat,0);
  goto out;
  /* Error in lower level routine.  */

error:
  *jstat = kistat;
  s6err("s2506",*jstat,0);
  goto out;


out:

  return;

}
