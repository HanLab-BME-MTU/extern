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
 * $Id: s2511.c,v 1.3 2001-06-12 11:07:34 jbt Exp $
 *
 */


#define S2511

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
   s2511(SISLSurf *surf, int ider, double derive[], double normal[],
      double *mehlum, int *jstat)
#else
 void s2511(surf, ider, derive, normal, mehlum, jstat)
      SISLSurf *surf;
      int ider;
      double derive[];
      double normal[];
      double *mehlum;
      int *jstat;
#endif
/*
***************************************************************************
*
***************************************************************************
*  PURPOSE      :  To compute the third order Mehlum curvature M(u,v) of a 
*                  surface for given values (u,v). This is a lower level 
*                  routine, used for evaluation of many M(u,v)'s.
*  INPUT        :
*          surf     - Pointer to the surface to evaluate.
*          ider     - Only implemented for ider=0 (derivative order).
*       derive      - Array containing derivatives from routine s1421().
*                     Size = idim*6.
*       normal      - Array containing the normal from routine s1421().
*                     Size = 3.
*
*  OUTPUT       :
*       mehlum      - Third order Mehlum curvature of the surface.
*        jstat      - Status messages
*
*                         = 0 : Ok.
*                         < 0 : Error.
*
*  METHOD       :  The third order Mehlum curvature is given by
*
*                      M(u,v) = (5G^3P^2 
*                             + (EG + 4F^2)
*                             *(9GQ^2 + 9ES^2 + 6GPS + 6EQT) 
*                             + 5E^3T^2 
*                             - 2F(3EG + 2F^2)(PT + 9QS)
*                             - 30F(G^2PQ + E^2ST))
*                             /(16(EG - F^2)^3).
*
*                  The variables E,F,G,P,Q,S and T
*                  are the coefficients of the first and third fundamental form.
*                  They are given by: 
*                  E = <Xu,Xu>, F = <Xu,Xv> and G = <Xv,Xv>. 
*                  P = <N,Xuuu> + 3(a*alpha + b*beta), 
*                  Q = <N,Xuuv> + c*alpha + d*beta + 2a*gamma + 2b*delta, 
*                  S = <N,Xuvv> + 2c*gamma + 2d*delta + a*epsilon + b*mu,
*                  T = <N,Xvvv> + 3(c*epsilon + d*mu), 
*                  where N is normalized, and
*                  a = Ff - Ge,
*                  b = Fe - Ef,
*                  c = Fg - Gf,
*                  d = Ff - Eg,
*                  e, f and g being the second fundamental form coefficients
*                  (e = <N,Xuu>, f = <N,Xuv> and g = <N,Xvv>), and
*                  alpha   = <Xuu,Xu>/||N||^2,
*                  beta    = <Xuu,Xv>/||N||^2,
*                  gamma   = <Xuv,Xu>/||N||^2,
*                  delta   = <Xuv,Xv>/||N||^2,
*                  epsilon = <Xvv,Xu>/||N||^2,
*                  mu      = <Xvv,Xv>/||N||^2.
*
*  REFERENCES   :  Differential Geometry of Curves and Surfaces,
*                    (Manfredo P. Do Carmo, Prentice Hall,
*                      ISBN: 0-13-212589-7).
*-
*  CALLS        : s2513()
*
*  LIMITATIONS  :
*                (i) If the surface is degenerated (not regular) at the point
*                    (u,v), it makes no sense to speak about the Mehlum
*                    M(u,v).
*               (ii) If the surface is closed to degenerate, the Mehlum
*                    M(u,v) can be numerical unstable.
*              (iii) The dimension of the space in which the surface lies must
*                    be 1,2 or 3, if not, jstat = -105 is returned.
*
*
* WRITTEN BY   :  Johannes Kaasa, SINTEF, Oslo, Norway.            Date: 1995-9
*****************************************************************************
*/
{
   double fundform[10]; /* The coefficients of the fundamental forms.
			   The sequence is: E, F, G, e, f, g, P, Q, S, T. */
   double length;       /* Square of normal length.                       */
   double a, b, c, d;   /* Utility coefficients.                          */
   double alpha;        /* Utility coefficient.                           */
   double beta;         /* Utility coefficient.                           */
   double gamma;        /* Utility coefficient.                           */
   double delta;        /* Utility coefficient.                           */
   double epsilon;      /* Utility coefficient.                           */
   double mu;           /* Utility coefficient.                           */
   double P, Q, S, T;   /* Third order coefficients.                      */
   double numerator;    /* Value of a numerator.                          */
   double denominator;  /* Value of a denominator.                        */

   if (ider != 0) goto err178;

   if (surf->idim == 1 || surf->idim == 3) /* 1D and 3D surface */
   {
      s2513(surf, ider, 3, 1, derive, normal, fundform, jstat);
      if (*jstat < 0) goto error;
      
      if (surf->idim == 3)
	 length = normal[0]*normal[0] + normal[1]*normal[1] +
	    normal[2]*normal[2];
      else if (surf->idim == 1)
	 length = 1. + derive[1]*derive[1] + derive[2]*derive[2];
	 
      
      alpha   = s6scpr(&derive[3*(surf->idim)], &derive[surf->idim],
		     surf->idim)/length;
      beta    = s6scpr(&derive[3*(surf->idim)], &derive[2*(surf->idim)],
		     surf->idim)/length;
      gamma   = s6scpr(&derive[4*(surf->idim)], &derive[surf->idim],
		     surf->idim)/length;
      delta   = s6scpr(&derive[4*(surf->idim)], &derive[2*(surf->idim)],
		     surf->idim)/length;
      epsilon = s6scpr(&derive[5*(surf->idim)], &derive[surf->idim],
		     surf->idim)/length;
      mu      = s6scpr(&derive[5*(surf->idim)], &derive[2*(surf->idim)],
		     surf->idim)/length;
      
      a = fundform[1]*fundform[4] - fundform[2]*fundform[3];
      b = fundform[1]*fundform[3] - fundform[0]*fundform[4];
      c = fundform[1]*fundform[5] - fundform[2]*fundform[4];
      d = fundform[1]*fundform[4] - fundform[0]*fundform[5];
      
      P = fundform[6] + 3*(a*alpha + b*beta);
      Q = fundform[7] + c*alpha + d*beta + 2*a*gamma + 2*b*delta;
      S = fundform[8] + 2*c*gamma + 2*d*delta + a*epsilon + b*mu;
      T = fundform[9] + 3*(c*epsilon + d*mu);
      
      numerator = 5*fundform[2]*fundform[2]*fundform[2]*P*P
	 + (fundform[0]*fundform[2] + 4*fundform[1]*fundform[1])
	 *(9*fundform[2]*Q*Q
	 + 9*fundform[0]*S*S
	 + 6*fundform[2]*P*S
	 + 6*fundform[0]*Q*T)
	 + 5*fundform[0]*fundform[0]*fundform[0]*T*T
	 - 2*fundform[1]*(3*fundform[0]*fundform[2] + 2*fundform[1]*fundform[1])
	 *(P*T + 9*Q*S)
	 - 30*fundform[1]*(fundform[2]*fundform[2]*P*Q
	 + fundform[0]*fundform[0]*S*T);
      
      denominator = fundform[0]*fundform[2] - fundform[1]*fundform[1];
      
      *mehlum = numerator/(16*denominator*denominator*denominator);
  }

  else if (surf->idim == 2) /* 2D surface */
  {
    /* The surface lies in a plane => K(u,v) = 0 */

    *mehlum = 0.0;
  }
  else /* When surf->idim != 1,2 or 3 */
  {
    goto err105;
  }

  /* Successful computations  */

  *jstat = 0;
  goto out;


   /* Error in input, surf->idim != 1,2 or 3 */
err105:
  *jstat = -105;
  s6err("s2511",*jstat,0);
  goto out;

  /* Illegal derivative requested. */
err178:
  *jstat = -178;
  s6err("s2511",*jstat,0);
  goto out;
  
error:
  s6err("s2511",*jstat,0);
  goto out;

out:

  return;

}
