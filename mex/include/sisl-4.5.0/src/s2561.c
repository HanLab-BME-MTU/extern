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
 * $Id:
 *
 */


#define S2561

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s2561( double     derive[],
       int        idim,
       double     p[],
       double     t[],
       double     n[],
       double     b[],
       int       *jstat )
#else
void s2561( derive, idim, p, t, n, b, jstat )
     double      derive[];
     int         idim;
     double      p[];
     double      t[];
     double      n[];
     double      b[];
     int        *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Evaluate the Frenet Frame (t,n,b) of a curve
*              based on the derivatives in "derive". We assume that
*              the curve lies in R, R^2 or R^3.
*
*
* INPUT      :
*            derive   - Double array of dimension [3*idim]
*                       containing the position and derivative vectors.
*                       These vectors are stored in the following order:
*                       First the position vector, then the tangent vector,
*                       then the idim components of the second derivative
*                       vector, and so on.
*
*
* OUTPUT     :
*
*         (t,n,b) - The Frenet Frame (in 3D) computed. Each of the vectors
*                   (t,n,b) are of dim. 3, and the data are
*                   stored like this: tx(parvalue), ty(parvalue), tz(parvalue).
*
*         p     - 3D curve posistions at parvalue.
*
*              jstat  - Status messages
*                                         > 0      : Warning.
*                                         = 0      : Ok.
*                                         < 0      : Error.
*
*  METHOD     : See formula in the book of Do Carmo
*              (Differential Geometry of Curves and Surfaces).
*
* REFERENCES : Differential Geometry of Curves and Surfaces,
*                    (Manfredo P. Do Carmo, Prentice Hall,
*                      ISBN: 0-13-212589-7).
*
*
*
* CALLS       : s6length(), s6crss()
*
* WRITTEN BY  :   Geir Westgaard, SINTEF, Oslo, November 1999
* MODIFIED BY :
* REVISED BY  :
*
*********************************************************************
*/
{
  int kstat = 0;        /* Local status variable      */
  double a;             /* Tmp. variable.             */
  double crpr[3];       /* Tmp. variable.             */
  double D[9];          /* Derivatives on canonical
			   R^3 form                   */


  /* Put the derivatives on canonical R^3 form */

  if ( idim == 1 )
  {
      D[0]  = 0.0;   D[1] = derive[0];  D[2]  = 0.0;
      D[3]  = 1.0;   D[4] = derive[1];  D[5]  = 0.0;
      D[6]  = 0.0;   D[7] = derive[2];  D[8]  = 0.0;
  }
  else if ( idim == 2 )
  {
      D[0]  = derive[0];   D[1] = derive[1];   D[2]  = 0.0;
      D[3]  = derive[2];   D[4] = derive[3];   D[5]  = 0.0;
      D[6]  = derive[4];   D[7] = derive[5];   D[8]  = 0.0;
  }
  else
  {
      D[0]  = derive[0];   D[1] = derive[1];   D[2]  = derive[2];
      D[3]  = derive[3];   D[4] = derive[4];   D[5]  = derive[5];
      D[6]  = derive[6];   D[7] = derive[7];   D[8]  = derive[8];
  }

  /* Set 3D curve pos. */

  p[0] = D[0];
  p[1] = D[1];
  p[2] = D[2];


  /* Compute the Frenet Frame (t,n,b). */

  a = s6length( D + 3, 3, &kstat );

  t[0] = D[3]/a;
  t[1] = D[4]/a;
  t[2] = D[5]/a;


  s6crss( D + 3, D + 6, crpr );

  a = s6length( crpr, 3, &kstat );

  if ( kstat != 0 )
  {
    b[0] = crpr[0]/a;
    b[1] = crpr[1]/a;
    b[2] = crpr[2]/a;

    s6crss( b, t, n );
  }
  else
  {
    n[0] = 0.0;
    n[1] = 0.0;
    n[2] = 0.0;

    b[0] = 0.0;
    b[1] = 0.0;
    b[2] = 0.0;

    goto war002;
  }



 /* Successful computations.  */

  *jstat = 0;
  goto out;


  /* Frenet Frame undefined, since the curvature = 0.0.  */
war002:
  *jstat = 2;
  goto out;


 out:

 return;

}
