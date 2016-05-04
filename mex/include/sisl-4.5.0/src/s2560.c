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

#define S2560

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s2560( SISLCurve *curve,
       double     parvalue,
       int       *leftknot,
       double     derive[],
       double     p[],
       double     t[],
       double     n[],
       double     b[],
       int       *jstat )
#else
void s2560( curve, parvalue, leftknot, derive, p, t, n, b, jstat )
     SISLCurve  *curve;
     double      parvalue;
     int        *leftknot;
     double      derive[];
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
*              at a given parameter value, from the right hand side.
*
*
*
* INPUT      : curve    - Pointer to the curve.
*              parvalue - The parameter value at which to compute
*                         the Frenet Frame.
*
*
*
* INPUT/OUTPUT : leftknot - Pointer to the interval in the knot vector
*                        where ax is located. If et is the knot vector,
*                        the relation
*
*                          et[ileft] < parvalue <= et[ileft+1]
*
*                        should hold. (If parvalue == et[ik-1] then ileft
*                        should be ik-1. Here in is the number of B-spline
*                        coefficients.)
*                        If ileft does not have the right value upon
*                        entry to the routine, its value will be changed
*                        to the value satisfying the above condition.
*
*
*
* OUTPUT     : derive   - Double array of dimension [3*idim]
*                       containing the position and derivative vectors.
*                       (idim is the number of components of each B-spline
*                       coefficient, i.e. the dimension of the Euclidean
*                       space in which the curve lies.)
*                       These vectors are stored in the following order:
*                       First the idim components of the position vector,
*                       then the idim components of the tangent vector,
*                       then the idim components of the second derivative
*                       vector, and so on.
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
*
* METHOD     : The derivatives are evaluated from the right hand
*              side by  s1221(), and the Frenet Frame are evaluated by s2561()
*
* REFERENCES :
*
*-
* CALLS      : s1221(), s2561()
*
* WRITTEN BY  :   Geir Westgaard, SINTEF, Oslo, November 1999
* MODIFIED BY :
* REVISED BY  :
*
*********************************************************************
*/
{
  int kdim = curve -> idim;   /* copy curve attribute to local parameter  */
  int kstat = 0;              /* local status variable                    */
  int kpos = 0;               /* local error position                     */





  /* Evaluate the derivatives */

  s1221( curve, 2, parvalue, leftknot, derive, &kstat );

  if ( kstat < 0 ) goto error;


  /* Evaluate the Frenet Frame based on the derivatives in "derive" */

  s2561( derive, kdim, p, t, n, b, &kstat );

  if ( kstat < 0 ) goto error;



 /* Successful computations.  */

  *jstat = 0;
  goto out;

  /* Error in lower level routine.  */

 error:  *jstat = kstat;
  s6err( "s2560", *jstat, kpos );
  goto out;


 out:
 return;

}
