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


#define S2555

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s2555( double     derive[],
       double    *torsion,
       int       *jstat )
#else
void s2555( derive, torsion, jstat )
     double      derive[];
     double     *torsion;
     int        *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Evaluate the torsion of a curve based on the
*              derivatives in "derive". We assume that the
*              curve lies in IR^3.
*
*
* INPUT      :
*            derive   - Double array of dimension [4*3]
*                       containing the position and derivative vectors.
*                       These vectors are stored in the following order:
*                       First the position vector, then the tangent vector,
*                       then the 3 components of the second derivative
*                       vector, and so on.
*
*
* OUTPUT     : torsion - The torsion value computed
*              jstat  - Status messages
*                                         > 0      : Warning.
*                                         = 0      : Ok.
*                                         < 0      : Error.
*
*
* METHOD     : See formula in the book of Do Carmo
*              (Differential Geometry of Curves and Surfaces).
*
* REFERENCES : Differential Geometry of Curves and Surfaces,
*                    (Manfredo P. Do Carmo, Prentice Hall,
*                      ISBN: 0-13-212589-7).
*
*
* CALLS       : s6length(), s6scpr()
*
* WRITTEN BY  :   Geir Westgaard, SINTEF, Oslo, November 1999
* MODIFIED BY :
* REVISED BY  :
*
*********************************************************************
*/
{
  int kstat = 0;        /* Local status variable       */
  double crpr[3], a;    /* Temp. variable.             */


  /* Evaluate torsion */

  s6crss( derive + 3, derive + 6, crpr );

  a = s6length( crpr, 3, &kstat );

  if ( a != 0.0 )
    {
      *torsion = ( s6scpr( derive + 9, crpr, 3 ) ) / ( a*a ) ;
    }
  else
    {
      *torsion = 0.0;

      goto war002;
    }



 /* Successful computations.  */

  *jstat = 0;
  goto out;


  /* Torsion undefined, since the curvature = 0.0.  */
war002:
  *jstat = 2;
  goto out;


 out:

 return;

}
