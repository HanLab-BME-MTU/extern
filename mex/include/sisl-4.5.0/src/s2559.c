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

#define S2559

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s2559( SISLCurve *curve,
       double     ax[],
       int        num_ax,
       double     p[],
       double     t[],
       double     n[],
       double     b[],
       int       *jstat )
#else
void s2559( curve, ax, num_ax, p, t, n, b, jstat )
     SISLCurve  *curve;
     double      ax[];
     int         num_ax;
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
*              at given parameter values ax[ 0 ],...,ax[ num_ax - 1 ].
*
*
*
* INPUT      : curve    - Pointer to the curve.
*              ax       - The parameter values
*          num_ax       - No. of parameter values
*
*
*
*
* OUTPUT     :
*         (t,n,b) - The Frenet Frame (in 3D) computed. Each of the arrays
*                   (t,n,b) are of dim. 3*num_ax, and the data are
*                   stored like this: tx(ax[0]), ty(ax[0]), tz(ax[0]),
*                   ...,tx(ax[num_ax-1]), ty(ax[num_ax-1]), tz(ax[num_ax-1]).
*
*         p     - 3D curve posistions at  ax[ 0 ],...,ax[ num_ax - 1 ]
*
*        jstat  - Status messages
*                                         > 0      : Warning.
*                                         = 0      : Ok.
*                                         < 0      : Error.
*
*
* METHOD     : Call the SISL rutine s2560() num_ax times
*
* REFERENCES :
*
*
* CALLS      : s2560()
*
* WRITTEN BY  :   Geir Westgaard, SINTEF, Oslo, November 1999
* MODIFIED BY :
* REVISED BY  :
*
*********************************************************************
*/
{
  int i, m;                   /* Index variable                      */
  int kstat = 0;              /* local status variable               */
  int kpos = 0;               /* local error position                */
  int leftknot = 0;           /* Knot index                          */
  double *derive = SISL_NULL;      /* Derivatives                         */



  /* Allocate local arrays */

  derive = newarray( 3*curve->idim, DOUBLE );

  if ( derive == SISL_NULL )      goto err101;


  /* Evaluate the Frenet Frame in all ax[i] positions. */

  m = 0;
  for ( i = 0; i < num_ax; i++ )
  {
    s2560( curve, ax[i], &leftknot, derive,
	   p + m, t + m, n + m, b + m, &kstat );

    m = m + 3;

    if ( kstat < 0 ) goto error;
  }




 /* Successful computations.  */

  *jstat = 0;
  goto out;

  /* Error in scratch allocation */
 err101:
  *jstat = -101;
  s6err( "s2559", *jstat, kpos );

  /* Error in lower level routine.  */

 error:  *jstat = kstat;
  s6err( "s2559", *jstat, kpos );
  goto out;


 out:
  /* Free local arrays */

 if ( derive != SISL_NULL ) freearray( derive );

 return;

}
