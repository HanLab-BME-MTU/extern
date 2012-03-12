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
 * $Id: s1514.c,v 1.4 2001-03-19 15:58:50 afr Exp $
 *
 */


#define S1514

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1514 (SISLSurf * ps1, double eyepoint[], int idim, double aepsco, double aepsge,
       double amax, SISLIntcurve * pintcr, int icur, int igraph, int *jstat)
#else
void
s1514 (ps1, eyepoint, idim, aepsco, aepsge, amax, pintcr, icur, igraph, jstat)
     SISLSurf *ps1;
     double eyepoint[];
     int idim;
     double aepsco;
     double aepsge;
     double amax;
     SISLIntcurve *pintcr;
     int icur;
     int igraph;
     int *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : To march the perspective silhouette curve described by an intersection
*              curve object, a B-spline surface and an eye point.
*
*
* INPUT      : ps1    - Pointer to surface.
*              eyepoint  - Eye point for perspective view
*              idim   - Dimension of the space in which the eyepoint
*                       lies.
*              aepsco - Computational resolution.
*              aepsge - Geometry resolution.
*              amax   - Maximal allowed step length. If amax <=aepsge
*                       amax is neglected.
*              icur   - Indicator telling if a 3-D curve is to be made
*                        0 - Don't make 3-D curve
*                        1 - Make 3-D curve
*                        2 - Make 3-D curve and curves in parameter plane
*              igraph - Indicator telling if the curve is to be outputted
*                       through function calls:
*                        0 - don't output curve through function call
*                        1 - output as straight line segments through
*                            s6move and s6line.
*
*
*
* INPUT/OUTPUT:pintcr - The intersection curve. When coming in as input
*                       only parameter values in the parameter plane
*                       exist. When coming as output the 3-D geometry
*                       and possibly the curve in the parameter plane
*                       of the surface is added.
*
* OUTPUT:      jstat  - status messages
*                         = 3      : Iteration stopped due to singular
*                                    point or degenerate surface. A part
*                                    of intersection curve may have been
*                                    traced out. If no curve is traced out
*                                    the curve pointers in the Intcurve
*                                    object point to SISL_NULL.
*                         = 0      : ok
*                         < 0      : error
*                         = -185   : No points produced on intersection curve.
*
*
* METHOD     : An implicit description of the problem is made and then
*              a routine for intersecting implicit represented geometry
*              by a B-spline surface is used.
*
* REFERENCES :
*
*-
* CALLS      : s6err, s1313
* WRITTEN BY : Mike Floater, SI, Oslo, Norway, 31 . January 1991
*                a modification of s1319
*
*********************************************************************
*/
{
  int kpos = 0;			/* Position of error                                  */
  int kdeg = 1004;		/* The degree of the implicit equation                */
  /* of the perspective silhouette                      */
  int kstat;			/* Local status variable                              */
  double simpli[3];		/* Array containing the implicit description          */
  /* of the perspective silhouette                      */



  if (idim != 3)
    goto err104;


  simpli[0] = eyepoint[0];
  simpli[1] = eyepoint[1];
  simpli[2] = eyepoint[2];

  /* Make intersection of implicit surface and B-spline surface */

  s1313 (ps1, simpli, kdeg, aepsco, aepsge, amax, pintcr, icur, igraph, &kstat);
  if (kstat == -185)
    goto err185;
  if (kstat < 0)
    goto error;

  *jstat = kstat;
  goto out;

  /* Dimension not 3 */

err104:
  *jstat = -104;
  s6err ("s1514", *jstat, kpos);
  goto out;

  /* Couldn't march */

err185:
  *jstat = -185;
  goto out;

  /* Error in lower level routine.  */

error:
  *jstat = kstat;
  s6err ("s1514", *jstat, kpos);
  goto out;

out:
  return;
}
