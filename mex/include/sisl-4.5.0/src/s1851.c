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
 * $Id: s1851.c,v 1.3 2001-03-19 15:58:54 afr Exp $
 *
 */


#define S1851

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void s1851(SISLSurf *ps1,double epoint[],double enorm[],int idim,
	   double aepsco,double aepsge,int *jpt,double **gpar,
	   int *jcrv,SISLIntcurve ***wcurve,int *jstat)
#else
void s1851(ps1,epoint,enorm,idim,aepsco,aepsge,
	   jpt,gpar,jcrv,wcurve,jstat)
     SISLSurf     *ps1;
     double   epoint[];
     double   enorm[];
     int      idim;
     double   aepsco;
     double   aepsge;
     int      *jpt;
     double   **gpar;
     int      *jcrv;
     SISLIntcurve ***wcurve;
     int      *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Find all intersections between a tensor-product surface
*              and a plane.
*
*
*
* INPUT      : ps1    - Pointer to surface.
*              epoint - SISLPoint in the plane.
*              enorm  - Normal to the plane.
*              idim   - Dimension of the space in which the plane lies.
*              aepsco - Computational resolution.
*              aepsge - Geometry resolution.
*
*
*
* OUTPUT     : *jpt   - Number of single intersection points.
*              gpar   - Array containing the parameter values of the
*                       single intersection points in the parameter
*                       plane of the surface. The points lie continuous.
*                       Intersection curves are stored in wcurve.
*              *jcrv  - Number of intersection curves.
*              wcurve  - Array containing descriptions of the intersection
*                       curves. The curves are only described by points
*                       in the parameter plane. The curve-pointers points
*                       to nothing. (See description of Intcurve
*                       in intcurve.dcl).
*              jstat  - status messages
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : The vertices of the surface are put into the equation of the
*              plane achieving a surface in the one-dimentional space.
*              Then the zeroes of this surface is found.
*
*
* REFERENCES : Main routine written by Vibeke Skytt, SI, 1988.
*
* CALLS      : sh1851, s6err.
*
* WRITTEN BY : Christophe Rene Birkeland, SINTEF, 93-06.
* Revised by : Paal Fugelli, SINTEF, Oslo, Norway, Oct. 1994.
*              Initialized 'jsurf'.
*
*********************************************************************
*/
{
  int kstat = 0;              /* Local status variable.                      */
  int kpos = 0;               /* Position of error.                          */
  int i;
  int trackflag = 0;
  int jtrack;
  SISLTrack **wtrack=SISL_NULL;
  int jsurf = 0;
  SISLIntsurf **wsurf=SISL_NULL;
  int *pretop=SISL_NULL;

  sh1851(ps1, epoint, enorm, idim, aepsco, aepsge,trackflag, &jtrack,
	 &wtrack,jpt, gpar,&pretop,jcrv,wcurve,&jsurf,&wsurf,&kstat);
  if(kstat < 0) goto error;

  if(pretop != SISL_NULL) freearray(pretop);

  for(i=0; i<jsurf; i++)
    freeIntsurf(wsurf[i]);
  if(wsurf != SISL_NULL) freearray(wsurf);

  if(jsurf > 0)
    *jstat=10;
  else
    *jstat = 0;
  goto out;

  /* Error in lower level routine.  */

  error :
    *jstat = kstat;
    s6err("s1851",*jstat,kpos);
    goto out;

  out:
    return;
}
