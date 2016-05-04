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
 * $Id: maketracks.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define MAKE_TRACKS
#include "sislP.h"

#if defined (SISLNEEDPROTOTYPES)
void make_tracks (SISLObject * po1, SISLObject * po2, int ideg,
		  double eimpli[], int icrv, SISLIntlist ** vlist,
		  int *jtrack,SISLTrack *** wcrv, double aepsge, int *jstat)
#else
void make_tracks (po1, po2, ideg,eimpli, icrv, vlist,
		  jtrack,wcrv, aepsge, jstat)
     SISLObject *po1;
     SISLObject *po2;
     int ideg;
     double eimpli[];
     int icrv;
     SISLIntlist **vlist;
     int *jtrack;
     SISLTrack ***wcrv;
     double aepsge;
     int *jstat;
#endif
/*
*********************************************************************
*
* PURPOSE    : An empty function, (to ensure similarity with other
*              versions).
*
*
* INPUT      : po1    - Pointer first object.
*              po2    - Pointer second object.
*                       internal format.
*              icrv   - Number of lists in vlist.
*              vlist  - Array representing intersection curves on the
*                       internal format.
*              ideg   - Type of track
*                           = 0, Bspline vs Bspline
*                           = 1, Bspline vs Plane
*                           = 2, Bspline vs Quadric surface
*                           = 1001 Bspline vs Torus surface
*                           = 1003 Bspline silhouette line, parallel projection
*                           = 1004 Bspline silhouette line, perspective projection
*                           = 1005 Bspline silhouette line, circular projection
*
*              eimpli[16]  Description of the implicit surface.
*              aepsge - Geometry tolerance
*
*OUTPUT:       jtrack - No of tracks made.
*              wtrack - Array containing pointers to tracks.
*              jstat - status messages
*                       >0:warning
*                       = 0:ok
*                       <0:error
*
*
*
*METHOD:   Refining datapoints and
*          Hermite interpolation using curvature radius.
*
*REFERENCES:
*
*
*-
*CALLS:     control_and_refine_ss(_si)
*          s1359 - Hermite interpolation of curve using curvature radius.
*
*WRITTEN BY:Ulf J.Krystad, SI, 30.06 .91.
*
*********************************************************************
*/
{

  *jstat  = 0;
  *jtrack = 0;

}

