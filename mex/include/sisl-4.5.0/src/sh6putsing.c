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
 * $Id: sh6putsing.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define SH6PUTSING
#include "sislP.h"                            

#if defined (SISLNEEDPROTOTYPES)
void
      sh6putsing( SISLIntpt *psource, SISLIntpt *pdest)
#else
	 
	 void sh6putsing(psource,pdest)
	    
	    SISLIntpt *psource,*pdest;
#endif

/*
*********************************************************************
*                                                                   
* PURPOSE    : To set some approximative values in a singular point
*              based on symmetry.
*
*
*
* INPUT      : psource - The previous point on the curves.
*              pdest   - The point that is to be fit in.
*                        Posisiton must be ok.
* OUTPUT     : pdest   - Geometric values, ie tangent is changed.

*
* METHOD     : For each curve we mirror the tangent in psource
*              around the difference vector pdest-psource.
*
* REFERENCES :
*
*-
* CALLS      : s6scpr, s6diff, s6norm
*
* WRITTEN BY : UJK
*
*********************************************************************
*/
{
  int kdim,kstat;                
  double alfa;
  double diffv[3];
  double delta[3];
  int ki;
  
  kdim = 3;
  s6diff(pdest->geo_track_3d,psource->geo_track_3d,kdim,diffv);
  s6norm(diffv,kdim,delta,&kstat);
  alfa = (double)2.0*s6scpr(delta,psource->geo_track_3d+kdim,kdim);
  for (ki=0;ki<kdim;ki++) 
    pdest->geo_track_3d[kdim+ki] = alfa*delta[ki] - psource->geo_track_3d[kdim+ki];

  pdest->geo_track_3d[9] = -(double) 1.0;
  
  kdim = 2;
  s6diff(pdest->geo_track_2d_1,psource->geo_track_2d_1,kdim,diffv);
  s6norm(diffv,kdim,delta,&kstat);
  alfa = (double)2.0*s6scpr(delta,psource->geo_track_2d_1+kdim,kdim);
  for (ki=0;ki<kdim;ki++) 
    pdest->geo_track_2d_1[kdim+ki] = alfa*delta[ki] - psource->geo_track_2d_1[kdim+ki];
  
  pdest->geo_track_2d_1[6] = -(double) 1.0;

  kdim = 2;
  s6diff(pdest->geo_track_2d_2,psource->geo_track_2d_2,kdim,diffv);
  s6norm(diffv,kdim,delta,&kstat);
  alfa = (double)2.0*s6scpr(delta,psource->geo_track_2d_2+kdim,kdim);
  for (ki=0;ki<kdim;ki++) 
    pdest->geo_track_2d_2[kdim+ki] = alfa*delta[ki] - psource->geo_track_2d_2[kdim+ki];

  pdest->geo_track_2d_1[6] = -(double) 1.0;

}
