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
 * $Id: sh6ptouchp.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define SH6PUTTOUCH
#include "sislP.h"                            

#if defined (SISLNEEDPROTOTYPES)
void
      sh6puttouch( SISLIntpt *psource, SISLIntpt *pdest, int seq)
#else
	 
	 void sh6puttouch(psource,pdest, seq)
	    
	    SISLIntpt *psource,*pdest;
	    int seq;
#endif

/*
*********************************************************************
*                                                                   
* PURPOSE    : To set right direction in a touch point (Singular, but no branch).
*              The invariant to hold is that the direction on the curve 
*              in the parameter space is parallel to the delta vector between
*              the to points.
*
*
*
* INPUT      : psource - The previous point on the curves.
*              pdest   - The point that is to be fit in.
*                        Posisiton must be ok.
*              seq     - The sequncing of the points
*                        +1 - psource comes before pdest.
*                        -1 - psource comes after pdest.
* OUTPUT     : pdest   - Geometric values, ie tangent is changed.

*
* METHOD     : 
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
  int kdim2=2;
  int kdim3=3;
  double diffv[3];
  int ki;
  double dot;
  /* ______________________ */
  
  if (psource->iinter == SI_ORD)
    sh6putsing(psource, pdest);
  else
    {
       kdim2 = 2;
       s6diff(pdest->geo_track_2d_1,psource->geo_track_2d_1,kdim2,diffv);
       dot = s6scpr(pdest->geo_track_2d_1+kdim2, diffv, kdim2);
       if (dot * seq < 0)
	 {
	    /* Turn direction */
	    for (ki=0;ki<kdim2;ki++) 
	      {
		 pdest->geo_track_2d_1[kdim2+ki] = -pdest->geo_track_2d_1[kdim2+ki];
		 pdest->geo_track_2d_2[kdim2+ki] = -pdest->geo_track_2d_2[kdim2+ki];
	      }  
	    
	    for (ki=0;ki<kdim3;ki++) 
	      pdest->geo_track_3d[kdim3+ki] = -pdest->geo_track_3d[kdim3+ki];
	 }
       
    }
}
