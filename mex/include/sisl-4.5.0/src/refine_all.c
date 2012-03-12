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
 * $Id: refine_all.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define REFINE_ALL

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
refine_all (SISLIntdat ** pintdat,
	    SISLObject * po1,
	    SISLObject * po2,
	    double eimpli[],
	    int ideg,
	    double aepsge,
	    int *jstat)

#else
void
refine_all (pintdat,
	    po1,
	    po2,
	    eimpli,
	    ideg,
	    aepsge,
	    jstat)

     SISLIntdat **pintdat;
     SISLObject *po1;
     SISLObject *po2;
     double eimpli[];
     int ideg;
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
* INPUT      : pintdat     - Pointer to pointer to the SISLIntdat data.
*              po1         - Pointer surface object.
*              po2         - Pointer surface object.
*              eimpli      - Array containing descr. of implicit surf
*	       ideg        - Type of impl surf.
              ang_tol     - Angle control tolerance ie ??
*              aepsge      - Absolute tolerance
*
*
* OUTPUT     :  jstat  - status messages
*                       = ?      : ?
*                       = 0      : ok
*                       < 0      : error
*
*
* METHOD     :
*
* USE        : The function is only working i 3-D
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Ulf J. Krystad, SI, Oslo, Norway, July-1990
*
*********************************************************************
*/
{
  *jstat = 0;
}
