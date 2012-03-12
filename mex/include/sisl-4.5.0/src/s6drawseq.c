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
 * $Id: s6drawseq.c,v 1.3 1994-12-19 16:58:19 pfu Exp $
 *
 */


#define S6DRAWSEQ

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
extern void s6move(DOUBLE[]);
extern void s6line(DOUBLE[]);
#else
extern void s6move();
extern void s6line();
#endif

#if defined(SISLNEEDPROTOTYPES)
void
  s6drawseq(double epoint[],int ipoint)
#else
void s6drawseq(epoint,ipoint)
     double epoint[];
     int ipoint;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Draw a broken line as a sequence of straight lines
*              described by the array epoint.
*
*
*
* INPUT      : epoint - Array describing the corners of the broken line
*                       to draw. Each corner is given as a 3-dimensional
*                       point. The corners is placed continuous in the
*                       array.
*              ipoint - Number of points in epoint.
*
*
* OUTPUT     :
*
*
* METHOD     :
*
*
* REFERENCES :
*
*
* REMARK     : This routine is machine-dependant and the internal of
*              this routine has to be reprogrammed if the functions
*              s6move() and s6line() are not defined for your graphics
*              sub-system.
*
*
*-
* CALLS      : s6move - Place pen at given position (empty dummy routine).
*              s6line - Draw a line from the current position to the given one
*                       (empty dummy routine).
*
* WRITTEN BY :
*
*********************************************************************
*/
{
  int ki;          /* Counter.                                    */
  double *spoint;  /* Pointer to corner point in the broken line. */

  /* Position pen at start of the broken line.  */

  s6move(epoint);

  /* Draw sequence of line-segments.  */

  for (ki=1,spoint=epoint+3; ki<ipoint; ki++,spoint+=3)
    s6line(spoint);

  return;
}
