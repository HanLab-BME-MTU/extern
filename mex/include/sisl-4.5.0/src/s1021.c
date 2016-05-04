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
 * $Id: s1021.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S1021

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1021(double bottom_pos[], double bottom_axis[], double ellipse_ratio, 
	 double axis_dir[], double height, SISLSurf **cyl, int *stat)
#else
void s1021(bottom_pos, bottom_axis, ellipse_ratio, axis_dir, 
	      height, cyl, stat) 
     double bottom_pos[];
     double bottom_axis[];
     double ellipse_ratio;
     double axis_dir[];
     double height;
     SISLSurf  **cyl;
     int    *stat;      
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To describe a truncated cylinder as a NURBS. The cylinder 
*              can be elliptic.
*             
*
* INPUT      : bottom_pos    - Center point of the bottom
*              bottom_axis   - One of the bottom axis (major or minor)
*              ellipse_ratio - Ratio betwwen the other axis and bottom_axis
*              axis_dir      - Direction of the cylinder axis
*              height        - Height of the cone, can be negative.
*
*
* OUTPUT     : 
*              stat          - status messages  
*                                             > 0      : warning
*                                             = 0      : ok
*                                             < 0      : error
*              cyl           - Pointer to the cylinder produced
*
* METHOD     : The cylinder is made as a rules surface between two NURBS ellipses.
*
*
* REFERENCES :
*
*-                                      
* CALLS      : 
*
* WRITTEN BY : Johannes Kaasa, SI, Oslo, Norway, Jan. 93
*
*********************************************************************
*/
{
   int kstat;                      /* Status variable.   */
   int kpos=0;                     /* Position of error. */
   double cone_angle = (double)0.; /* No conical slope.  */

   s1022(bottom_pos, bottom_axis, ellipse_ratio, axis_dir, 
	     cone_angle, height, cyl, &kstat);
   if (kstat < 0) goto error;
  
   *stat = 0;
   goto out;
  
   /* Error in lower level routine. */
      
   error:
      *stat = kstat;
      s6err("s1021", *stat, kpos);
      goto out;
  
   out:  
      return;
}
    
