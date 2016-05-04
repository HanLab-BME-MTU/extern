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
 * $Id: s6equal.c,v 1.2 2001-03-19 15:59:01 afr Exp $
 *
 */


#define S6EQUAL

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
int 
s6equal(double a1,double a2,double aref)
#else
int s6equal(a1,a2,aref)
     double a1;
     double a2;
     double aref;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Test if two numbers are equal.
*
*
*
* INPUT      : a1     - First number.
*              a2     - Second number.
*              aref   - Reference value.
*
*
*
* OUTPUT     : s6equal - Tells if numbers are equal.
*                        = 1 : Equal.
*                        = 0 : Not equal.
*
*
* METHOD     :
*
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Vibeke Skytt, SI, 88-11.
*
*********************************************************************
*/
{
  double tval;   /* Number used to test equality.  */
  
  tval = a1 - a2;
  tval += aref;
  tval -= aref;
  
  return(DEQUAL(tval,DZERO));
}


