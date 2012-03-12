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
 * $Id: s1792.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S1792

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
double 
s1792(double et[],int ik,int in)
#else
double s1792(et,ik,in)
     double et[];
     int    ik;
     int    in;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Finding a good subdividing parametric value.
*
*
*
* INPUT      : et[]     - The knot vector.
*              ik       - The order of the spline function.
*              in       - The number of basic spline function.
*
*
*
* OUTPUT     : s1792    - The subdivision value.
*
*
* METHOD     : This function is written to match s1791().
*
*
* REFERENCES :
*
*
* WRITTEN BY : Arne Laksaa, SI, 89-04.
*
*********************************************************************
*/                                     
{
  if (in > ik)
    {
      int kpar = (in + ik)/2;
      
      if (DNEQUAL(et[ik-1],et[kpar]) || DNEQUAL(et[in],et[kpar]))
        return  et[kpar];
    }
  
  return (et[ik-1]+et[in])*(double)0.5;
}