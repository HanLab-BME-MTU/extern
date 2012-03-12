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
 * $Id: s6schoen.c,v 1.2 2001-03-19 15:59:02 afr Exp $
 *
 */



#define S6SCHOEN

#include "sislP.h"
#if defined(SISLNEEDPROTOTYPES)
double
s6schoen(double et[], int ik, int index)
#else
double s6schoen(et,ik,index)
     double et[];
     int    ik;
     int    index;
#endif

/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : To determine the knot value of a specified vertice.
*
*
* INPUT      : et     - Double array of dimension [in+ik] containing
*                       the knot vector.
*              ik     - The polynomial order of the B-splines associated
*                       with et.
*              index  - The vertice index at where the knot values are to 
*                       be computed.
*
*                
*
* INPUT/OUTPUT :
*              s6schoen - The knot value at index.
*
*
* METHOD     : The aim is to calculate the knot value of a vertice, using the
*              Schoenberg spline expression:
*
*               *
*              t  = (t    + .............. +t     )/k-1.
*               i     i+1                    i+k-1
*
* REFERENCES :
*
*-
* CALLS      : 
*
* WRITTEN BY : Per Evensen,SI, August 1991.
*
*********************************************************************
*/                                     
{
  int i;             /* Loop variable                                   */
  double kval=DZERO; /* knot value variable                             */
  
  for (i=index+1;i<index+ik;i++) kval+=et[i];
  kval = kval/(ik-1);

return(kval);
}

