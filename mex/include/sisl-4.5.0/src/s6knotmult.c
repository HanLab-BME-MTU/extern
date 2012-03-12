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
 * $Id: s6knotmult.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S6KNOTMULT

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
int 
s6knotmult(double et[],int ik,int in,int *ileft,double ax,int *jstat)
#else
int s6knotmult(et,ik,in,ileft,ax,jstat)
     double et[];
     int    in;
     int    ik;
     int    *ileft;
     double ax;
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : To determine the multiplicity of a knot at a specified
*              parameter value.
*
*
* INPUT      : et     - Double array of dimension [in+ik] containing
*                       the knot vector.
*              ik     - The polynomial order of the B-splines associated
*                       with et.
*              in     - The dimension of the spline space associated with
*                       the knot vector et.
*              ax     - The point at which the B-spline values and derivatives
*                       are to be computed.
*
*                
*
* INPUT/OUTPUT :
*              ileft - Pointer to the interval in the knot vector
*                       where ax is located, check the relations above.
*              jstat  - Status messages  
*                                         > 0      : Warning.
*                                         = 0      : Ok.
*                                         < 0      : Error.
*
*
* METHOD     : The aim is to do as little work as possible in the cases
*              where ileft has the right or almost the right value.
*              First of all we make sure that ileft has a legal value
*              (a value in the range ik-1 to in-1). Then we check
*              if the current value is OK.
*              If it is not we check that ax is in the interior of et
*              or if the right value is obtained by either increasing
*              or decreasing ileft by 1. If the right value still has
*              not been found we do a binary search.
*
*
* REFERENCES :
*
*-
* CALLS      : s1219
*
* WRITTEN BY : Knut Moerken, University of Oslo, August 1988.
*
*********************************************************************
*/                                     
{
  int kpos=0;         /* The position of the error.                      */
  int kstat;          /* Local status variable                           */
  int kmult=0;        /* Multiplicity of knot                            */
  int ki;             /* Loop variable                                   */
  
  /* Localize knot interval */
  
  s1219(et,ik,in,ileft,ax,&kstat);
  if (kstat<0) goto error;
  
  if (et[*ileft] == ax)
    {
      kmult = 1;
      ki    = *ileft-1;
      for (ki=(*ileft)-1; 0 <= ki; ki--)
        if (et[ki] == ax) kmult++;
    }
  if (et[in] == ax)
    {
      for (ki=in ; ki<in+ik;ki++)
        if (et[ki] == ax) kmult++;
    }
  
  *jstat = 0;
  goto out;

/* Error in lower level function */

error:  *jstat = kstat;
        s6err("s6knotmult",*jstat,kpos);
        goto out;

out:
return(kmult);
}
