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
 * $Id: s1399.c,v 1.2 2001-03-19 15:58:49 afr Exp $
 *
 */


#define S1399

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1399(SISLCurve *pc,double astart,double astop)
#else
void s1399(pc,astart,astop)
     SISLCurve  *pc;
     double astart;
     double astop;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Change the knotvector to go from astart to astop.
*
* INPUT      : pc      - The curve.
*              astart  - Parametervalue at new startpoint.
*              astop   - Parametervalue at new endpoint.
*
*-
* CALLS      :
*
* WRITTEN BY : Morten Daehlen, SI, 88-09.
*
********************************************************************/
{
  int  kk= pc->ik;             /* Order of the input curve.             */
  int  kn= pc->in;             /* Number of vertices in the input curve.*/
  double *st=SISL_NULL;             /* Pointers used in loop.                */ 
  double a,b;
  int ii, kpos=0, kstat=0;
  if (!pc) goto out;
  
  if((st = newarray(kk+kn,DOUBLE)) == SISL_NULL) goto err101;
  
  a = pc -> et[kk-1];
  b = pc -> et[kn];
  
  for (ii=0;ii<kn+kk;ii++)
    st[ii] = astart+(((pc -> et[ii])-a)/(b-a))*(astop-astart);
  for (ii=0;ii<kn+kk;ii++)
    pc -> et[ii] = st[ii];
  goto out;
  
  /* Error in scratch allocation */
  
  err101: 
    kstat = -101;
    s6err("s1399",kstat,kpos);
    goto out;
      
  out:
    if (st != SISL_NULL) freearray(st);
    return;
}
