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
 * $Id: s1393.c,v 1.2 2001-03-19 15:58:49 afr Exp $
 *
 */


#define S1393

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1393(int n1,SISLCurve *pc1[],SISLCurve *sc1[],SISLCurve *ec1[],int *jstat)
#else
void s1393(n1,pc1,sc1,ec1,jstat)
     int   n1;
     SISLCurve *pc1[];
     SISLCurve *sc1[];
     SISLCurve *ec1[];
     int   *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* Purpose :  Split curve at midpoint, turn second curve, 
*            and normalize parameterinterval.
*
* Input     : pc1       - Pointers to first boundary curves
*             n1        - Number of curves
*
* Output    : sc1       - Pointers to first part of curve.
*             ec1       - Pointers to second part of curve.
*
*             jstat     - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*                      
*-
* Calls      : s6err - error messages.
*              s1710 - split curve at given parameter value.
*              s1706 - Turn parametrization of curve.
*              s1399 - Normalize parameterinterval.
*
* Written by : Mortend Daehlen, SI, Aug. 88.
*
*********************************************************************
*/                                     
{
  int kpos = 0;
  int ki;
  int kstat = 0;
  double ax,astart,astop;
  SISLCurve *h1,*h2;
  
  astart = DZERO;
  astop  = (double)1.0;
  
  /* For each curve in pc1 split/turn and normalize. */
  for (ki=0;ki<n1;ki++)
    {
      
      /* Split */
      
      ax=(pc1[ki]->et[pc1[ki]->in]-(pc1[ki]->et[(pc1[ki]->ik)-1]))/(double)2.0;
      s1710(pc1[ki],ax,&h1,&h2,&kstat); 
      if (kstat < 0) goto error;
      
      /* Turn */
      
      s1706(h2);
      if (kstat < 0) goto error;
      
      /* Normalize */
      
      s1399(h1,astart,astop);
      if (kstat < 0) goto error;
      s1399(h2,astart,astop);
      if (kstat < 0) goto error;
      sc1[ki] = h1;
      ec1[ki] = h2;
    }
  
  *jstat=0;
  goto out;
  
  /* Error in lower level routine.   */
  
 error: *jstat = kstat;
  s6err("s1393",*jstat,kpos);
  goto out;
  
 out: return;
}
