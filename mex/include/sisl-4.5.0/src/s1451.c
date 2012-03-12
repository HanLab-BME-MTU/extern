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
 * $Id: s1451.c,v 1.2 2005-02-28 09:04:48 afr Exp $
 *
 */


#define S1451

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1451(SISLCurve *pc1,double aepsge,int *jdgen,int *jstat)
#else
void s1451(pc1,aepsge,jdgen,jstat)
     SISLCurve  *pc1;
     double aepsge;
     int    *jdgen;
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : To check if a curve is degenerated.
*                              
* INPUT      : pc1    - Pointer to the curve to be tested.
*              aepsge - The curve is degenerate if all vertices lie
*                       within the distance aepsge from each other
*
* OUTPUT     : jdgen  - Degenerate indicator
*                        0 - The curve is not degenerate
*                        1 - The curve is degenerate
*
*              jstat  - status messages  
*                                         > 0      : Warning
*                                         = 0      : ok
*                                         < 0      : error
*                      
*
* METHOD     : Testing  distance between all vertices
*
* USE:         int kdgen,kstat;
*              double aepsge;
*              SISLCurve *qc1;
*               .
*               .
*              s1451(qc1,aepsge,&kdgen,&kstat);
*               .
*
*
* REFERENCES :
*
*-
* CALLS      : s6dist, s6err; 
*
* WRITTEN BY : Tor Dokken, SI, Norway, 1988-11
*
*********************************************************************
*/                                     
{
  int kpos=0;         /* Position of error.                              */
  int kdim;           /* Dimension of the space in which the surface lies. */
  int kn;             /* Number of vertices                              */
  int kk;             /* Polynomial order                                */
  
  int ki,kj;          /* Control variables in for loop                   */
  double *scoef;      /* Vertices                                        */
  double *sv1,*sv2;   /* Pointer to vertices                             */
  
  
  if (aepsge < (double)0.0) goto err184;
  
  /* Initiate degenerate indicator */
  
  *jdgen = 1;
  
  kn  = pc1->in;
  kk  = pc1->ik;
  kdim = pc1 -> idim;
  scoef  = pc1 -> ecoef;
  
  for (kj=0,sv2=scoef ; kj < kn ; kj++,sv2+=kdim)
    {
      for (ki=kj+1,sv1=sv2+kdim; ki < kn ; ki++,sv1+=kdim)
	{
	  if (s6dist(sv1,sv2,kdim) > aepsge)
	    {
	      *jdgen = 0;
	      ki = kn;
	      kj = kn;
	    }
	}
    }
  
  *jstat = 0;
  goto out;
  
  /* Negative absolute tolerance.   */
  
 err184: *jstat = -184;
  s6err("s1451",*jstat,kpos);
  goto out;
  
 out:
  return;
}

