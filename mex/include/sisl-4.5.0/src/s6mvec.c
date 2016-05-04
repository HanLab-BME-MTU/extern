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
 * $Id: s6mvec.c,v 1.2 2001-03-19 15:59:02 afr Exp $
 *
 */


#define S6MVEC

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s6mvec(double emat[],double evec1[],int inbvec,double evec2[])
#else
void s6mvec(emat,evec1,inbvec,evec2)
     double emat[];
     double evec1[];
     int    inbvec;
     double evec2[];
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To post multiply a matrix by inbvec 3-D vectors.
*
* INPUT      : emat    - The matrix (16 elements, Homogenous coordinates)
*              evec1   - The input vectors (3-D coordinates)
*              inbvec  - The number of input vectors
*
* OUTPUT     : evec2   - The resulting vectors (3-D coordinates)
*
*-  
* CALLS      :
*
* WRITTEN BY : Tor Dokken, SI, Oslo, Norway. 1988-may-24
*                                  
*********************************************************************
*/
{
  double svec2[3];           /* Temporary vector                        */
  double *svec1;             /* Pointer to current vector               */
  double *svec3;             /* pointer to current result vector        */
  register double tdum;      /* Temporary storage of real               */
  register int ki,kj,kl,kp;  /* Control variables in loops              */
  int kstop;                 /* Stop condition for loop                 */
  
  /* Multiply matrix by all vectors */
  
  kstop = 3*inbvec;
  
  for (kl=0;kl<kstop;kl=kl+3)
    {
      /* Multiply rotational part of the matrix by the vector */
      
      svec1 = evec1 + kl;
      svec3 = evec2 + kl;
      
      for (ki=0;ki<3;ki++)
        {
	  kp = ki;
	  
	  tdum = DZERO;
	  for (kj=0;kj<3;kj++)
            {
	      tdum += emat[kp]*svec1[kj];
	      kp += 4;
            }
	  
	  /* Add translation part */
	  svec2[ki]  = tdum + emat[kp];
        }
      /*  Check if the bottom row is 0,0,0,1 */                     
      
      if (DNEQUAL(emat[3],DZERO) || DNEQUAL(emat[7],DZERO) ||
	  DNEQUAL(emat[11],DZERO) || DNEQUAL(emat[15],(double)1.0))
        {
	  /* Compute last element of vector */
	  
	  tdum = evec1[0]*emat[3] + evec1[1]*emat[7] + evec1[2]*emat[11];
	  if (DNEQUAL(tdum,DZERO))
            {
	      for (ki=0;ki<3;ki++)
		svec2[ki] /= tdum;
            }
	}
      svec3[0] = svec2[0];    
      svec3[1] = svec2[1];
      svec3[2] = svec2[2];
    }
  
  return;
}
