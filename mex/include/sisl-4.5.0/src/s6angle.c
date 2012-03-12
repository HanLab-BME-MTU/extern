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
 * $Id: s6angle.c,v 1.2 2001-03-19 15:59:00 afr Exp $
 *
 */


#define S6ANGLE

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
double
     s6angle(double evec1[],double evec2[],double enorm[],int idim,int *jstat) 
#else
double s6angle(evec1,evec2,enorm,idim,jstat)
     double evec1[];
     double evec2[];
     double enorm[];
     int    idim;
     int *jstat;
#endif     
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Compute the angle (in radians) between two vectors
*              projected down into a given plane.
*
*
*
* INPUT      : evec1   - First vector 
*              evec2   - Second vector 
*              enorm   - Normal of plane
*              idim    - Dimension of the space in which the vectors lie.
*
*
*
* OUTPUT     : s6ang   - Angle in radians between vectors
*              jstat   - Status messages
*
*
* METHOD     : Make cosine of the angle by computing the scalar product,
*              then divide by the length of the two vectors.
*
* REFERENCES :
*
*-
* CALLS      : s6scpr   - Scalar product between two vector.
*              s6length - Length of vector.
*              s6crss   - Cross product between two vector.
*
* WRITTEN BY : Tor Dokken SI, 88-07.
*              Arne Laksaa SI, 89-07.
*              Vibeke Skytt, SI, 90-08.
*
*********************************************************************
*/                                     
{
  double sa[3];
  double sb[3];
  double sn[3];
  
  double tscpr1,tscpr2,tang,tlength1,tlength2,tcos;
  int    kstat1,kstat2,ki;
    
  if (idim != 3) goto err104;
  
  tscpr1 = s6scpr(evec1,enorm,idim);
  tscpr2 = s6scpr(evec2,enorm,idim);
  
  for (ki=0; ki<idim; ki++)
    {
      sa[ki] = evec1[ki] - tscpr1*enorm[ki];
      sb[ki] = evec2[ki] - tscpr2*enorm[ki];
    }
  
  tscpr1 = s6scpr(sa,sb,idim);

  tlength1 = s6length(sa,idim,&kstat1);
  tlength2 = s6length(sb,idim,&kstat2);
  
  if (!kstat1 || !kstat2)
    tang = DZERO;
  else
    {
      tcos = tscpr1/(tlength1*tlength2);
      tcos = MIN((double)1.0,tcos);
      tcos = MAX(-(double)1.0,tcos);
      tang = acos(tcos);
    }
  
  s6crss(sa,sb,sn);
  if (s6scpr(sn,enorm,idim) < DZERO) tang = TWOPI - tang;
  
  *jstat = 0;
  goto out;
  
 err104:
  *jstat = -104;
  goto out;
  
 out:
  return(tang);
}
