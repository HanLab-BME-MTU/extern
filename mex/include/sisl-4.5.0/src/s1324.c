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
 * $Id: s1324.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S1324

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1324(double ecentr[],double aradiu,double enorm[],int idim,
	   double carray[],int *jstat)
#else
void s1324(ecentr,aradiu,enorm,idim,carray,jstat)
     double ecentr[];
     double aradiu;
     double enorm[];
     int    idim;
     double carray[];
     int    *jstat;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To make two matrix of dimension 4x4
*              describing a 3-D circle as two implicit functions.
*
*
* INPUT      : ecentr - Center of the circle
*              aradiu - Radius of the circle
*              enorm  - Normal vector of circle plane
*              idim   - The dimension of the space the cirle lies in
*
*
*
* OUTPUT     : carray - The description of the circle. Outside
*                       this function the space for this array must be
*                       allocated. The need is 32 double variables.
*                       First the matrix for the sphere is stored,
*                       then the matrix of the plane.
*              jstat  - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : The circle is described as an intersection between a
*              cylinder and the plane. The matrix describing the
*              cylinder is put first in the output array, the matrix
*              describing the plane follows then.
*              
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Tor Dokken, SI, Oslo, Norway, 29-June-1988
*                                    
*********************************************************************
*/
{
  int kdimp1;         /* Dimension of matrix kdimp1 = idim + 1         */
  int kdimp2;         /* idim + 2                                      */
  int kstop;          /* Stop condition for for loop                   */
  int ki;             /* Running variables in loop                     */
  int kpos=0;         /* Position of error                             */
  int kstat;          /* Status variable                               */
  
  
  
  /* Test i legal input */
  if (idim != 3) goto err104;
  
  kdimp1 = idim + 1;
  kdimp2 = idim + 2;
  kstop  = 2*kdimp1*kdimp1;
  
  for (ki=0;ki<kstop;ki++)
    {
      carray[ki] = (double)0.0;
    }
  
  /* Make description of cylinder */
  
  s1322(ecentr,enorm,aradiu,idim,1,carray,&kstat);
  if (kstat<0) goto error;
  
  
  /* Make description of plane, element (1,4), (2,4) and (3,4) */
  
  carray[28] = enorm[0];
  carray[29] = enorm[1];
  carray[30] = enorm[2];
  
  /* Make element (4,4) */
  
  carray[31] = -s6scpr(enorm,ecentr,idim);
  
  *jstat = 0;
  goto out;
  
  /* Dimension not 3 */
 err104: *jstat = -104;
  s6err("s1324",*jstat,kpos);
  goto out;
  
  /* Error in lower level routine */
 error: *jstat = kstat;
  goto out;
  
  
 out:
  return;
}
