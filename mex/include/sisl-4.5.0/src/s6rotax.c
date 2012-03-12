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
 * $Id: s6rotax.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S6ROTAX

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s6rotax(double ep[],double eaxis[],double expnt[],double emat[],int *jstat)
#else
void s6rotax(ep,eaxis,expnt,emat,jstat)
     double ep[];
     double eaxis[];
     double expnt[];
     double emat[];
     int    *jstat;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To make a matrix that transforms the z-axis ont to
*              the axis specified by ep and eaxis, and transforms
*              the point (1,0,0) onto the point expnt. The matrix is
*              prepared for post multiplication of vectors.
*
* INPUT      : ep      - SISLPoint on axis
*              eaxis   - Direction of axis 
*              expnt   - The point (1,0,0) is trnasformed on to.
*
* OUTPUT     : jstat   - Status message
*                                        >0      : Warning
*                                        =0      : ok
*                                        <0      : Error
*              emat    - The 4x4 matrix produced.
*
*
* METHOD     : 
*
*-
* CALLS      :
*
* WRITTEN BY : Tor Dokken, SI, Oslo, Norway. 1988-may-23
*
*********************************************************************
*/
{
  double sdiff[3];            /* Vector for storage of  expnt-ep        */
  double saxis[3];            /* Storage of normalized eaxis            */
  double tradius;             /* Distance from expnt to axis            */
  double sxdir[3];            /* Transformation direction of x-axis     */
  double sydir[3];            /* Transformation direction of y-axis     */
  double strans[3];           /* Translation vector for the origin      */
  double zfak;                /* Length of projection                   */
  int    kdim=3;              /* The dimension of the space we work in  */
  int    kstat;               /* Local status varaible                  */
  int    ki;                  /* Variable used in for loop              */
  
  
  /* Normalize direction of axis */
  
  (void)s6norm(eaxis,kdim,saxis,&kstat);
  
  
  /*  Make difference of expnt  and ep  */
  
  for (ki=0;ki<3;ki++)
    sdiff[ki] = expnt[ki] - ep[ki];
  
  /* Find projection of expnt-ep onto saxis */
  
  zfak = s6scpr(sdiff,saxis,kdim);
  
  
  /* Make transformation of the vector (1,0,0) */
  
  for (ki=0;ki<3;ki++)
    sxdir[ki] = sdiff[ki] - zfak*saxis[ki];

  /* Normalize sxdir */
  
  tradius = s6norm(sxdir,kdim,sxdir,&kstat);
  
  /* Make the vector (0,1,0) is to be projected onto */
  
  s6crss(saxis,sxdir,sydir);
  (void)s6norm(sydir,kdim,sydir,&kstat);
  
  /* Make translation of origo */
  
  for (ki=0;ki<3;ki++)
    strans[ki] = ep[ki] + zfak*saxis[ki];
  
  /* Build transformation matrix for post multiplication of vectors */
  
  emat[0]  = tradius*sxdir[0];
  emat[1]  = tradius*sxdir[1];
  emat[2]  = tradius*sxdir[2];
  emat[3]  = (double)0.0;
  emat[4]  = tradius*sydir[0];
  emat[5]  = tradius*sydir[1];
  emat[6]  = tradius*sydir[2];
  emat[7]  = (double)0.0;
  emat[8]  = tradius*saxis[0];
  emat[9]  = tradius*saxis[1];
  emat[10] = tradius*saxis[2];
  emat[11] = (double)0.0;
  emat[12] = strans[0];
  emat[13] = strans[1];
  emat[14] = strans[2];
  emat[15] = (double)1.0;
  
  /* Matrix made */
  
  *jstat = 0;
  goto out;
  
 out: 
  return;
}
