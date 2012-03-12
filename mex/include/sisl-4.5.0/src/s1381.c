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
 * $Id: s1381.c,v 1.2 2001-03-19 15:58:48 afr Exp $
 *
 */


#define S1381

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1381(double et[],int in,int ik,double **gt,int *jkn,int jkk,int *jstat)           
#else
void s1381(et,in,ik,gt,jkn,jkk,jstat)
     double et[];
     int    in;
     int    ik;
     double **gt;
     int    *jkn;
     int    jkk;
     int    *jstat;           
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Make the knot vector for the representing a spline
*              basis of order jkk with the same knot values as et.
*              This basis is used e.g. for representing the silhouette line
*              functions of a surface. 
                   jkk = 2*ik-2  for ordinary silhouette
                   jkk = 3*ik-3  for perspective silhouette
                   jkk = 3*ik-3  for circular silhouette
*
* INPUT      : et     - Knots of input spline basis
*              in     - Number of vertices in input basis
*              ik     - Order of input basis
*              jkk    - Order of B-spline basis produced
*
* OUTPUT     : gt     - Pointer to array of knots. The array is allocated
*                       inside this routine.
*              jkn    - Number of vertices
*
*              jstat  - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
* METHOD     : 
*
* REFERENCES :
*
*-                                   
* CALLS      : 
*
* WRITTEN BY : Tor Dokken, SI, 88-11.
* REVISED BY : Mike Floater, SI, 91-04 jkk is input not output.
*
*********************************************************************
*/                                                               
{                                                                     
  double tval;     /* Value of knot                                 */
  double *sdum;    /* Pointer to knot array                         */
  int ki,kl;       /* Variable in loop                              */
  int knumb;       /* Number of intervals                           */
  int kstop;       /* Loop stop variable                            */
  int kpos=0;      /* Position of error                             */
  
  /* Run through the knot vector to decide how many intervals exist */
  knumb = 0;       
  tval = et[ik-1];
  
  for (ki=ik ; ki<=in ; ki++)
    {
      if (tval < et[ki])
        {
	  /*      New knot value found */
	  knumb = knumb + 1;
	  tval = et[ki];
        }
    }
  
  *jkn = (jkk-1)*(knumb-1) + jkk;
  
  sdum = newarray(*jkn+jkk,DOUBLE);
  if (sdum == SISL_NULL) goto err101;
  
  *gt  = sdum; 
  
  /* Make knot values */
  
  tval = et[ik-1];
  
  /* Make jkk first knot values */
  
  for (kl=0;kl<jkk;kl++)
    {        
      sdum[kl] = tval;
    }
  
  /* kl points to the array entry where the next knot value is to be stored
   */
  
  for (ki=ik ; ki<=in ; ki++)
    {
      if (tval < et[ki])
        {
	  /* New knot value, remember this and make knots */
	  tval = et[ki];
	  kstop = kl + jkk-1;
	  for (;kl<kstop;kl++)
            sdum[kl] = tval;
        }   
    }
  
  /* Make last knot value */
  
  sdum[kl] = tval;
  
  *jstat = 0;
  goto out;
  
  /* Error in space allocation */
 err101: *jstat = -101;
  s6err("s1381",*jstat,kpos);
  goto out;
 out:
  
  return;
}                                               
