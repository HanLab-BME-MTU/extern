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
 * $Id: s1994.c,v 1.2 2001-03-19 15:58:59 afr Exp $
 *
 */


#define S1994

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1994(SISLSurf *s1,int *jstat)
#else
void s1994(s1,jstat)
     SISLSurf *s1;
     int  *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Check if a point-surface intersection in one dimention
*              is a simple case,
*              i.e. the intersection will result in one single point.
*
*
*
* INPUT      : s1    - Surface in the intersection problem.
*
*
*
* OUTPUT     : jstat  - status messages  
*                                         = 1      : simpel case.
*                                         = 0      : not simpel case.
*                                         < 0      : error.
*
*
* METHOD     :
*
*
* REFERENCES :
*
* CALLS      : 
*
* WRITTEN BY : TDO SI, 89-08.
*
*********************************************************************
*/
{
  register int ki,kj;
  int kk1, kk2, kn1, kn2;
  int kbez;
  
  double tmaxt, tmaxs;
  double tmint, tmins;
  double tdiff;
  double *scoef=SISL_NULL;
  double noice = (double)100.0 * REL_COMP_RES;   /* Noice killer */ 
  
  /* Init to  simple case. */
  *jstat = 1;
  
  tmaxt = tmaxs = - HUGE;
  tmint = tmins =   HUGE;
  
  /* Get surface attributes. */
  kk1  = s1->ik1;
  kk2  = s1->ik2;
  kn1  = s1->in1;
  kn2  = s1->in2;
  kbez = (kk1 == kn1) && (kk2 == kn2); 
  
  
  /* If the surface is linear in some direction it is simpel case. */
  if ((kk1 == 2 && kn1 == 2) || (kk2 == 2 && kn2 == 2)) goto out;
  
  
  /* Run through vertices in first parameter direction to find
     intervall of first derivative. */
  
  for (kj = 0, scoef = s1->ecoef;kj < kn2; kj++,scoef++)
    for (ki = 1; ki < kn1; ki++,scoef++ )
      {
	tdiff = *(scoef + 1) - *scoef;
	tmint = min(tmint,tdiff);
	tmaxt = max(tmaxt,tdiff);
      }
  
  /* Run through vertices in second parameter direction to find
     intervall of first derivative. */
  for (ki = 0 ;ki < kn1; ki++)
    for (kj = 1 , scoef = s1->ecoef + ki; kj < kn2; kj++, scoef +=kn1 )
      {
	tdiff = *(scoef + kn1) - *scoef;
	tmins = min(tmins,tdiff);
	tmaxs = max(tmaxs,tdiff);
      }

  /* ALA and UJK 30.10.90, remove noice near by zero */
  
  if (fabs(tmint) < noice) tmint = DZERO; 
  if (fabs(tmaxt) < noice) tmaxt = DZERO; 
  if (fabs(tmins) < noice) tmins = DZERO; 
  if (fabs(tmaxs) < noice) tmaxs = DZERO; 


  
  
  /* The first derivatives decide directions of possible intersection curves. */
  if (kbez && (tmint*tmaxt >=DZERO || tmins*tmaxs >=DZERO)) 
    *jstat = 1;
  else if (tmint*tmaxt > DZERO || tmins*tmaxs > DZERO) 
    *jstat = 1;
  else if (tmint == tmaxt  || tmins == tmaxs) 
    *jstat = 1;
  else
    /* Not a simple case. */
    *jstat = 0;
  
  goto out;
 out: ;
}

