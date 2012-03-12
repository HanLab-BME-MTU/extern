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
 * $Id: s1993.c,v 1.2 2001-03-19 15:58:59 afr Exp $
 *
 */


#define S1993

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1993(SISLCurve *c1,int *jstat)
#else
void s1993(c1,jstat)
     SISLCurve *c1;
     int   *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Check if a point-curve intersection in one dimention
*              is a simple case,
*              i.e. the intersection will result in one single point.
*
*
*
* INPUT      : c1    - SISLCurve in the intersection problem.
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
* WRITTEN BY : Arne Laksaa, SI, 89-06.
* Revised by : ALA and UJK 01.11.90, totale rewritten, same strategy
*              as in s1994 (surface case).
*********************************************************************
*/
{
  register int ki;

  int kk,kn;
  int kbez;
  double tmax;
  double tmin;
  double tdiff;
  double *scoef=SISL_NULL;
  double noice = (double)100.0 * REL_COMP_RES;   /* Noice killer */ 
  /* ----------------------------------------------------------- */
  
  /* Init to  simple case. */
  *jstat = 1;
  
  tmax = - HUGE;
  tmin =   HUGE;
  
  /* Get curve attributes. */
  kk  = c1->ik;
  kn  = c1->in;
  kbez = (kk == kn);
  
  /* Run through vertices to find
     intervall of first derivative. */
  
  for (ki = 1,scoef= c1->ecoef;ki < kn; ki++,scoef++ )
      {
	tdiff = *(scoef + 1) - *scoef;
	tmin = min(tmin,tdiff);
	tmax = max(tmax,tdiff);
      }
  
  if (fabs(tmin) < noice) tmin = DZERO; 
  if (fabs(tmax) < noice) tmax = DZERO; 
  
  
  /* Simple case when no genuin zero's of first derivative. */
  if (kbez && (tmin*tmax >=DZERO)) 
    *jstat = 1;
  else if (tmin*tmax > DZERO) 
    *jstat = 1;
  else if (tmin == tmax)
    *jstat = 1;
  else
    /* Not a simple case. */
    *jstat = 0;

}



