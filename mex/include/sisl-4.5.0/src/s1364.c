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
 * $Id: s1364.c,v 1.2 2001-03-19 15:58:47 afr Exp $
 *
 */


#define S1364

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1364(SISLCurve *pc,double aepsge,int *jstat)
#else
void s1364(pc,aepsge,jstat)
     SISLCurve  *pc;
     double aepsge;
     int    *jstat;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To decide if a B-spline curve is closed within a
*              tolerance
*
* INPUT      : pc     - The B-spline curve.   
*              aepsge - Geometric tolerance
*
* OUTPUT     : 
*              jstat  - status messages  
*                                         = 1      : SISLCurve closed
*                                         = 0      : SISLCurve open
*                                         < 0      : error
*
* METHOD     : 
*
*
* REFERENCES :
*
*-                                                 
* CALLS      : s1221, s6dist, s6err
*              
*
* WRITTEN BY : Tor Dokken, SI, Oslo, Norway. Nov 1988
*
*********************************************************************
*/
{
  int kn;             /* The number of B-splines, i.e., the dimension of
			 the spline space associated with the knot
			 vector.                                         */
  int kk;             /* The polynomial order of the curve.              */
  int kdim;           /* Dimension of space                              */
  int kleft=0;        /* Pointer to knots                                */
  int kder=0;         /* Derivatives to be calculated                    */
  int kstat;          /* Local status variable                           */
  int kpos=0;         /* Position of error                               */
  double *st;         /* Pointer to the first element of the knot vector
			 of the curve. The knot vector has [kn+kk]
			 elements.                                       */
  double sdum1[3];    /* Arrays for calculation of points                */ 
  double sdum2[3];    /* Arrays for calculation of points                */
  double *sder1 = SISL_NULL; /* Pointers to points                            */
  double *sder2 = SISL_NULL; /* Pointers to points                            */
  double tdist;       /* Distance between points                         */
  
  /* Check if curve is correct */
  
  s1707(pc,&kstat);
  if (kstat<0) goto error;
  
  
  /* Copy curve attributes to local parameters.  */
  
  kn = pc -> in;
  kk = pc -> ik;
  kdim = pc -> idim;
  st = pc -> et;
  
  if (kdim>3)
    {
      sder1 = newarray(kdim,DOUBLE);
      sder2 = newarray(kdim,DOUBLE);
    }
  else
    {
      sder1 = sdum1;
      sder2 = sdum2;
    }
  
  /* Calculate start point of curve */
  
  s1221(pc,kder,st[kk-1],&kleft,sder1,&kstat);
  if (kstat<0) goto error;
  
  /* Calculate end point of curve */
  
  s1221(pc,kder,st[kn],&kleft,sder2,&kstat);
  if (kstat<0) goto error;
  
  tdist = s6dist(sder1,sder2,kdim);
  
  if (tdist>aepsge)
    *jstat = 0;
  else
    *jstat = 1;
  goto out;
  
  /* Error in lower level function */
  
 error:  *jstat = kstat;
  s6err("s1364",*jstat,kpos);
  goto out;
 out:
  
  if (kdim>3)
    {
      if (sder1 != SISL_NULL) freearray(sder1);
      if (sder2 != SISL_NULL) freearray(sder2);
    }
  
  return;
}          
