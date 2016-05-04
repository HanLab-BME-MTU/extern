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
 * $Id: s1240.c,v 1.2 2001-03-19 15:58:43 afr Exp $
 *
 */


#define S1240

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1240(SISLCurve *pcurve,double aepsge,double *clength,int *jstat)
#else
void s1240(pcurve,aepsge,clength,jstat)
     SISLCurve  *pcurve;
     double aepsge;
     double *clength;
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Calculate the length of a B-spline curve. The length
*              calculated will not deviate more than (aepsge/the
*              length calculated) from the real length of the curve.
*
*
*
* INPUT      : pcurve - Pointer to curve.
*              aepsge - Geometry resolution.
*
*
*
* OUTPUT     : clength - The length of the curve.
*              jstat   - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     :
*
*
* REFERENCES :
*
*-
* CALLS      : s6dist,s1251,make_cv_kreg.
*
* WRITTEN BY : Vibeke Skytt, SI, 88-11.
*
*********************************************************************
*/
{
  int kstat = 0;  /* Local status variable.                          */
  int kpos = 0;   /* Position of error.                              */
  int ki;         /* Counter.                                        */
  int kdim;       /* Dimension of the space in which the curve lies. */
  int kn;         /* Number of vertices of curve.                    */
  int kcalc;      /* Indicates if correct length of curve is found.  */
  double tlength; /* Length of curve.                                */
  double tprev;   /* Previous length of curve calculated.            */
  double teps;    /* Local tolerance.                                */
  double *s1;     /* Pointer used to traverse real array.            */
  SISLCurve *qc=SISL_NULL;  /* k-regular local curve.                     */
  
  if (pcurve->cuopen == SISL_CRV_PERIODIC)
    {
       /* Make curve k-regular. */
       
       make_cv_kreg(pcurve,&qc,&kstat);
       if (kstat < 0) goto error;
    }
  else qc = pcurve;
       
  /* Copy curve information to local parameters. */
  
  kdim = qc -> idim;
  kn   = qc -> in;
  
  /* Calculate length of control polygon.  */
  
  tlength = 0;
  for (ki=1,s1=qc->ecoef+kdim; ki<kn; ki++,s1+=kdim)
    tlength += s6dist(s1-kdim,s1,kdim);
  
  /* Set up local tolerance.  */
  
  teps = aepsge*100;
  
  kcalc = 0;
  while (kcalc == 0)
    {
      teps = teps/2.0;
      tprev = tlength;
      
      /* Compute length of curve.  */
      
      s1251(qc,teps,&tlength,&kstat);
      if (kstat < 0) goto error;
      
      /* Test if the error is within the tolerance. */
      
      if (fabs(tprev-tlength)/MAX(tprev,tlength) < aepsge) kcalc = 1;
      
    }
  
  /* Length of curve calculated. */
  
  *clength = tlength;
  *jstat = 0;
  goto out;
  
  /* Error in lower level routine. */
  
  error : *jstat = kstat;
  s6err("s1240",*jstat,kpos);
  goto out;
  
 out: 
    
    /* Free local curve.  */
    
    if (qc != SISL_NULL && qc != pcurve) freeCurve(qc);
    
    return;
}

