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
 * $Id: s1989.c,v 1.2 2001-03-19 15:58:58 afr Exp $
 *
 */


#define S1989

#include "sislP.h"                                                 


#if defined(SISLNEEDPROTOTYPES)
void s1989(SISLSurf *ps,double **emax,double **emin,int *jstat)
#else
void s1989(ps,emax,emin,jstat)
     SISLSurf *ps;
     double **emax;
     double **emin;
     int   *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Find the bounding box of the SISLSurf. NB. The geometric
*              bounding box is returned also in the rational case, that
*              is the box in homogenous coordinates is NOT computed.
*
*
* INPUT      : ps        - SISLSurface to treat.
*
* OUTPUT     : emin      - Array of dimension idim containing
*                          the minimum values of the bounding box,
*                          i.e. down-left corner of the box.
*              emax      - Array of dimension idim containing
*                          the maximum values of the bounding box,
*                          i.e. top-right corner of the box.
*              jstat     - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*                                                                     
*
*
* METHOD     :
*
*
* REFERENCES :
*
* CALLS      :
*
* WRITTEN BY : Christophe Rene Birkeland, SINTEF Oslo, July 1993.
*
*********************************************************************
*/                                     
{
  int i,j;                          /* Loop control variables    */
  int kpos = 0;                     /* Position of error.        */
  int bsdim;
  int len;
  int in = ps->in1 * ps->in2;
  double *coeff;
  double *minim=SISL_NULL;
  double *maxim=SISL_NULL;

  /* initialize variables */

  bsdim = ps->idim;
  coeff = ps->ecoef;
  len = bsdim;

  minim = newarray(bsdim, DOUBLE);
  maxim = newarray(bsdim, DOUBLE);
  if(minim == SISL_NULL || maxim == SISL_NULL) goto err101;

  for(j=0; j<bsdim; j++)
    {
      minim[j] = coeff[j];
      maxim[j] = coeff[j];
    }
  for(i=1, len=bsdim; i<in; i++, len+=bsdim)
    for(j=0; j<bsdim; j++)
      {
	minim[j] = MIN(minim[j], coeff[len+j]);
	maxim[j] = MAX(maxim[j], coeff[len+j]);
      }
  *emin = minim;
  *emax = maxim;      

  /* Success ! */
  
  *jstat = 0;
  goto out;
  
  
  /* Error in space allocation.  */
  
  err101: 
    *jstat = -101;
    s6err("s1989",*jstat,kpos);
    goto out;
  
  out: 
    return;
}
