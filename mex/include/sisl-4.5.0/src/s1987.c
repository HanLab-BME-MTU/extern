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

#define S1987

#include "sislP.h"                                                 


#if defined(SISLNEEDPROTOTYPES)
void s1987(SISLSurf *ps, double aepsge, int *jgtpi, double **gaxis,
	   double *cang,int *jstat)
#else
void s1987(ps,aepsge,jgtpi,gaxis,cang,jstat)
     SISLSurf *ps;
     double aepsge;
     int   *jgtpi;
     double **gaxis;
     double *cang;
     int   *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Find the direction cone of a surface.
*
*
* INPUT      : ps        - Surface to treat.
*              aepsge    - Geometry tolerance.
*
* OUTPUT     : jgtpi     - To mark if the angle of the direction cone is
*                          greater than pi.
*                           0 - The direction cone of the surface
*                               is not greater than pi in any
*                               parameter direction.
*                           1 - The direction cone of the surface
*                               is greater than pi in the first
*                               parameter direction.
*                           2 - The direction cone of the surface is greater
*                               than pi in the second parameter direction.                          
*                          10 - The direction cone of a boundary curve of
*                               the surface is greater than pi in the first
*                               parameter direction.
*                          20 - The direction cone of a boundary curve of
*                               the surface is greater than pi in the second
*                               parameter direction.                      
*              gaxis     - Allocated array containing the coordinates of the
*                          center of the cone. It is only computed if
*                          *jgtpi = 0.
*              cang      - The angle from the center to the boundary of the
*                          cone. It is only computed if *jgtpi = 0.
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
* WRITTEN BY : Vibeke Skytt, SINTEF, 9403.
*
*********************************************************************
*/                                     
{
   int kstat = 0;        /* Local status variable.  */
   int kpos = 0;
   int kdim = ps->idim;
   
   /* Allocate scratch for the output array. */
   
   if ((*gaxis = newarray(kdim, DOUBLE)) == SISL_NULL) goto err101;
   
   /* Let s1990 compute the cone. */
   
   s1990(ps, aepsge, &kstat);
   if (kstat < 0) goto error;
   
   /* Copy the resulting cone to output parameters. */
   
   *jgtpi = ps->pdir->igtpi;
   *cang = ps->pdir->aang;
   memcopy(*gaxis, ps->pdir->ecoef, kdim, DOUBLE);
   
  /* Success ! */
  
  *jstat = 0;
  goto out;
  
  
  /* Error in space allocation.  */
  
  err101: 
     *jstat = -101;
  s6err("s1987",*jstat,kpos);
  goto out;
  
  /* Error in lower level routine. */
  
  error:
     *jstat = kstat;
  s6err("s1987",*jstat,kpos);
  goto out;
     
  out: 
    return;
}
