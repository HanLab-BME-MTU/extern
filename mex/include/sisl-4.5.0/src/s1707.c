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
 * $Id: s1707.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S1707

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void s1707(SISLCurve *pc,int *jstat)
#else
void s1707(pc,jstat)
     SISLCurve *pc;
     int   *jstat;
#endif
/*
********************************************************************
*
*********************************************************************
*
* PURPOSE    : Check if a B-spline curve is correct.
*
* INPUT      : pc     - SISLCurve to treat.
*
* OUTPUT     : jstat     - status messages
*                        > 0      : warning (1&2 only used when 
*                                            cuopen=SISL_CRV_PERIODIC)
*                                      = 1: Cyclic but not full freedom.
*                                      = 2: Not cyclic.
*                                      = 8: Non-positive rational weights.
*                       = 0      : ok
*                       < 0      : error
*
* CALLS      :
*
* WRITTEN BY : Arne Laksaa, SI, 88-06.
* REVISED BY : Christophe Rene Birkeland, SI-SINTEF, May 1993.
*
**********************************************************************/
{

  int kpos=0;              /* Position of error. */
  int kstat=0;
  int step = 0;
  register double *s1,*s2; /* Pointers used in loop. */
  
  if (!pc) goto err150;

  if (pc->ik > pc->in) goto err111;
  
  if (pc->ik <= 0) goto err110;
  
  if (pc->in <= 0) goto err159;
  
  if (pc->idim <= 0) goto err102;
  
  if (pc->et[pc->in+pc->ik-1] <= *pc->et) goto err112;
  
  for (s1=pc->et,s2=pc->et+pc->in+pc->ik-1; s1<s2; s1++)
    if (s1[1] < *s1) goto err112;

  /* Check rational coefficients */
  if(pc->ikind == 2 || pc->ikind == 4)
    {
      step = pc->idim + 1;
      for (s1 = pc->rcoef + pc->idim, s2 = pc->rcoef + pc->in*step; 
	   s1 < s2; 
	   s1+= step)
	if (*s1 <= 0) goto war08;
    }

  /* Check if curve really is cyclic */
  if(pc->cuopen == SISL_CRV_PERIODIC)
    {
      test_cyclic_knots(pc->et,pc->in,pc->ik,&kstat);
      if (kstat < 0) goto error;
      if (kstat == 0) goto war02;
      if (kstat == 1) goto war01;
    }
      

  /* Updating output. No errors ! */
  
  *jstat = 0;
  goto out;
  
  /* Warning: Cuopen = SISL_CRV_PERIODIC, but knotvector does not give
   * full freedom. */
  
  war01:
    *jstat = 1;
    goto out;
  
  /* Warning: Cuopen = SISL_CRV_PERIODIC, but knotvector not cyclic. */
  
  war02:
    *jstat = 2;
    goto out;
  
  /* Warning: Non-positive rational coefficients. */
  
  war08:
    *jstat = 8;
    goto out;
  
  /* Dimension less than 1. */
  
  err102:
    *jstat = -102;
    s6err("s1707",*jstat,kpos);
    goto out;
  
  /* Error. Order less than 1. */
  
  err110:
    *jstat = -110;
    s6err("s1707",*jstat,kpos);
    goto out;
  
  /* Error. Order greater than number of vertices. */
   
  err111:
    *jstat = -111;
    s6err("s1707",*jstat,kpos);
    goto out;

  /* Error. Error in knotvector. */
  
  err112:
    *jstat = -112;
    s6err("s1707",*jstat,kpos);
    goto out;

  /* Error. Null pointer. */
  
  err150:
    *jstat = -150;
    s6err("s1707",*jstat,kpos);
    goto out;
  
  /* Error. Number of vertices less than 1. */
  
  err159:
    *jstat = -159;
    s6err("s1707",*jstat,kpos);
    goto out;
  
  /* Error in lower level routine */
      
  error:
    *jstat = kstat;
    s6err("s1707",*jstat,kpos);
    goto out;

  out: 
    return;
}
