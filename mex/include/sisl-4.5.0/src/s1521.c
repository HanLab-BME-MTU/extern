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
 * $Id: s1521.c,v 1.2 2001-03-19 15:58:50 afr Exp $
 *
 */


#define S1521

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
SISLCurve* 
s1521(SISLCurve *pc,int *jstat)
#else
SISLCurve* s1521(pc,jstat)
     SISLCurve *pc;
     int   *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : To convert a (non-rational) B-spline curve
*              to rational form.
*
*
* INPUT      : pc    - Pointer to input curve.
*
*
* OUTPUT     : s1521  - Pointer to new output curve.
*              jstat  - status messages  
*                     = 1      : OK. The input curve is already rational.
*                                Output curve is simplya copy.
*                     = 0      : OK.  New curve made.
*                     < 0      : error
*
*
* METHOD     : The weights are set to 1 and a new curve is made.
*
*
* REFERENCES :
*
* CALLS      : newCurve
*
* WRITTEN BY : Michael Floater, SI, 91-09.
* Revised by : Christophe Rene Birkeland, SINTEF Oslo, May 1993.
*              rcoef is free'd only when ikind == 1 or 3
*
*********************************************************************
*/
{
  int in;                /* Number of vertices of curve.            */
  int ik;                /* Order of curve.                         */
  int ikind;             /* Kind of curve.                          */
  int idim;              /* Dimension of the curve's space.         */
  double *et;            /* Pointer to knot vector.                 */
  double *ecoef;         /* Pointer to vertices.                    */
  double *rcoef=SISL_NULL;    /* Pointer to homogeneous vertices.        */
  SISLCurve *ratpc=SISL_NULL; /* Output SISLCurve.                       */
  int i,j,J,jj;          /* loop variables                          */
 
  *jstat=0;
  in    = pc -> in;
  ik    = pc -> ik;
  ikind = pc -> ikind;
  idim  = pc -> idim;
  et    = pc -> et;
  ecoef = pc -> ecoef;
  
  if (ikind == 2 || ikind == 4)
  {
      *jstat=1;

      /* Pointer to the homogeneous vertices. */

      rcoef = pc -> rcoef;
  }
  else
  {
      /* Calculate the homogeneous vertices. */
    
      rcoef = newarray(in*(idim+1),DOUBLE);
      if(rcoef == SISL_NULL) goto err101;
    
      for(i=0, j=0, J=0; i<in; i++, j++)
      {
          for (jj=0; jj<idim; jj++, j++, J++) 
          {
              rcoef[j]=ecoef[J];
          }
    
          rcoef[j]=(double)1.0;
      } 
      ikind++;
  }

  /* Make new rational curve. */

  ratpc=newCurve(in,ik,et,rcoef,ikind,idim,1);
  if(ratpc == SISL_NULL) goto err101;


  /* New curve made. */
  
  goto out;
  
  /* Error in space allocation.  */
  
 err101: *jstat = -101;
  s6err("s1521",*jstat,0);
  goto out;
  
 out: 
  /* Free space occupied by local array.  */
  
  if(pc->ikind == 1 || pc->ikind == 3)
    if (rcoef != SISL_NULL) freearray(rcoef);
  
  return ratpc;
}
