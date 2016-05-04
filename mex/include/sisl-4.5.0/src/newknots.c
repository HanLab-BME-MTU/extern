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
 * $Id: newknots.c,v 1.2 2001-03-19 15:58:40 afr Exp $
 *
 */


#define NEWKNOTS

#include "sislP.h"


#if defined(SISLNEEDPROTOTYPES)
void
  newknots(double et[],int in,int ik,double epar[],int inpar,double aeps,
	   double **ginsert,int *jinsert,int *jstat)
#else
void newknots(et,in,ik,epar,inpar,aeps,ginsert,jinsert,jstat)
   double et[];
   int in;
   int ik;
   double epar[];
   int inpar;
   double aeps;
   double **ginsert;
   int *jinsert;
   int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : To set up an array of parameter values to insert into
*              a B-spline curve when the curve is to have knot 
*              multiplisity equal to the order at specified parameter
*              values.
*
*
*
* INPUT      : et      - Knot vector of the B-spline curve.
*              in      - Number of vertices of the curve.
*              ik      - Order of the curve.
*              epar    - Parameter values at which the refined curve
*                        is going to have an order-multiple knot.
*              inpar   - Number of parameter values in epar.
*              aeps    - The smallest legal distance between new and
*                        current knots.
*
*
*
* OUTPUT     : ginsert - Array containing new knots to the curve.
*              jinsert - Number of knots in ginsert.
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
* CALLS      : s1219  - Find position of parameter value in the knot
*                       vector of a curve.
*
* WRITTEN BY : Vibeke Skytt, SI, 11.90.
*
*********************************************************************
*/
{
   int kstat = 0;        /* Local status variable.               */
   int ki,kj;            /* Counters.                            */
   int knpar;            /* Number of new knots found currently. */
   int kleft = 0;        /* Position in knot vector.             */
   double tpar;          /* Value of current new knot.           */
   
   /* Allocate maximum needed scratch for output array.  */
   
   *jinsert = 0;
   if ((*ginsert = newarray(ik*inpar,DOUBLE)) == SISL_NULL) goto err101;
   
   /* Traverse the array of parameter values and compute how
      may knots is to be inserted at each value.  */
   
   for (ki=0; ki<inpar; ki++)
   {
      tpar = epar[ki];
      
      /* Find position in knot vector.  */
      
      s1219(et,ik,in,&kleft,tpar,&kstat);
      if (kstat < 0) goto error;
      
      /* Check if the parameter value is close to a knot.  */
      
      if (tpar-et[kleft] < aeps)
      {
	 tpar = et[kleft];
	 for (kj=1, knpar=ik-1; 
	  kj<=kleft && DEQUAL(et[kleft],et[kleft-kj]);
	  kj++, knpar--);
      }
      else if (et[kleft+1]-tpar < aeps)
      {
	 tpar = et[kleft+1];
	 for (kj=2, knpar=ik-1;
	  kleft+kj<in+ik && DEQUAL(et[kleft+1],et[kleft+kj]);
	  kj++, knpar--);
      }
      else knpar = ik;
      
      /* Register the new knots.  */
      
      for (kj=0; kj<knpar; kj++)
	 (*ginsert)[*jinsert + kj] = tpar;
      (*jinsert) += knpar;
   }
   
   /* Set correct size of output array.  */
   
   if (*jinsert != ik*inpar)
   {
      if ((*ginsert = increasearray(*ginsert,MAX(1,*jinsert),DOUBLE)) == SISL_NULL)
	 goto err101;
   }
   
   /* New knots found.  */
   
   *jstat = 0;
   goto out;
   
   /* Error in scratch allocation.  */
   
   err101 : *jstat = -101;
   goto out;
   
   /* Error in lower level routine.  */
    
   error : *jstat = kstat;
   goto out;
   
   out :
      return;
}
