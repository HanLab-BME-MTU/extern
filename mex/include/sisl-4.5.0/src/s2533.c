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
 * $Id: s2533.c,v 1.2 2001-03-19 15:58:59 afr Exp $
 *
 */


#define S2533

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void s2533(double *et, int ik, int in, int multinc, int newik, int *newin,
	   double **newet, int *stat)
#else
void s2533(et, ik, in, multinc, newik, newin, newet, stat)
      double *et;
      int ik;
      int in;
      int multinc; 
      int newik;
      int *newin;
      double **newet;
      int *stat;
#endif
/*
*********************************************************************
*
* PURPOSE : To derive a new knot vector from an existing one, given a new
*           order and a new internal knot multiplicity.
*
*           We assume that the input knot vectors are k-regular, and that the
*           knot multiplicity is prechecked to avoid interior knot multiplicity
*           equal to or larger than the order.
*
*
*
* INPUT   : et        - The original knot vector.
*           ik        - The original order.
*           in        - The original number of coefficients.
*           multinc   - The multiplicity increment.
*                       In addition the multiplicity is increased by the order
*                       increase.
*           newik     - The new order.
*
*
*
* OUTPUT   : newin    - The new number of coefficients.
*            newet    - The new knot vector.
*            stat     - Status messages
*                       > 0      : Warning
*                       = 0      : Ok
*                       < 0      : Error
*
*
* METHOD   : 
*
*
* REFERENCES :
*
*-
* CALLS      : 
*
* WRITTEN BY :  Johannes Kaasa, SINTEF, Oslo, Norway.    Date: 1995-8
*
*********************************************************************
*/
{
   int ki, kj, kl;          /* Indices.                         */
   int add_knot;            /* Number of interior knots to add. */
   int knot_pos;            /* Number of knot positions.        */
   int *new_mult = SISL_NULL;    /* Array of new multiplicities.     */
   double *knot_par = SISL_NULL; /* Array of knot parameters.        */


   /* Check input. */

   if (et == SISL_NULL || multinc < 0 || newik < (multinc + 2)) 
      goto err150;

   /* Initiation and allocation of utility arrays. */

   if (in > ik)
   {
      /* Inner knots. */
      
      add_knot = (newik - ik) + multinc;
      
      if ((new_mult = newarray(in - ik, INT)) == SISL_NULL) goto err101;
      if ((knot_par = newarray(in - ik, DOUBLE)) == SISL_NULL) goto err101;
   }
   
   /* Examine the original knot multiplicity. */
   
   *newin = newik;
   knot_pos = 0;
   for (ki = ik, kl = 0; ki < in; kl++)
   {
      knot_par[kl] = et[ki];
      
      new_mult[kl] = add_knot + 1;
      kj = ki + 1;
      while (DEQUAL(et[kj], et[ki]))
      {
	 new_mult[kl]++;
	 kj++;
      }
      
      if (new_mult[kl] >= newik)
	 goto err150;
      
      *newin += new_mult[kl];
      knot_pos++;
      ki = kj;
   }
   
   /* Allocate the output array. */
   
   if ((*newet = newarray((*newin + newik), DOUBLE)) == SISL_NULL) goto err101;
   
   /* Fill in the new values. */
   
   for (kl = 0; kl < newik; kl++)
      (*newet)[kl] = et[ik - 1];
   
   for (ki = 0; ki < knot_pos; ki++)
   {
      for (kj = 0; kj < new_mult[ki]; kj++, kl++)
	 (*newet)[kl] = knot_par[ki];
   }
   
   for (ki = 0; ki < newik; ki++, kl++)
      (*newet)[kl] = et[in];
   
   goto out;
  
  
  
   /* ---------------------- ERROR EXITS ------------------------------- */

   /* Error in space allocation */
   
 err101: 
   *stat = -101;
   s6err("s2533", *stat, 0);
   goto out;

   /* Error in input. */
   
 err150:
   *stat = -150;
   s6err("s2533", *stat, 0);
   goto out;

   /* ---------------------- NORMAL EXIT ------------------------------- */

 out:
   if (new_mult != SISL_NULL) freearray(new_mult); 
   if (knot_par != SISL_NULL) freearray(knot_par); 
   
   return;
}
