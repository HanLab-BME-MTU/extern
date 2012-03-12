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
 * $Id: s1903.c,v 1.2 2001-03-19 15:58:55 afr Exp $
 *
 */


#define S1903

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1903 (double epar[], int in, int ik, int cuopen, double *eknots[], int *jstat)

#else
void
s1903 (epar, in, ik, cuopen, eknots, jstat)
     double epar[];
     int in;
     int ik;
     int cuopen;
     double *eknots[];
     int *jstat;

#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    :	To produce the knot vector of a B-spline basis satisfying the
*		interpolation requirements reflected in the epar array.
*
* INPUT      :	epar	- Array containing a parametrization of the
*			  interpolation conditions. Each interpolation
*			  condition has got a distinct parameter value
*			  expect form the cases where several conditions
*			  are conflicting. In that case a multiple parameter
*                         value indicates the need of a multiple knot. The
*                         parameter values are sorted in increasing order.
*                         The dimension of the array is 'in' if the curve is
*			  open and 'in+1' if it is closed.
*               in     	- Number of interpolation conditions.
*               ik      - Order of B-spline basis.
*		cuopen	- Open/closed curve.
*
* OUTPUT     : eknots - The produced knot vector. The dimension of
*                       the array is in+ik if the curve is open.
*			If the curve is closed the dimension of the array
*			is in+2*ik-1 if the curve is even and in+2*ik if it is
*			odd.
*              jstat  - status messages
*
* METHOD     :
*
*
* REFERENCES :
*
*
* CALLS      :
*
*
* WRITTEN BY :	Vibeke Skytt, SI, 91-03
* REVISED BY :	Trond Vidar Stensby, SI, 91-06
*
*********************************************************************
*/
{
  int kpos = 0;
  int ki;			/* Counter used to traverse the knot vector.	*/
  int kpar;			/* Counter used to traverse the parametrization
	                      	   array.					*/
  int kn;			/* The number of conditions in epar. (closed)   */
  int kk2;			/* Half the order.				*/
  int kstop;			/* Control variable of loop.            	*/
  int kmult;			/* Multiplisity of knot.                       	*/
  double tprev;			/* Value of previous knot.                     	*/
  double curr;			/* Value of current knot.                      	*/
  double tval1;			/* Start parameter value.			*/
  double tval2;			/* End parameter value.            		*/
  double tparint;		/* The parameter interval. (closed)		*/
  double tdum;			/* Help parameter used for parameter interval.	*/

  *jstat = 0;


  /* Check if curve is closed or open. */

  if (cuopen)
    {
      /* O P E N   C U R V E */

      *eknots = newarray (in +ik, DOUBLE);
      if (*eknots == SISL_NULL)
	goto err101;

      kk2 = ik / 2;
      tval1 = epar[0];
      tval2 = epar[in -1];

      /* Store a knot of multiplisity equal to the order at the start of the
	 curve The value of the knot is equal to the value of the start
	 parameter value. */

      for (ki = 0; ki < ik; ki++)
	(*eknots)[ki] = tval1;

      if (ik % 2 == 0)
	{
	  /* The order is even.
	     Place the internal knots at the parameter values.  */

	  for (kpar = kk2, kstop = in -kk2; kpar < kstop; kpar++, ki++)
	    (*eknots)[ki] = epar[kpar];
	}
      else
	{
	  /* The order is odd.
	     Place the internal knots between the parameter values.  */

	  for (kpar = kk2, kstop = in -kk2 - 1; kpar < kstop; kpar++, ki++)
	    (*eknots)[ki] = (double) 0.5 *(epar[kpar] + epar[kpar + 1]);
	}

      /* Store a knot of multiplisity equal to the order at the end of
	 the curve. The value of the knot is equal to the value of the
	 end parameter value.  */

      for (ki = 0; ki < ik; ki++)
	(*eknots)[in +ki] = tval2;
    }
  else
    {
      /* C L O S E D   C U R V E */

      *eknots = newarray (in +2 * ik, DOUBLE);
      if (*eknots == SISL_NULL)
	goto err101;

      kn = in +1;
      kk2 = ik / 2;
      kstop = in +2 * ik - 1;
      tparint = epar[in] -epar[0];

      if (ik % 2 == 0)
	{
	  /* The order of the B-spline curve is even.
	     Make the ik-1 first knots as a shift of the last knots.  */

	  for (ki = 0, kpar = in -ik + 1; ki < ik - 1; ki++, kpar++)
	    (*eknots)[ki] = epar[kpar] - tparint;

	  /* Make the knots corresponding to the data points. */

	  for (kpar = 0; kpar < kn; ki++, kpar++)
	    (*eknots)[ki] = epar[kpar];

	  /* Make the ik-1 last knots.  */

	  for (kpar = 1; ki < kstop; ki++, kpar++)
	    {
	      tdum = tparint;

	      /* We may risk that a double cyclic use of the parameter
	         values may result.          */

	      if (kpar > kn)
		{
		  tdum += tparint;
		  kpar -= in;
		}
	      (*eknots)[ki] = epar[kpar] + tdum;
	    }
	}
      else
	{
	  /* The order of the B-spline curve is odd.
	     Make the ik-1 first knots.             */

	  for (ki = 0, kpar = in -ik + 1; ki < ik - 1; ki++, kpar++)
	    (*eknots)[ki] = (double) 0.5 *(epar[kpar] + epar[kpar + 1]) - tparint;

	  /* Make the in next knots.  */

	  for (kpar = 0; kpar < in; ki++, kpar++)
	    (*eknots)[ki] = (double) 0.5 *(epar[kpar] + epar[kpar + 1]);

	  /* Make the ik remaining knots.  */

	  for (kpar = 0; ki < kstop; ki++, kpar++)
	    {
	      tdum = tparint;

	      /* We may risk that a double cyclic use of the parameter
	         values may result.        */

	      if (kpar > kn)
		{
		  tdum += tparint;
		  kpar -= in;
		}
	      (*eknots)[ki] = (double) 0.5 *(epar[kpar] + epar[kpar + 1]) + tdum;
	    }
	}
    }


  /* Check that the produced knots are in increasing order and that
     the multiplicity is not greater than ik.                       */

  if (cuopen)
    kstop = in +ik;

  for (ki = 1, tprev = (*eknots)[0], kmult = 0; ki < kstop; ki++, tprev = curr)
    {
       curr = (*eknots)[ki];
      kmult++;
      if (tprev > curr)
	goto err112;		/* Decreasing parameter value. */
      if (tprev < curr)
	kmult = 1;
      if (kmult > ik)
	goto err112;		/* Knot multiplisity greater than order. */
    }

  /* The knot vector is produced.  */

  goto out;


  /* Error in scratch allocation. */

err101:
  *jstat = -101;
  s6err ("s1903", *jstat, kpos);
  goto out;

  /* Error in the knot vector.  */

err112:
  *jstat = -112;
  s6err ("s1903", *jstat, kpos);
  goto out;

out:
  return;
}
