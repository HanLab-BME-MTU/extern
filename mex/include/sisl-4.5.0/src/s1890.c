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
 * $Id: s1890.c,v 1.3 2005-02-28 09:04:49 afr Exp $
 *
 */


#define S1890

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1890 (double oknots[], int oik, int oin, double *par[], int *der[],
       int *jstat)
#else
void
s1890 (oknots, oik, oin, par, der, jstat)
     double oknots[];
     int oik;
     int oin;
     double *par[];
     int *der[];
     int *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    :  To produce a set of parameter values and their corresponding
*		derivative indicatores.
*		The parameters will later be used in an interpolation problem.
*
* INPUT      :  oknots	- The original knot vector.
*		oik	- The order of the basis.
*		oin	- The number of degrees of freedom in the basis given
*			  by the knot vector.
*
* OUTPUT     :  par	- The computed parameter values.
*		der	- The derivate indicators. (all 0.0)
*		jstat	- Status variable:
*						> 0	: warning
*						= 0	: ok
*						< 0 	: error
*
* METHOD     :  First one parameter value is inserted at the start parameter
*		value, then ik+1 knots are added successively and divided by
*		ik+1 to produce internal parameter values. Then the end
*		parameter value is inserted. This procedure can result in
*		parameter values outside the legal range, these are then
*		corrected. All parameter values are assigned derivative
*		indicator 0.
*
* REFERENCES :  Fortran version:
*		Tor Dokken, SI, 1981-10
*
* CALLS      :  s6err.
*
* WRITTEN BY	: Christophe R. Birkeland & Trond Vidar Stensby, SI, 1991-06
*
*********************************************************************
*/
{
  int kpos = 0;
  int count1, count2;		/* Loop control variables     */
  int start, stop;
  int numb;			/* Number of wrong parameters */

  double sum;			/* Sum of knot values         */
  double pvl;			/* Single parameter value     */
  double delta;			/* Used for correcting wrong
				 * parameter values           */

  *jstat = 0;


  /* Test if legal input. */

  if (oik <= 1 || oin < oik)
    goto err112;


  /* Test if input knot vector degenerate. */

  if (oknots[oik - 1] >= oknots[oin])
    goto err112;


  /* Allocate arrays par and der. */

  *par = newarray (oin, DOUBLE);
  if (*par == SISL_NULL)
    goto err101;
  *der = new0array (oin, INT);
  if (*der == SISL_NULL)
    goto err101;


  /* P R O D U C E  P A R A M E T E R   V A L U E S.
   * First we produce parameter values by a simple algorithm.
   * The parameter values calculated in a wrong way are then corrected. */

  (*par)[0]       = oknots[oik - 1];
  (*par)[oin - 1] = oknots[oin];
  
  for (count1 = 2; count1 < oin; count1++)
    {
      stop = count1 + oik;
      sum = (double) 0.0;
      for (count2 = count1; count2 <= stop; count2++)
	sum += oknots[count2 - 1];
      (*par)[count1 - 1] = sum / (oik + 1);
    }

  /* Find second distinct knot value. */

  pvl = oknots[oik - 1];
  for (count1 = oik; oknots[count1] <= pvl; count1++) ;


  /* Find number of parameter values with wrong value at start of curve. */

  pvl = (oknots[oik - 1] + oknots[count1]) / (double)2.0;
  for (numb = 0, start = 1; (*par)[start] <= pvl; start++, numb++) ;

  if (numb > 0)
    {
      delta = (pvl - (*par)[0]) / (numb + 1);

      /* Fill inn missing parameter values. */

      pvl = (*par)[0] + delta;

      for (count1 = 1; count1 <= numb; count1++)
	{
	  (*par)[count1] = pvl;
	  pvl += delta;
	}
    }

  /* Find last but one distinct knot value. */

  pvl = oknots[oin];
  for (count1 = oin - 1; oknots[count1] >= pvl; count1--) ;


  /* Find end parameters in wrong interval. */

  pvl = (oknots[count1] + oknots[oin + 1]) / (double) 2.0;
  for (numb = 0, stop = oin - 2; (*par)[stop] >= pvl; stop--, numb++) ;

  if (numb > 0)
    {
      delta = ((*par)[oin - 1] - pvl) / (numb + 1);
      pvl = (*par)[oin - 1] - delta;
      for (count1 = 1; count1 <= numb; count1++)
	{
	  (*par)[oin - 1 - count1] = pvl;
	  pvl -= delta;
	}
    }

  /* Make derivative indicators */

  /* We used new0array which initializes all elements with zeroes 
   * and then this code is redundant.
   *
   * for (count1 = 0; count1 < oin; count1++)
   *  (*der)[count1] = 0;
   */
  /* Knots produced */

  goto out;


  /* Not enough memory. */

err101:
  *jstat = -101;
  s6err ("s1890", *jstat, kpos);
  goto out;

  /* Error in description of B-spline. */

err112:
  *jstat = -112;
  s6err ("s1890", *jstat, kpos);
  goto out;

out:
  return;
}
