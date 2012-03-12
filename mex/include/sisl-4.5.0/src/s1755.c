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
 * $Id: s1755.c,v 1.2 2005-02-28 09:04:48 afr Exp $
 *
 */


#define S1755

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1755 (double orknt[], int in, int ik, double extknt[],
       int *inh, int *jstat)

#else
void
s1755 (orknt, in, ik, extknt, inh, jstat)
     double orknt[];
     int in;
     int ik;
     double extknt[];
     int *inh;
     int *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : 	To produce a knot vector of one order higher than an other
*		knot vector reflecting the same continuity as the original
*		knot vector.
*
*
* INPUT      : 	orknt		- The original knot vector.
*				  Dimension (1:in+ik).
*		in		- The number of degrees of freedom in the
*				  B-spline given by the knot vector.
*		ik		- The order of the basis.
*
* OUTPUT     :  extknt		- The 'extended' knot vector.
*				  Dimension (1:inh+ik+1).
*		inh		- The number of degrees of freedom in the knot
*				  vector produced.
*               jstat           - Output status:
*                                  < 0: Error.
*                                  = 0: Ok.
*                                  > 0: Warning.
*
* METHOD     :	At orknt[ik-1], (ik+1) knots are inserted. For internal knot values,
*		the multiplicity 'M' is found, the continuity determined
*		'ik-M-1', and the new multiplicity 'M+1' deternined.
*		At orknt[in], (ik+1) knots are inserted.
*
* REFERENCES :
*
* CALLS      : s6err.
*
* WRITTEN BY : 	Christophe R. Birkeland, SI, 1991-07
* REWRITTEN BY :
* REVISED BY :
*
*********************************************************************
*/
{
  int ki;			/* Loop control parameters 		*/
  int kstart, kstop;
  int numb;
  int kpos = 0;			/* Position indicator for errors	*/
  double prev, par;		/* Parameters used to find consecutive
				   distinct knotvector values		*/
  double tstart, tstop;		/* tstart=orknt[ik-1], tstop=orknt[in]	*/

  *jstat = 0;


  /* Test if legal input. */

  if ((ik < 1) || (in <ik))
    goto err112;


  /* Test if input knot vector degenerate. */

  if (orknt[ik - 1] >= orknt[in])
    goto err112;

  kstop = in +ik;


  /* PRODUCTION OF KNOTS: First we fill in extra knots at each
     distinct knot value, then we remove the superfluous knots.	*/

  numb = 0;
  prev = orknt[0] - 1;
  for (ki = 0; ki < kstop; ki++)
    {
      par = orknt[ki];
      if (par < prev)
	goto err112;

      if (par != prev)
	{
	  /* New distinct knot value, insert additional knot. */

	  extknt[numb] = par;
	  numb++;
	}
      extknt[numb] = par;
      prev = par;
      numb++;
    }

  /* Remove superfluous knots at start. Find greatest start knot. */

  kstart = 0;
  tstart = orknt[ik - 1];
  while (extknt[kstart] <= tstart)
    kstart++;
  kstart--;


  /* Find smallest end knot 		*/

  kstop = numb - 1;
  tstop = orknt[in];
  while (extknt[kstop] >= tstop)
    kstop--;
  kstop++;


  /* The knots from kstart-ik up to
   * kstop+ik are the knots to be kept	*/

  *inh = kstop - kstart + ik;


  /* Copy the knots to be kept to the start of the knot array */

  memcopy (extknt, &extknt[kstart - ik], *inh + ik + 1, DOUBLE);
  goto out;


  /* Error in description of B-spline */

err112:
  *jstat = -112;
  s6err ("s1755", *jstat, kpos);
  goto out;

out:
  return;
}
