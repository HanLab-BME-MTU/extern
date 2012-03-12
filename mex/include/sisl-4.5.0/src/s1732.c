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
 * $Id: s1732.c,v 1.2 1994-10-19 14:55:31 pfu Exp $
 *
 */


#define S1732

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1732(SISLCurve *pc,int icont,double *cstart,double *cend,double *gcoef,int *jstat)
#else
void s1732(pc,icont,cstart,cend,gcoef,jstat)
     SISLCurve  *pc;
     int    icont;
     double *cstart;
     double *cend;
     double *gcoef;
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : To pick out the next Bezier curve of a B-spline curve.
*              This function requere a B_spline curve that is the
*              result of s1730. This rutine do not check that the
*              curve is correct.
*
*
*
* INPUT      : pc      - B-spline curve to convert.
*              icont   - The number of the Bezier curve to pick.
*
*
*
* OUTPUT     : cstart  - The start parameter value to the Bezier curve.
*              cend    - The end parameter value to the Bezier curve.
*              gcoef   - The vertices to the Bezier curve.
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
* CALLS      :
*
* WRITTEN BY : Arne Laksaa, SI, 88-11.
* REVISED BY : Johannes Kaasa, SI, May 1992 (Introduced NURBS)
* REVISED BY : Christophe Birkeland, July 1992 (Test line 97)
* Revised by : Christophe Rene Birkeland, SINTEF Oslo, May 1993.
*              *jstat = 0  in begining
* Revised by : Paal Fugelli, SINTEF, Oslo, Norway, Oct. 1994.
*              Added check on 'icont' lower bound.
*
*********************************************************************
*/
{
  int kpos=0;      /* Position of error.                  */
  int kfi,kla;     /* Index to the first and last element
		      in the knot-vector with the start
		      value to the Bezier curve.          */
  double *rcoef;   /* Potential rational vertices.        */
  int kdim;        /* Potential rational dimension.       */

  *jstat = 0;

  /* Check if this is a rational curve. */

  if (pc->ikind == 2 || pc->ikind ==4)
    {
       rcoef = pc->rcoef;
       kdim = pc->idim + 1;
    }
  else
    {
       rcoef = pc->ecoef;
       kdim = pc->idim;
    }

  /* Check that we have a Bezier curve to treat. */

  if ( icont >= 0  &&  icont < pc->in/pc->ik )
    {
      /* The first and last element in pc->et with the
	 start value. */

      kfi = icont*pc->ik;
      kla = kfi + pc->ik;

      /* Updating the start and the end parameter value
	 to the curve. */

      *cstart = pc->et[kfi];
      *cend = pc->et[kla+1];

      /* Updating the vertices to the Bezier curve. */

      memcopy(gcoef,&rcoef[kfi*kdim],kdim*pc->ik,double);
    }
  else
    {
      /* Error, no curve to return. */

      *jstat = -151;
      s6err("s1732",*jstat,kpos);
    }
}
