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
 * $Id: s1733.c,v 1.2 1994-10-19 15:30:55 pfu Exp $
 *
 */


#define S1733

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1733(SISLSurf *ps,int icont1,int icont2,double *cstart1,double *cend1,
	   double *cstart2,double *cend2,double *gcoef,int *jstat)
#else
void s1733(ps,icont1,icont2,cstart1,cend1,cstart2,cend2,gcoef,jstat)
     SISLSurf   *ps;
     int    icont1;
     int    icont2;
     double *cstart1;
     double *cend1;
     double *cstart2;
     double *cend2;
     double *gcoef;
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : To pick out the next Bezier patch of a B-spline surface.
*              This function requere a B_spline surface that is the
*              result of s1731. This rutine do not check that the
*              surface is correct.
*
*
*
* INPUT      : ps         - B-spline surface to convert.
*              icont1     - The horisontal number of the Bezier patch to pick.
*              icont2     - The vertical number of the Bezier patch to pick.
*
*
*
* OUTPUT     : cstart1    - The horisontal start parameter value to
*                           the Bezier patch.
*              cend1      - The horisontal end parameter value to
*                           the Bezier patch.
*              cstart2    - The vertical start parameter value to
*                           the Bezier patch.
*              cend2      - The vertical end parameter value to
*                           the Bezier patch.
*              gcoef      - The vertices to the Bezier curve.
*              jstat      - status messages
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
* REVISED BY : Johannes Kaasa, SI, May 1992 (Introduced NURBS).
* REVISED BY : Christophe Birkeland, SI, July 92
* Revised by : Paal Fugelli, SINTEF, Oslo, Norway, Oct. 1994. Added lower
*              bound checks on 'icont1' and 'icont2' and corrected upper
*              bound checks from "<=" to "<" - caused memory problems.
*
**********************************************************************/
{
  int kpos=0;       /* Position of error.                  */
  int ki;           /* Control variable in loop.           */
  int kfi1,kla1;    /* Index to the first and last element
		       in the knot-vector with the start
		       value to the Bezier curve.          */
  int kfi2,kla2;    /* Index to the first and last element
		       in the knot-vector with the start
		       value to the Bezier curve.          */
  int kdim;         /* Dimension of the geometry.          */
  double *scoef;    /* Vertices.                           */

  *jstat = 0;

  /* Check if the surface is rational. */

  if (ps->ikind == 2 || ps->ikind == 4)
  {
     kdim = ps->idim + 1;
     scoef = ps->rcoef;
  }
  else
  {
     kdim = ps->idim;
     scoef = ps->ecoef;
  }

  /* Check that we have a Bezier patch to treat. */

  if ( icont1 >= 0  &&  icont2 >= 0  &&
       icont1 < ps->in1/ps->ik1  &&  icont2 < ps->in2/ps->ik2 )
    {
      /* Finding the first and last element in ps->et1 with the
	 start1 value. */

      kfi1 = icont1*ps->ik1;
      kla1 = kfi1 + ps->ik1;

      /* Updating the start1 and the end1 parameter value
	 to the patch. */

      *cstart1 = ps->et1[kfi1];
      *cend1 = ps->et1[kla1+1];

      /* Finding the first and last element in ps->et2 with the
	 start2 value. */

      kfi2 = icont2*ps->ik2;
      kla2 = kfi2 + ps->ik2;

      /* Updating the start2 and the end2 parameter value
	 to the patch. */

      *cstart2 = ps->et2[kfi2];
      *cend2 = ps->et2[kla2+1];

      /* Updating the vertices to the Bezier curve. */

      for (ki=0;ki < ps->ik2;ki++)
	memcopy(gcoef+ki*kdim*ps->ik1,
		&scoef[((kfi2+ki)*ps->in1 + kfi1)*kdim],
		kdim*ps->ik1,double);
    }
  else
    {
      /* Error, no patch to return. */

      *jstat = -152;
      s6err("s1733",*jstat,kpos);
    }
}
