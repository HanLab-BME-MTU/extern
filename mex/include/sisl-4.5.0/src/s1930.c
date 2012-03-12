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
 * 
 *
 */


#define S1930

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1930 (int inbcrv, SISLCurve ** vpcrv, double **gknot2,
       double **gcoef2, int *jn2, int *jord2, int *jstat)
#else
void
s1930 (inbcrv, vpcrv, gknot2, gcoef2, jn2, jord2, jstat)
     int inbcrv;
     SISLCurve **vpcrv;
     double **gknot2;
     double **gcoef2;
     int *jn2;
     int *jord2;
     int *jstat;
#endif
/*
*********************************************************************
*
* PURPOSE    : Given a set of curves where some may be rational, 
*              put these on a common basis.
*              (Common knot-vector of length (jn2+jord2).
*              The vertices are recomputed according to this new basis.
*
* INPUT      : inbcrv - Number of curves in the curve-set.
*              vpcrv  - Array (length inbcrv) of pointers to the
*                       curves in the curve-set.
*
* OUTPUT     : gknot2 - Common knot-vector (new basis) for the curves.
*                       (jn2+jord2).
*              gcoef2 - The vertices of the inbcrv curves
*                       expressed in the new basis. Stored in sequence
*                       first curve, second curve,...
*                       Also the weight is stored. For non-rational curves
*                       the weight is one
*              jn2    - The no. of vertices in each of the inbcrv curves.
*              jord2  - The order of the new representation of the curves.
*              jstat  - Output status:
*                        < 0: Error.
*                        = 0: Ok.
*                        > 0: Warning.
*
* NOTE	     : This routine is an interface to s1931 using a rational format
*
*
* WRITTEN BY : Vibeke Skytt, SINTEF, 2009-08
*
*********************************************************************
*/
{
  int kstat = 0;
  int kpos = 0;
  SISLCurve **ratcurves=SISL_NULL;
  int ki, kj;
  int kdim = vpcrv[0]->idim;
  int rat;
  int kn;

  ratcurves = newarray(inbcrv,SISLCurve*);
  for(ki=0; ki<inbcrv; ki++)
  {
    /* Set pointers to the homogeneous coordinates or create array. */
    double *sc = SISL_NULL;
    rat = (vpcrv[ki]->ikind == 2 || vpcrv[ki]->ikind == 4);
    if (!rat)
      {
	kn = vpcrv[ki]->in;
	sc = newarray((kdim+1)*kn, double);
	if (!sc)
	  goto err101;

	for (kj=0; kj<kn; kj++)
	  {
	    memcopy(sc+ki*(kdim+1), vpcrv[ki]->ecoef+ki*kdim, kdim, double);
	    sc[ki*(kdim+1)+kdim] = 1.0;
	  }

      }
    ratcurves[ki] = newCurve(vpcrv[ki]->in,vpcrv[ki]->ik,vpcrv[ki]->et,
			     (rat) ? vpcrv[ki]->rcoef : sc, 1, kdim+1, 1);
    if (sc) freearray(sc);

    if (!ratcurves[ki])
      goto err101;
  }

  /* Put the curves into common basis. */

  s1931 (inbcrv, ratcurves, gknot2, gcoef2, jn2, jord2, &kstat);
  if (kstat < 0)
    goto error;

  *jstat = 0;
  goto out;

  err101:
    *jstat = -101;
    s6err("s1930",*jstat,kpos); 
    goto out;
  
  /* Error in lower level routine.  */

  error:
    *jstat = kstat;
    s6err ("s1930", *jstat, kpos);
    goto out;

  out:
    /* Release allocated curve pointer array and curves */
    for (ki=0; ki<inbcrv; ki++)
      freeCurve(ratcurves[ki]);
    freearray(ratcurves);
    return;
  }
