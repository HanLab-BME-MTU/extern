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
 * $Id: s1439.c,v 1.2 1994-12-05 15:46:49 pfu Exp $
 *
 */


#define S1439

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void s1439(SISLSurf *ps1,double apar,int idirec,SISLCurve **rcurve,int *jstat)
#else
void s1439(ps1,apar,idirec,rcurve,jstat)
     SISLSurf   *ps1;
     double apar;
     int idirec;
     SISLCurve  **rcurve;
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Pick a curve along a constant parameter line in a NURBS
*              surface.
*              The constant parameter value used is apar and is in the
*              idirec parameter direction.
*              This routine replaces s1436() and s1437().
*
*
*
* INPUT      : ps1    - Surface.
*              apar   - Parameter value to use when picking out constant
*                       parameter curve.
*              idirec - Parameter direction in which to pick (must be 1 or 2)
*
*
* OUTPUT     : rcurve - Constant parameter curve.
*              jstat  - status messages
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     :
*
* REFERENCES :
*
* CALLS      : s1436, s1437 - These two routines do the job, which
*                             one is called depends on what parameter
*                             direction to pick from.
*
* WRITTEN BY : Christophe Rene Birkelan, SINTEF Oslo, July 1993.
*
*********************************************************************
*/
{
  int kpos = 0;      /* Position of error.                            */

  if(idirec == 1)
    {
      s1437(ps1, apar, rcurve, jstat);
      if(*jstat < 0) goto error;
    }
  else if(idirec == 2)
    {
      s1436(ps1, apar, rcurve, jstat);
      if(*jstat < 0) goto error;
    }
  else
    goto err115;

  /* Success !  Curve picked */

  goto out;


  /* Error in input parameter idirec.  */

  err115:
    *jstat = -115;
    s6err("s1439",*jstat,kpos);
    goto out;

  /* Error in lower level routine.  */

  error:
    s6err("s1439",*jstat,kpos);
    goto out;

  out:
    return;
}
