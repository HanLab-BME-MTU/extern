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
 * $Id: s1440.c,v 1.3 2001-03-19 15:58:49 afr Exp $
 *
 */


#define S1440

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1440(SISLSurf *ps1,SISLSurf **rs2,int *jstat)
#else
void s1440(ps1,rs2,jstat)
     SISLSurf *ps1;
     SISLSurf **rs2;
     int  *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Interchange the two parameter directions used in the
*              mathematical description of a surface and thereby
*              change the direction of the normal vector of the surface.
*
*
*
* INPUT      : ps1    - Pointer to the original surface.
*
*
*
* OUTPUT     : rs2    - Pointer to the surface with interchanged
*                       parameter directions.
*              jstat  - status messages
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
* CALLS      : s6chpar  - Change parameter directions of the vertices
*                         of the surface.
*              newSurf  - Create new surface.
*
* WRITTEN BY : Vibeke Skytt, SI, 88-11.
* REVISED BY : Johannes Kaasa, SI, 91-09 (Introduced NURBS).
* Revised by : Paal Fugelli, SINTEF, Oslo, Norway, Nov. 1994.  Added
*              handling of 'cuopen' flags.
*
*********************************************************************
*/
{
  int kpos = 0;          /* Position of error.                  */
  double *ssurf = SISL_NULL;  /* Pointer to vertices of new surface. */
  int kdim;              /* Local (rational) dimension.         */
  double *vert;          /* Pointer to vertices.                */

  /* Check for rational surface. */

  if (ps1->ikind == 2 || ps1->ikind == 4)
    {
      kdim = ps1->idim + 1;
      vert = ps1->rcoef;
    }
  else
    {
      kdim = ps1->idim;
      vert = ps1->ecoef;
    }

  /* Allocate scratch for vertices of new surface.  */

  ssurf = newarray(ps1->in1*ps1->in2*kdim,double);
  if (ssurf == SISL_NULL) goto err101;

  /* Change parameter directions of vertices.  */

  s6chpar(vert,ps1->in1,ps1->in2,kdim,ssurf);

  /* Create output surface.  */

  *rs2 = SISL_NULL;
  if ((*rs2 = newSurf(ps1->in2,ps1->in1,ps1->ik2,ps1->ik1,ps1->et2,
		      ps1->et1,ssurf,ps1->ikind,ps1->idim,1)) == SISL_NULL) goto err101;

  /* Set periodicity flag */

  (*rs2)->cuopen_1 = ps1->cuopen_2;
  (*rs2)->cuopen_2 = ps1->cuopen_1;

  /* Parameter directions changed.  */

  *jstat = 0;
  goto out;

  /* Error in space allocation.  */

 err101: *jstat = -101;
  s6err("s1440",*jstat,kpos);
  goto out;

 out:

  /* Free space occupied by local array.  */

  if (ssurf != SISL_NULL) freearray(ssurf);

  return;
}
