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
 * $Id: s1531.c,v 1.3 2005-02-28 09:04:48 afr Exp $
 *
 */


#define S1531

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void s1531(double ea[],int idim,int in1,int in2,double **eb,int *jstat)
#else
void s1531(ea,idim,in1,in2,eb,jstat)
     double ea[];
     int idim;
     int in1;
     int in2;
     double **eb;
     int    *jstat;
#endif
/*
************************************************************************
*
* PURPOSE: To compute the transpose in the last two indices, of the
*          matrix given by ea and to output it as eb.
*
* INPUT:
*          ea     - Array of dimension idim*in1*in2 containing the
*                   matrix to be transformed.
*          idim   - The length of the first index of ea.
*          in1    - The length of the second index of ea.
*          in2    - The length of the third index of ea.
*
* OUTPUT:
*          eb     - Pointer to the output array containing the
*                   transposed matrix (dimension idim*in1*in2).
*          jstat  - Status variable
*                    < 0 - Memory allocation error.
* METHOD:
*
* REFERENCES :
*
* CALLS:
*
* WRITTEN BY: Michael Floater, SI, June 1992.
*
*********************************************************************
*/
{
  int i,j,jbase;       /* Loop variable                             */
  int ki,kj,kk;        /* Loop variable                             */
  int idiff;           /*                                           */
  int kpos=0;          /* Position of error                         */
  double *mat=SISL_NULL;    /* Temporary output matrix                   */


  mat = newarray(idim*in1*in2, DOUBLE);
  if(mat == SISL_NULL) goto err101;

  i = 0;
  jbase = 0;
  idiff = (in1 - 1) * idim;

  for(ki=0; ki<in1; ki++,jbase+=idim)
  {
      for(kj=0,j=jbase; kj<in2; kj++,j+=idiff)
      {
          for(kk=0; kk<idim; kk++,i++,j++)
          {
	      mat[i] = ea[j];
          }
      }
  }

  (*eb) = mat;

  /* Calculation completed */

  *jstat = 0;
  goto out;


  /* Error in space allocation */
 err101: *jstat = -101;
  s6err("s1531",*jstat,kpos);
  goto out;


 out:

  return;
}
