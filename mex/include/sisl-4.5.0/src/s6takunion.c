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
 * $Id: s6takunion.c,v 1.2 2001-03-19 15:59:02 afr Exp $
 *
 */


#define S6TAKEUNION

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void s6takeunion(double evec1[],int ielem1,double evec2[],int ielem2,
		 double **gunion,int *jnmbelem,int *jstat)
#else	 
void s6takeunion(evec1,ielem1,evec2,ielem2,gunion,jnmbelem,jstat)
     int ielem1,ielem2,*jnmbelem,*jstat;
     double evec1[],evec2[],**gunion;
#endif     
/*
*********************************************************************
*                                                                   
* PURPOSE    : Take the union between two sorted double vectors. Identical 
*              elements in the two vectors are only represented once. 
*              If several elements in one array are equal, all are
*              represented.
*
*
* INPUT      : evec1    - First vector.
*              ielem1   - Number of elements of evec1.
*              evec2    - Second vector.
*              ielem2   - Number of elements of evec2.
*                       
*
* OUTPUT     : gunion   - Union vector.
*              jnmbelem - Number of elements of gunion.
*              jstat    - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : 
*
* REFERENCES : 
*
* USE        : 
*
*-
* CALLS      : 
*
* WRITTEN BY : Vibeke Skytt, SI, 06.90.
*
*********************************************************************
*/
{
  int knelem;
  int knunion;
  double *sunion = SISL_NULL;
  double *s1,*s1stop;
  double *s2,*s2stop;
  
  /* Make local array to store the union of the vectors.  */

  knelem = ielem1 + ielem2;
  if ((sunion = newarray(knelem,DOUBLE)) == SISL_NULL) goto err101;
  
  /* Produce union vector. */

  for (s1=evec1,s1stop=s1+ielem1,s2=evec2,s2stop=s2+ielem2,knunion=0; 
       s1<s1stop && s2<s2stop;)
      {
	if (*s1 < *s2)
	  sunion[knunion++] = *s1++;
	else if (*s2 < *s1)
	  sunion[knunion++] = *s2++;
	else
	  {
	    sunion[knunion++] = *s1++;
	    s2++;
	  }
      }

  for (; s1<s1stop; s1++,knunion++)
    sunion[knunion] = *s1;

  for (; s2<s2stop; s2++,knunion++)
    sunion[knunion] = *s2;
  
  /* Allocate scratch for output union vector.  */

  *gunion = SISL_NULL;
  if ((*gunion = newarray(knunion,DOUBLE)) == SISL_NULL) goto err101;
  
  /* Copy union vector to output vector.  */

  memcopy(*gunion,sunion,knunion,DOUBLE);
  *jnmbelem = knunion;
  
  /* Union found.  */

  *jstat = 0;
  goto out;
  
  /* Error in scratch allocation.  */

  err101 :
    *jstat = -101;
  goto out;
  
  out :
    /* Free scratch occupied by local array.  */

    if (sunion != SISL_NULL) freearray(sunion);
  
  return;
}
