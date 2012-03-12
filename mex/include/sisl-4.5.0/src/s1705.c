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
 * $Id: s1705.c,v 1.2 2005-02-28 09:04:48 afr Exp $
 *
 */


#define S1705

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1705(SISLCurve *pc,int *jstat)
#else
void s1705(pc,jstat)
     SISLCurve *pc;
     int   *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Remove unnessesary knots and vertices of a B-spline curve.
*
*
*
* INPUT      : pc     - SISLCurve to treat.
*
*
*
* OUTPUT     : jstat     - status messages
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : Traverse the knot vector and remove multiplisety higher than
*              the order of the curve, which give zero basis.
*
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Arne Laksaa, SI, 88-06.
* REVISED BY : Johannes Kaasa, SI, 92-04 (Introduced NURBS).
* REVISED BY : Christophe Birkeland, SINTEF, 93-05 (*jstat = 0 in start).
*
**********************************************************************/
{
  int kk=pc->ik;           /* Order of the input curve.                  */
  int kn=pc->in;           /* Number of the vertices in input curves.    */
  int kdim=pc->idim;       /* Dimensjon of the space in whice curve lies.*/
  int rdim = kdim + 1;     /* Rational dimension.                        */
  int knnew=0;             /* Number of vertices in the new curves.      */
  register int ki;                  /* Control variable in loop.         */
  register double *s1,*s2,*s3,*s4;  /* Pointers used in loop.            */
  register double *st=pc->et;       /* The first new knot-vector.        */
  register double *scoef=pc->ecoef; /* The first new vertice.            */
  double *rcoef = pc->rcoef; /* The rational vertices.                    */
  int kind = pc->ikind;    /* The type of curve, 2 and 4 are rational.   */

  *jstat = 0;
  
  /* s1 is used to traverse scoef, s2 is used to traverse st.
     We just remove unnecesary knots by kompressing the arraies. */
  
  s4 = rcoef;
  for (s1=scoef,s2=st,s3=s2+kn; s2<s3; s1+=kdim,s4+=rdim,s2++)
    if (s2[kk] > *s2)
      {
	/* Here we copies nessecary vertices to compress the vector. */
	
	for (ki=0; ki<kdim; ki++)  
	  scoef[knnew*kdim+ki] = s1[ki];
	
        /* Here we copy rational vertices. */
	
	if (kind == 2 || kind == 4)
	  {
            for (ki=0; ki<rdim; ki++)  
	      rcoef[knnew*rdim+ki] = s4[ki];
	  }
	
	/* Here we copies nessecary knots to compress the vector. */
	
	st[knnew] = *s2;
	
	/* Updating number of copied knots. */
	
	knnew++;
      }
  
  /* At last we have to copy the last kk knots. */
  
  for (ki=0; ki<kk; ki++) 
    st[knnew+ki] = s3[ki];
  
  /* If some knots are removed we have to update the size of the vectors.*/
  
  if (knnew == 0) 
    goto err111;
  else
    if (knnew < kn)
      {
	pc->in = knnew;
      }
  goto out;
  
  err111:
    *jstat = -111;
  
  out: 
    return;
}
