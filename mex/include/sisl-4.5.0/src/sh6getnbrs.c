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
 * $Id: sh6getnbrs.c,v 1.2 2001-03-19 15:59:07 afr Exp $
 *
 */


#define SH6GETNHBRS

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
      sh6getnhbrs(SISLIntpt *pt,SISLIntpt **pt1,SISLIntpt **pt2,int *jstat)
#else
void sh6getnhbrs(pt,pt1,pt2,jstat)
   SISLIntpt *pt;
   SISLIntpt **pt1;
   SISLIntpt **pt2;
   int       *jstat;
#endif   
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Given the point pt, find both its neighbours if
*              they are unique.
*              If pt is help point, look for both type of neighbours.
*              If main point, look only for main points.
*
* INPUT      : pt       - SISLIntpt point.
*
* OUTPUT     : pt1      - One neighbour.
               pt2      - Second neighbour.
*              jstat    - Error flag.
*                         jstat =  0  => successful, 2 unique neighbours
*                         jstat =  1  => pt is end point, pt2 SISL_NULL
*                         jstat =  2  => pt is junction point, both SISL_NULL
*                         jstat =  3  => pt is isolated, both SISL_NULL
*                         jstat = -1  => error in data structure.
*                         jstat <  0  => error in lower level routine
*
*
* METHOD     : 
*
*
* REFERENCES :
*
* WRITTEN BY : Michael Floater, SI, Oslo, Norway. June 91.
*
*********************************************************************
*/
{
  int num;              /* count number of pointers    */
  int i;                /* Loop variable. */
  
   *pt1 = SISL_NULL;
   *pt2 = SISL_NULL;
   *jstat = 0;
  
  if(sh6ismain(pt))  /* pt is main point. */
  {
      num=0;

      for(i=0; i < pt->no_of_curves; i++)
      {
	  if(sh6ismain(pt->pnext[i]))
	  {
	      if(num == 0) *pt1 = pt->pnext[i];
	      else *pt2 = pt->pnext[i];
	      num++;
	  }
      }

      if(num == 0) *jstat = 3; /* pt is an isolated point. */
      else if(num == 1) *jstat = 1; /* pt is an end point. */
      else if(num > 2) /* pt is a junction point. */
      {
	  *pt1 = SISL_NULL;
	  *pt2 = SISL_NULL;
          *jstat = 2;
      }
  }
  else  /* pt is help point. */
  {
      num=pt->no_of_curves;

      if(num == 0) *jstat = 3; /* pt is an isolated point. */
      else
      {
          *pt1=pt->pnext[0];
          if(num == 1) *jstat = 1; /* pt is an end point. */
	  else
          {
              *pt2=pt->pnext[1];
              /* UJK; Oh, yeah ?, don't discriminate help points. */
	      /* if(num > 2) goto err1; Error in data structure. */
	      if (num > 2)
	      {
		 *pt1 = SISL_NULL;
		 *pt2 = SISL_NULL;
		 *jstat = 2;
	      }
          }
      }
  }
  
  goto out;
  

/* Error in data structure. */
  /*
err1: *jstat = -1;
      s6err("sh6getnhbrs",*jstat,0);
      goto out; */

   out:
      return;
}

