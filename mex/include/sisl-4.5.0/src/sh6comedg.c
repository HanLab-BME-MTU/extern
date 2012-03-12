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
 * $Id: sh6comedg.c,v 1.2 2001-03-19 15:59:07 afr Exp $
 *
 */


#define SH6COMEDG

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
sh6comedg (SISLObject * po1, SISLObject * po2, SISLIntpt *pt1, SISLIntpt *pt2,
	   						int *jstat)
#else
void
sh6comedg (po1, po2, pt1, pt2, jstat)
     SISLObject *po1;
     SISLObject *po2;
     SISLIntpt *pt1;
     SISLIntpt *pt2;
     int *jstat;
#endif

/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Test if two intersection points is on a common edge and
*		connected along this edge on one
*		of the objects.
*
*
* INPUT      : po1      - Pointer to the first object in the intersection.
*              po2      - Pointer to the second object in the intersection.
*              pt1      - Intersection point.
*              pt2      - Intersection point.
*
*
* OUTPUT     : jstat    - status messages
*                        < 0  : Error.
*			   0  : Not on a common edge.
*			   1  : On a common edge in first object.
*			   2  : On a common edge in second object.
*			   3  : On a common edge in bouth objects.
*
*
* METHOD     :
*
* CALLS      :
*
* REFERENCES :
*
* WRITTEN BY : Arne Laksaa, SI, Oslo, sep. -91.
*********************************************************************
*/
{
 int kstat=0;
 double minpar[4];
 double maxpar[4];
 int nrpar, np1, np2;
 int common_edg, i, j;
 int is_inside = 1;
 int on_edge1 = 0;
 int on_edge2 = 0;
 int ind1, ind2;
 /* ---------------------------------------------------------------- */
 
 *jstat = 0;
   
 if (pt1 != SISL_NULL && pt2 != SISL_NULL)
 {
    /* Making the parametric boarders */
 
    if (po1->iobj == SISLSURFACE)
    {
       nrpar = 2;
       np1 = 4;

       minpar[0] = po1->s1->et1[po1->s1->ik1-1];
       minpar[1] = po1->s1->et2[po1->s1->ik2-1];
       maxpar[0] = po1->s1->et1[po1->s1->in1];
       maxpar[1] = po1->s1->et2[po1->s1->in2];
    }
    else if (po1->iobj == SISLCURVE)
    {
       nrpar = 1;
       np1 = 2;

       minpar[0] = po1->c1->et[po1->c1->ik-1];
       maxpar[0] = po1->c1->et[po1->c1->in];
    }
    else    /* SISLPOINT */
       np1 = nrpar = 0;

    if (po2->iobj == SISLSURFACE)
    {
       minpar[nrpar]   = po2->s1->et1[po2->s1->ik1-1];
       minpar[nrpar+1] = po2->s1->et2[po2->s1->ik2-1];
       maxpar[nrpar]   = po2->s1->et1[po2->s1->in1];
       maxpar[nrpar+1] = po2->s1->et2[po2->s1->in2];
       nrpar += 2;
       np2 = 4;
    }
    else if (po2->iobj == SISLCURVE)
    {
       minpar[nrpar] = po2->c1->et[po2->c1->ik-1];
       maxpar[nrpar] = po2->c1->et[po2->c1->in];
       nrpar++;
       np2 = 2;
    }
    else   np2 = 0;     /* SISLPOINT */

    /* Testing. */
    
    /* UJK, aug.92 */
    /* for (i = 0; i < nrpar || !is_inside; i++) */
    for (i = 0; i < nrpar && is_inside; i++)
      {
	 if (pt1->epar[i] <= maxpar[i] + REL_PAR_RES &&
	     pt1->epar[i] >= minpar[i] - REL_PAR_RES)
	   {
	      /* pt1 is inside. */
	      
	      if (pt1->epar[i] >= maxpar[i] - REL_PAR_RES)
		on_edge1 +=  (1 << (2*i));	/* On edge/end */
	      if (pt1->epar[i] <= minpar[i] + REL_PAR_RES)
		on_edge1 +=  (1 << (2*i+1));	/* On edge/end */
	      
	   }
	 else  is_inside = 0;
	 
	 if (pt2->epar[i] <= maxpar[i] + REL_PAR_RES &&
	     pt2->epar[i] >= minpar[i] - REL_PAR_RES)
	   {
	      /* pt2 is inside. */
	      
	      if (pt2->epar[i] >= maxpar[i] - REL_PAR_RES)
		on_edge2 +=  (1 << (2*i));	/* On edge/end */
	      if (pt2->epar[i] <= minpar[i] + REL_PAR_RES)
		on_edge2 +=  (1 << (2*i+1));	/* On edge/end */
	      
	   }
	 else  is_inside = 0;
      }
    
    common_edg = on_edge1 & on_edge2;
    (*jstat) = 0;
    
    if(is_inside && common_edg)
    {
       if (np1 > 0)
       {
	  j = (15>>(4-np1));
	  if (common_edg & j)
	  {
	     sh6getlist(pt1,pt2,&ind1,&ind2,&kstat);
             if (kstat < 0) goto err106;
             if (kstat == 0)
	     {
		if (common_edg & 3) i = 2;
		else i = 0;
		if (common_edg & (3<<2)) i+= 4;
		if (pt1->curve_dir[ind1] & i)  (*jstat) = 1;
	     }
	  }
       }
       if (np2 > 0)
       {
	  j = (15>>(4-np2));
	  j <<= np1;
	  if (common_edg & j)
	  {
	     sh6getlist(pt1,pt2,&ind1,&ind2,&kstat);
             if (kstat < 0) goto err106;
             if (kstat == 0)
	     {
		if (common_edg & (3<<np1)) i = 8;
		else i = 0;
		if (common_edg & (3<<(np1+2))) i+= 16;
		if (pt1->curve_dir[ind1] & i)  (*jstat) += 2;
	     }
	  }
       }
    }
    else
       (*jstat) = 0;
 }
 else goto err108;

 
  /* Done. */

  goto out;

  /* Error in input. Conflicting dimensions.  */

err106:*jstat = -106;
  s6err("sh6comedg",*jstat,0);
  goto out;

  /* Error in input. No points */

err108:*jstat = -108;
  s6err("sh6comedg",*jstat,0);
  goto out;

out:
  return;
}

