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
 * $Id: sh6isinsid.c,v 1.2 2001-03-19 16:06:03 afr Exp $
 *
 */


#define SH6ISINSIDE

#include "sislP.h"


#if defined(SISLNEEDPROTOTYPES)
void
sh6isinside (SISLObject * po1, SISLObject * po2, SISLIntpt *intpt, int *jstat)
#else
void
sh6isinside (po1, po2, intpt, jstat)
     SISLObject *po1;
     SISLObject *po2;
     SISLIntpt *intpt;
     int *jstat;
#endif

/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Test if an intersection point is inside the parametric space
*		of the objects.
*
*
* INPUT      : po1      - Pointer to the first object in the intersection.
*              po2      - Pointer to the second object in the intersection.
*              intpt    - Intersection point.
*
*
* OUTPUT     : jstat    - status messages
*                        < 0  : Error.
*			   0  : outside (UJK:at least one?) all objects
*			   1  : in the inner of all objects
*			   2  : on at least one edge/end of the objects
*			   3  : on one corner of one objects
*			   4  : on one corner of each objects
*                          5  : on one edge/end of both objects
*                               (not corner)
*
*
* METHOD     :
*
* CALLS      :
*
* REFERENCES :
*
* WRITTEN BY : Arne Laksaa, SI, Oslo, sep. -91.
* CORRECTED BY: UJK
* REVISED BY : ALA and VSK, 09.92.  Introduced output status *jstat = 5.
*********************************************************************
*/
{
 double minpar[4];
 double maxpar[4];
 int nrpar;
 int i;
 int is_inside = 1;
 int on_edge = 0;

   
 if (intpt != SISL_NULL)
 {
    /* Making the parametric boarders */
 
    if (po1->iobj == SISLSURFACE)
    {
       nrpar = 2;

       minpar[0] = po1->s1->et1[po1->s1->ik1-1];
       minpar[1] = po1->s1->et2[po1->s1->ik2-1];
       maxpar[0] = po1->s1->et1[po1->s1->in1];
       maxpar[1] = po1->s1->et2[po1->s1->in2];
    }
    else if (po1->iobj == SISLCURVE)
    {
       nrpar = 1;

       minpar[0] = po1->c1->et[po1->c1->ik-1];
       maxpar[0] = po1->c1->et[po1->c1->in];
    }
    else    /* SISLPOINT */
       nrpar = 0;

    if (po2->iobj == SISLSURFACE)
    {
       minpar[nrpar]   = po2->s1->et1[po2->s1->ik1-1];
       minpar[nrpar+1] = po2->s1->et2[po2->s1->ik2-1];
       maxpar[nrpar]   = po2->s1->et1[po2->s1->in1];
       maxpar[nrpar+1] = po2->s1->et2[po2->s1->in2];
       nrpar += 2;
    }
    else if (po2->iobj == SISLCURVE)
    {
       minpar[nrpar] = po2->c1->et[po2->c1->ik-1];
       maxpar[nrpar] = po2->c1->et[po2->c1->in];
       nrpar++;
    }
    /*else    SISLPOINT */
    
    if (nrpar != intpt->ipar) goto err106;

    /* Testing. */
    
    for (i = 0; (i < nrpar) && is_inside; i++)
    {
       if ((intpt->epar[i] <= maxpar[i] + REL_PAR_RES || 
	    DEQUAL(intpt->epar[i], maxpar[i])) &&
	    (intpt->epar[i] >= minpar[i] - REL_PAR_RES ||
	     DEQUAL(intpt->epar[i], minpar[i])))
	{
	   /* Int point is inside. */
	   
	   if (intpt->epar[i] >= maxpar[i] - REL_PAR_RES)
	      on_edge +=  (1 << (2*i));	/* On edge/end */
	   if (intpt->epar[i] <= minpar[i] + REL_PAR_RES)
	      on_edge +=  (1 << (2*i+1));	/* On edge/end */
	   
	}
	else  is_inside = 0;
    }
    
    if (is_inside)
    {
       (*jstat) = 1;
       if(on_edge)
       {
	  (*jstat) += 1;
	  if(on_edge > 1)
	  {
    	     if (po1->iobj == SISLSURFACE)
     	     {
		if ((on_edge & 1 || on_edge & 2) &&
		    (on_edge & 4 || on_edge & 8))
		   (*jstat) += 1;
	     }

	     if (po2->iobj == SISLSURFACE )
     	     {
		int ui = 2*(nrpar - 2);
		if ((on_edge & (1 << (ui)) || on_edge & (1 << (ui+1))) &&
		    (on_edge & (1 << (ui+2)) || on_edge & (1 << (ui+3))))
		   (*jstat) += 1;
	     }
	  }
       }
       
       /* Test if the intersection point lies at an edge in both
	  objects and is not registered as a corner point.       */
       
       if (*jstat == 2 && (on_edge & 15) && (on_edge & 240))
	  *jstat = 5;
    }
    else
       (*jstat) = 0;
 }
 else goto err108;

 
  /* Done. */

  goto out;

  /* Error in input. Conflicting dimensions.  */

err106:*jstat = -106;
  goto out;

  /* Error in input. No points */

err108:*jstat = -108;
  goto out;

out:
  return;
}
