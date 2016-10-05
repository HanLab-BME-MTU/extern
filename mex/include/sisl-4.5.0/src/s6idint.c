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
 * $Id: s6idint.c,v 1.2 2001-03-19 15:59:01 afr Exp $
 *
 */


#define S6IDINT

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s6idint(SISLObject *po1,SISLObject *po2,SISLIntdat *pintdat,SISLIntpt **rpt,int iob)
#else
void s6idint(po1,po2,pintdat,rpt,iob)
     SISLObject *po1;
     SISLObject *po2;
     SISLIntdat *pintdat;
     SISLIntpt  **rpt;
     int    iob;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : To find an internal intersection point in object iob
*              from pintdat.
*
*
*
* INPUT      : pintdat  - Pointer to intersection data.
*              po1      - Pointer to first object
*              po2      - Pointer to second object
*              iob      - Number of object to find internal 
*                         intersection poin in.
*
*
* OUTPUT     : rpt      - Pointer to an internal intersection point.
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
* WRITTEN BY : Arne Laksaa, 05.89.
*
*********************************************************************
*/                                     
{
  register int  ki,kj;
  int  kpar1,kpar2;
  double sstart1[2],send1[2];
  double sstart2[2],send2[2];
  
  
  /* Initiate to emty list. */
  
  *rpt = SISL_NULL;
  
  
  /* We have to be sure that we have an intdat structure. */
  
  if (pintdat == SISL_NULL)
    goto out;
  
  
  if (po1 == SISL_NULL || po1->iobj == SISLPOINT)
    kpar1 = 0;
  else if (po1->iobj == SISLCURVE)
    {
      kpar1 = 1;
      sstart1[0] = po1->c1->et[po1->c1->ik-1];
      send1[0] = po1->c1->et[po1->c1->in];
    }
  else if (po1->iobj == SISLSURFACE)
    {
      kpar1 = 2;
      sstart1[0] = po1->s1->et1[po1->s1->ik1-1];
      send1[0] = po1->s1->et1[po1->s1->in1];
      sstart1[1] = po1->s1->et2[po1->s1->ik2-1];
      send1[1] = po1->s1->et2[po1->s1->in2];
    }
  
  
  if (po2 == SISL_NULL || po2->iobj == SISLPOINT)
    kpar2 = 0;
  else if (po2->iobj == SISLCURVE)
    {
      kpar2 = 1;
      sstart2[0] = po2->c1->et[po2->c1->ik-1];
      send2[0] = po2->c1->et[po2->c1->in];
    }
  else if (po2->iobj == SISLSURFACE)
    {
      kpar2 = 2;
      sstart2[0] = po2->s1->et1[po2->s1->ik1-1];
      send2[0] = po2->s1->et1[po2->s1->in1];
      sstart2[1] = po2->s1->et2[po2->s1->ik2-1];
      send2[1] = po2->s1->et2[po2->s1->in2];
    }
  
  
  if (iob == 1 && kpar1 == 0)
    goto out;
  
  if (iob == 2 && kpar2 == 0)
    goto out;
  
  
  /* We have to go trough all intersection points to search for internal
     intersection points. */
  
  for (ki=pintdat->ipoint-1; ki>=0; ki--)
    {
      for (kj=0; kj<kpar1; kj++)
        if (sstart1[kj] > pintdat->vpoint[ki]->epar[kj]  ||
	    send1[kj] < pintdat->vpoint[ki]->epar[kj])
	  goto end;
      for (kj=0; kj<kpar2; kj++)
        if (sstart2[kj] > pintdat->vpoint[ki]->epar[kpar1+kj]  ||
	    send2[kj] < pintdat->vpoint[ki]->epar[kpar1+kj])
	  goto end;
      
      if (iob == 1)
        {
	  for (kj=0; kj<kpar1; kj++)
	    if (DEQUAL(sstart1[kj],pintdat->vpoint[ki]->epar[kj]) ||
	        DEQUAL(send1[kj],pintdat->vpoint[ki]->epar[kj]))
	      goto end;
        }
      else
        {
	  for (kj=0; kj<kpar2; kj++)
	    if (DEQUAL(sstart2[kj],pintdat->vpoint[ki]->epar[kpar1+kj]) ||
	        DEQUAL(send2[kj],pintdat->vpoint[ki]->epar[kpar1+kj]))
	      goto end;
        }
      
      
      (*rpt) = pintdat->vpoint[ki];
      goto out;
    end:;
    }
 out:;
}
