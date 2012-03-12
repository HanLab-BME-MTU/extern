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
 * $Id: s6idget.c,v 1.2 2001-03-19 15:59:01 afr Exp $
 *
 */


#define S6IDGET

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s6idget(SISLObject *po1,SISLObject *po2,int ipar,double apar,SISLIntdat *pintdat,
	     SISLIntdat **rintdat,int *jstat)
#else
void s6idget(po1,po2,ipar,apar,pintdat,rintdat,jstat)
     SISLObject *po1;
     SISLObject *po2;
     int    ipar;
     double apar;
     SISLIntdat *pintdat;
     SISLIntdat **rintdat;
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : To insert all points in one intdat with one less number
*              of parameters into rintdat. New copies is made with
*              the missing parameter. All listes is also to
*              be uppdated.
*
*
*
* INPUT      : po1     - First object in the intersection.
*              po2     - Second object in the intersection.
*              ipar    - Number of the parameter that is missing in rintdat.
*              apar    - Parameter value of the missing parameter.
*              pintdat - Pointer to intersection data on the mother problem.
*
*
* OUTPUT     : rintdat  - Pointer to a pointer to intersection data.
*              jstat    - status messages  
*                               = 0      : Inserting done.
*                               < 0      : error
*
*
* METHOD     : 
*
*
* REFERENCES :
*
*-
* CALLS      : s6err      - Gives error message.
*              s6idnpt    - Insert a new intpt structure.
*              copyIntpt  - Copy an intpt structure.
*              newIntdat  - Create new intdat structure.
*
* WRITTEN BY : Arne Laksaa, 06.89.
* Revised by : Christophe Rene Birkeland, SINTEF Oslo, May 1993.
*              Error message corrected
*
*********************************************************************
*/                                     
{
  int kstat;                    /* Local status variable.                 */
  int kpos=0;                   /* Position of error.                     */
  int ki,kj,kn;
  double tstart[4];
  double tend[4];
  double spar[4];               /* Storing uppdated parametervalues.      */
  SISLIntpt *qpt;
  
  
  if (pintdat == SISL_NULL)
    goto out;
  
  if (po1->iobj == SISLCURVE)
    {
      tstart[0] = po1->c1->et[po1->c1->ik - 1];
      tend[0]   = po1->c1->et[po1->c1->in];
    }
  else if (po1->iobj == SISLSURFACE)
    {
      tstart[0] = po1->s1->et1[po1->s1->ik1 - 1];
      tend[0]   = po1->s1->et1[po1->s1->in1];
      tstart[1] = po1->s1->et2[po1->s1->ik2 - 1];
      tend[1]   = po1->s1->et2[po1->s1->in2];
    }
  
  if (po2->iobj == SISLCURVE)
    {
      tstart[po1->iobj]   = po2->c1->et[po2->c1->ik - 1];
      tend[po1->iobj]     = po2->c1->et[po2->c1->in];
    }
  else if (po2->iobj == SISLSURFACE)
    {
      tstart[po1->iobj]   = po2->s1->et1[po2->s1->ik1 - 1];
      tend[po1->iobj]     = po2->s1->et1[po2->s1->in1];
      tstart[po1->iobj+1] = po2->s1->et2[po2->s1->ik2 - 1];
      tend[po1->iobj+1]   = po2->s1->et2[po2->s1->in2];
    }
  
  tstart[ipar] = tend[ipar] = apar;
  
  
  /* Uppdate the array. */
  
  for (ki=0; ki<pintdat->ipoint; ki++)
    {
      for (kj=0; kj<pintdat->vpoint[ki]->ipar; kj++)
	if((DNEQUAL(pintdat->vpoint[ki]->epar[kj],tstart[kj]) &&
	    pintdat->vpoint[ki]->epar[kj] < tstart[kj]) ||
	   (DNEQUAL(pintdat->vpoint[ki]->epar[kj],tend[kj]) &&
	    pintdat->vpoint[ki]->epar[kj] > tend[kj]))
	  break;
      
      if (kj == pintdat->vpoint[ki]->ipar)
	{
	  for(kn=0; kn<ipar; kn++) spar[kn] = pintdat->vpoint[ki]->epar[kn];
	  for(; kn<pintdat->vpoint[ki]->ipar-1; kn++)
	    spar[kn] = pintdat->vpoint[ki]->epar[kn+1];
	  
	  /* Than we can insert all new intersection points in rintdat. */
	  
	  qpt = newIntpt(pintdat->vpoint[ki]->ipar-1,spar,DZERO);
	  if (qpt == SISL_NULL) goto err101;
	  
	  s6idnpt(rintdat,&qpt,1,&kstat);
	  if (kstat < 0) goto error;
	}
    }
  
  *jstat = 0;
  goto out;
  
  /* Error in space allocation.  */

  err101: 
    *jstat = -101;
    s6err("s6idget",*jstat,kpos);
    goto out;

  /* Error in sub function.  */

  error: 
    *jstat = kstat;
    s6err("s6idget",*jstat,kpos);
    goto out;

  out: ;
}
