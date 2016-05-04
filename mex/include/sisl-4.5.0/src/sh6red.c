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
 * $Id: sh6red.c,v 1.3 2005-02-28 09:04:50 afr Exp $
 *
 */


#define SH6RED

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
sh6red (SISLObject * po1, SISLObject * po2,
	SISLIntdat * pintdat, int *jstat)
#else
void
sh6red (po1, po2, pintdat, jstat)
     SISLObject *po1;
     SISLObject *po2;
     SISLIntdat *pintdat;
     int *jstat;
#endif

/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Reduse main int point to help point if they are not legal.
*
*
* INPUT      : po1      - Pointer to the first object in the intersection.
*              po2      - Pointer to the second object in the intersection.
*              pintdat  - Intersection data structure.
*
*
* OUTPUT     : jstat    - status messages
*                                > 0   : Warning.
*                                = 0   : Ok.
*                                < 0   : Error.
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
  int kstat, i, j;
  double tepsge = (double)10000.0*REL_COMP_RES;
  double weight = (double) 0.5;
  int changed;
  SISLIntpt *pcurr,*pstart,*plast;	/* to traverse list of points.     */
  int indstart,indlast,inddum;		/* Indexes used in lists           */
  int log_1, log_2;

  /* Remove all internal points in a list when along a
     constant parameter direction */
  
  if (((po1->iobj == SISLSURFACE && po2->iobj == SISLPOINT
        && po1->s1->idim == 1) ||
       (po2->iobj == SISLSURFACE && po1->iobj == SISLPOINT
        && po2->s1->idim == 1) ||
       (po1->iobj == SISLSURFACE && po2->iobj == SISLSURFACE
        && po1->s1->idim == 3)) &&
        pintdat != SISL_NULL)
     for (j = 0; j < pintdat->ipoint; j++)
     {
	
	pcurr = pintdat->vpoint[j];
	sh6isinside (po1, po2, pcurr, &kstat);
	if (kstat < 0)
	   goto error;
	
	/* VSK && ALA. 01.93. Do not remove points at corners. */
	if (kstat != 1 && kstat != 2) continue;
	
	sh6getnhbrs (pcurr, &pstart, &plast, &kstat);
	if (kstat < 0)
	   goto error;
	
	if (kstat == 0)
	{
	   /* Two neighbours, check */
	   sh6getlist (pcurr, pstart, &indstart, &inddum, &kstat);
	   if (kstat < 0)
	      goto error;		/* Error. */
	   if (kstat == 1)
	      goto errinconsist;	/* pcurr and pstart are not linked. */
	   
	   sh6getlist (pcurr, plast, &indlast, &inddum, &kstat);
	   if (kstat < 0)
	      goto error;		/* Error. */
	   if (kstat == 1)
	      goto errinconsist;	/* pcurr and plast are not linked. */
	   
	   log_1 = pcurr->curve_dir[indstart];
	   log_1 = log_1>>1;
	   log_1 &= 15;
	   log_2 = pcurr->curve_dir[indlast];
	   log_2 = log_2>>1;
	   log_2 &= 15;
	   	   
	   if (log_1 & log_2 )
	   {
	      sh6idkpt (&pintdat, &pcurr, 1, &kstat);
	      if (kstat < 0)
		 goto error;
	      /* Recursive nature : */
	      j = -1;
	   }
	   
	   
	}
     }
   
  
  if (pintdat != SISL_NULL)
    {
      /* Weight value in 3D sf vs sf case is one */
      if (pintdat->vpoint[0]->ipar == 4)
	weight = (double) 1.0;

      /* Reduce an illegal trim_curve to one point. */

      for (i = 0; i < pintdat->ipoint; i++)
	{
	  if (pintdat->vpoint[i]->iinter == SI_TRIM)
	    {
	      SISLIntpt **trim = SISL_NULL;
	      int no_trim = 0;
	      int no_alloc = 0;
	      sh6trimlist (pintdat->vpoint[i], &trim, &no_trim, &no_alloc);
	      for (j = 0; j < no_trim; j++)
		{
		  sh6isinside (po1, po2, trim[j], &kstat);
		  if (kstat < 0)
		    goto error;
		  if (kstat != 1)
		    break;
		}
	      if (j == no_trim)
		{
		  /* Internal trim area. */
		  for (j = 1; j < no_trim; j++)
		    {
		       /* sh6idunite (&pintdat, &trim[0], &trim[j], weight, &kstat);
			  */
		       /* VSK. 01.93. */
		       sh6idnewunite(po1, po2, &pintdat, &trim[0], &trim[j], 
				     weight, tepsge, &kstat);
		      if (kstat < 0)
			goto error;

		      /* We now need to correct the intpoint. */
		    }
		  trim[0]->iinter = SI_SING;
		}
	      if (trim)
		freearray (trim);
	    }
	}

      /* Reduse ilegal main points to help points. */
      do
	{
	  changed = 0;
	  for (i = 0; i < pintdat->ipoint; i++)
	    {
	      sh6isinside (po1, po2, pintdat->vpoint[i], &kstat);
	      if (kstat < 0)
		goto error;
	      if (kstat == 1)
		{
		  if (sh6ismain (pintdat->vpoint[i]) &&
		      sh6nmbmain (pintdat->vpoint[i], &kstat) == 1)
		    {
		      sh6tohelp (pintdat->vpoint[i], &kstat);
		      if (kstat < 0)
			goto error;
		      changed = 1;
		    }
		}
	    }
      } while (changed);
 
      /*UJK, 12.08.93 */
      /* Disconnect trim pts with 3 neighbours */
      do
	{
	   int ind_1,ind_2;
	   SISLIntpt *p_neighb[3];
	   int log_check[3];
	   changed = 0;
	   for (i = 0; i < pintdat->ipoint; i++)
	   {
	      pcurr = pintdat->vpoint[i];
	      sh6isinside (po1, po2, pcurr, &kstat);
	      if (kstat < 0)
		 goto error;
	      if (kstat &&
		  pcurr->iinter == SI_TRIM &&
		  sh6nmbmain (pcurr, &kstat) == 3)
	      {
		 for (ind_1=ind_2=0;ind_1<pcurr->no_of_curves;ind_1++)
		    if (pcurr->pnext[ind_1]->iinter == SI_TRIM)
		    {
		       sh6isinside (po1, po2, pcurr->pnext[ind_1], &kstat);
		       if (kstat < 0)
			  goto error;
		       if (kstat)
		       {
			  p_neighb[ind_2]  = pcurr->pnext[ind_1];
			  log_check[ind_2] = pcurr->curve_dir[ind_1];
			  log_check[ind_2] = log_check[ind_2]>>1;
			  log_check[ind_2] &= 15;
			  ind_2++;
		       }
		    }
		 
		 if (ind_2 == 3)
		 {
		    if (log_check[0] & log_check[1])
		       ind_2 = 2;
		    else if (log_check[0] & log_check[2])
		       ind_2 = 1;
		    else if (log_check[1] & log_check[2])
		       ind_2 = 0;
		    
		    if (ind_2 < 3)
		    {
		       changed = TRUE;
		       sh6disconnect(pcurr,p_neighb[ind_2],&kstat);
		       if (kstat < 0) goto error;
		       /* afr: Changed line below from an empty if-statement. */
		       sh6nmbmain (p_neighb[ind_2], &kstat);
		       if (kstat < 0) goto error;
		       sh6idkpt (&pintdat, &p_neighb[ind_2], 0, &kstat);
		       if (kstat < 0) goto error;
		    }
		 }
	      }
	   }
	} while (changed);
    }


  /* Reduction done. */

  (*jstat) = 0;
  goto out;

errinconsist:
  *jstat = -500;
  s6err ("sh6red", *jstat, 0);
  goto out;
  
error:(*jstat) = kstat;
  s6err ("sh6red", *jstat, 0);
  goto out;

out:
  return;
}

