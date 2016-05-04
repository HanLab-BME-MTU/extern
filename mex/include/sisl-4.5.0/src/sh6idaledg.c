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
 * $Id: sh6idaledg.c,v 1.5 2005-02-28 09:04:50 afr Exp $
 *
 */


# define SH6ALLEDG
#include "sislP.h"
#if defined(SISLNEEDPROTOTYPES)
void
    sh6idalledg (SISLObject * pob1, SISLObject * pob2, SISLIntdat * pintdat,
	       SISLEdge * wedge[], int *jstat)
#else
       void
   sh6idalledg (pob1, pob2, pintdat, wedge, jstat)
     SISLObject *pob1;
     SISLObject *pob2;
     SISLIntdat *pintdat;
     SISLEdge *wedge[];
     int *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE     : Update edge structure.
*
*
*
* INPUT      : pob1    - Pointer to first object.
*              pob2    - Pointers to second object.
*              pintdat - intersection data.
*
*
*
* OUTPUT     : wedge[] - SISLEdge structures to update.
*              jstat   - status messages
*                               = 0     : OK!
*                               < 0     : error.
*
*
* METHOD     :
*
*
* REFERENCES :
*
*
* WRITTEN BY : Arne Laksaa, SI, 89-05.
* REVISED BY: Ulf J. Krystad, SI, 91-07.
* REVISED BY : Vibeke Skytt, SI, 92-09.
* Revised by : Paal Fugelli, SINTEF, Oslo, Norway, Nov.1994. Initialized
*              'kleft1' and 'kleft2' to avoid memory problems.
*********************************************************************
*/
{
  int kpos = 0;			/* Position of error.       */
  int kstat = 0;		/* Local error status.      */
  int ki,kj, kn;		/* Counters.                */
  int kj1, kn1;			/* Counters.                */
  int kndir;                    /* Number of par. dir.      */
  int kedg, kedg1;		/* Number of edges.         */
  int kpar;			/* Parameter number.        */
  int kleft1 = 0;               /* Index of knot.          */
  int ln[4];                    /* Number of vertices in each par. dir. */
  int lk[4];                    /* Order in each par. dir.  */
  double tpar;			/* Parameter value at edge. */
  double tparmain;              /* Parameter of main point. */
  double tparhelp;              /* Parameter of help point. */
  double *st[4];                /* Pointer to knot vector in each par. dir. */
  SISLObject *qob1, *qob2;	/* Help pointer to object.  */
  SISLObject *qob11, *qob21;	/* Help pointer to object.  */
  SISLPtedge *pte = SISL_NULL;
  SISLPtedge *prev = SISL_NULL;
  SISLPtedge *pte1 = SISL_NULL;
  SISLIntpt *pmain;
  int notfound;

  /* Set up information about the knot vector in each parameter direction. */

  if (pob1->iobj == SISLCURVE)
  {
     ln[0] = pob1->c1->in;
     lk[0] = pob1->c1->ik;
     st[0] = pob1->c1->et;
  }
  else if (pob1->iobj == SISLSURFACE)
  {
     ln[0] = pob1->s1->in1;
     lk[0] = pob1->s1->ik1;
     st[0] = pob1->s1->et1;
     ln[1] = pob1->s1->in2;
     lk[1] = pob1->s1->ik2;
     st[1] = pob1->s1->et2;
  }

  if (pob2->iobj == SISLCURVE)
  {
     ln[pob1->iobj] = pob2->c1->in;
     lk[pob1->iobj] = pob2->c1->ik;
     st[pob1->iobj] = pob2->c1->et;
  }
  else if (pob2->iobj == SISLSURFACE)
  {
     ln[pob1->iobj] = pob2->s1->in1;
     lk[pob1->iobj] = pob2->s1->ik1;
     st[pob1->iobj] = pob2->s1->et1;
     ln[pob1->iobj+1] = pob2->s1->in2;
     lk[pob1->iobj+1] = pob2->s1->ik2;
     st[pob1->iobj+1] = pob2->s1->et2;
  }

  for (kn = 0, qob1 = pob1, qob2 = pob2; kn < 2; kn++, qob1 = pob2, qob2 = pob1)
    {
      kedg = (qob1->iobj == SISLPOINT ? 0 : (qob1->iobj == SISLCURVE ? 2 : 4));

      if (kedg)
	wedge[kn]->ipoint = 0;

      for (kj = 0; kj < kedg; kj++)
	{
	  if (qob1->iobj == SISLCURVE)
	    {
	      tpar = (kj == 0 ? qob1->c1->et[qob1->c1->ik - 1] :
		      qob1->c1->et[qob1->c1->in]);
	      kpar = 1;
	    }
	  else if (kj == 0)
	    {
	      tpar = qob1->s1->et2[qob1->s1->ik2 - 1];
	      kpar = 2;
	    }
	  else if (kj == 1)
	    {
	      tpar = qob1->s1->et1[qob1->s1->in1];
	      kpar = 1;
	    }
	  else if (kj == 2)
	    {
	      tpar = qob1->s1->et2[qob1->s1->in2];
	      kpar = 2;
	    }
	  else
	    {
	      tpar = qob1->s1->et1[qob1->s1->ik1 - 1];
	      kpar = 1;
	    }

	  s6idedg ((kn == 0 ? qob1 : qob2), (kn == 0 ? qob2 : qob1),
		   kn + 1, kpar, tpar, pintdat,
		   &(wedge[kn]->prpt[kj]), &(wedge[kn]->ipoint), &kstat);
	  if (kstat < 0)
	    goto error;
	}
    }

  /* UJK newi; remove helppoints if main point is present */
  for (kn = 0, qob1 = pob1, qob2 = pob2; kn < 2; kn++, qob1 = pob2, qob2 = pob1)
    {
      kedg = (qob1->iobj == SISLPOINT ? 0 : (qob1->iobj == SISLCURVE ? 2 : 4));

      for (kj = 0; kj < kedg; kj++)
	{
	  for (pte = wedge[kn]->prpt[kj], prev = pte;
	       pte != SISL_NULL;)
	    {
	      notfound = TRUE;

	      if ((pmain = sh6getmain (pte->ppt)))
		{
		  /* pte is a help point and pmain is the
	             main point connected to it */

		   /* Check if the help point and main point lie
		      in different knot intervals in any parameter
		      direction. In that case, keep the help point.  */

		   for (kndir=pob1->iobj+pob2->iobj, ki=0;
		    ki<kndir; ki++)
		   {
		      tparmain = pmain->epar[ki];
		      tparhelp = pte->ppt->epar[ki];

		      /* Find position of parameter value in
			 to the knot vector.                  */

		      /* __________________________________ */
		      /* UJK, sept 93, this did not work
			 for left hand help pts. */
		      /*s1219(st[ki],lk[ki],ln[ki],&kleft1,tparmain,&kstat);
			 if (kstat < 0) goto error;

			 s1219(st[ki],lk[ki],ln[ki],&kleft2,tparhelp,&kstat);
			 if (kstat < 0) goto error;

			 if (kleft1 != kleft2) break;*/

		      s6fndintvl(st[ki],lk[ki],ln[ki],&kleft1,
				 tparmain,tparhelp,0,&kstat);
		      if (kstat < 0) goto error;

		      if (kstat) break;

		      /* UJK, sept 93, END */
		      /* __________________________________ */
		   }

		   if (ki == kndir)
		   {
		      /* Search for pmain */
		      for (kn1 = 0, qob11 = pob1, qob21 = pob2;
		       kn1 < 2 && notfound;
		       kn1++, qob11 = pob2, qob21 = pob1)
		      {
			 kedg1 = (qob11->iobj == SISLPOINT ?
				  0 : (qob11->iobj == SISLCURVE ? 2 : 4));

			 for (kj1 = 0; kj1 < kedg1 && notfound; kj1++)
			 {
			    for (pte1 = wedge[kn1]->prpt[kj1];
			     pte1 != SISL_NULL && notfound;
			     pte1 = pte1->pnext)
			       if (pte1->ppt == pmain)
				  notfound = FALSE;
			 }
		      }
		   }
		}

	      if (notfound == FALSE)
		{
		  /* Main point is present, remove help point. */
		  if (prev == pte)
		    {
		      wedge[kn]->prpt[kj] = pte->pnext;
		      freePtedge (pte);
		      pte = wedge[kn]->prpt[kj];
		      prev = pte;
		      wedge[kn]->ipoint--;
		    }
		  else
		    {
		      prev->pnext = pte->pnext;
		      freePtedge (pte);
		      pte = prev->pnext;
		      wedge[kn]->ipoint--;
		    }
		}
	      if (notfound == TRUE)
		{
		  sh6tomain(pte->ppt, &kstat);
		  prev = pte;
		  pte = pte->pnext;
		}

	    }
	}
    }


  *jstat = 0;

  goto out;

/* Error in lower level routine.      */

error:*jstat = kstat;
  s6err ("sh6idalledg", *jstat, kpos);
  goto out;

out:;
}

