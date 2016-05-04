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
 * $Id: sh1925.c,v 1.2 2001-03-19 15:59:06 afr Exp $
 *
 */

#define SH1925

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
      sh1925(SISLCurve *pc1,SISLCurve *pc2,int idim,double ea[],int nstart[],
	     int nstop[],double emxerr[],double el2err[],int ilend,
	     int irend,int *jstat)
#else
void sh1925(pc1,pc2,idim,ea,nstart,nstop,emxerr,el2err,
	    ilend,irend,jstat)
   SISLCurve *pc1;
   SISLCurve *pc2;
   int idim;
   double ea[];
   int nstart[];
   int nstop[];
   double emxerr[];
   double el2err[];
   int ilend;
   int irend;
   int *jstat;
#endif     
/*
*********************************************************************
* 
* PURPOSE    : To finish off the calculation of a discrete least squares
*              spline approximation to the spline curve pc1, see the 
*              routine sh1926, and to calculate the discrete max-error
*              and L2-error of the approximation.
* 
* 
* INPUT      : pc1    - The spline curve to be approximated.
*              pc2    - The curve approximating pc1 on a subset of the
*                       knotvector.
*              idim   - The dimension of the geometry space.
*              ea     - Real array of dimension (pc1->in*pc1->ik) containing 
*                       the B-spline refinement matrix from the knot vector
*                       pc2->et to the knot vector pc1->et. This matrix has
*                       dimension pc1->in*pc2->inm byt subce at most
*                       pc1->ik entries are nonzero in each row, it can
*                       be stored in a pc1->in*pc1->ik array together
*                       with two integer arrays indicating the position
*                       of the first and last nonzero elements in each
*                       row.
*              nstart - Integer array of dimension (pc1->in) containing 
*                       pointers to the first nonzero element of each row 
*                       of the B-spline refinement matrix from pc2->et to
*                       pc1->et.
*              nstop  - Integer array of dimension (pc1->in) containing 
*                       pointers to the last nonzero element of each row 
*                       of the B-spline refinement matrix from pc2->et to
*                       pc1->et.
*              ilend  - The number of derivatives that have been kept fixed
*                       at the left end of the spline. This means that the 
*                       first ilend coefficients are not to be multiplied
*                       by the constants induced from the pc2->et knot
*                       vector.
*              irend  - The number of derivatives that have been kept fixed
*                       at the right end of the spline. This means that the 
*                       last irend coefficients are not to be multiplied
*                       by the constants induced from the pc2->et knot
*                       vector.
*
* 
* OUTPUT     : emxerr - Real array of dimension (idim) containing the 
*                       absolute value of the largest B-spline coefficient
*                       of f-g (see method below) in each component when
*                       f-g is expressed as a spline on the pc1->et knot
*                       vector.
*              el2err - Real array of dimension (idim) containing a
*                       weighted L2-norm of the B-spline coefficients of
*                       f-g in each component when f-g is expressed as a
*                       spline on the knot vector pc1->et.
*              jstat     - status messages 
*                          > 0 : warning 
*                          = 0 : ok 
*                          < 0 : error 
*             
* 
* METHOD     : First the B-spline coefficients pc2->ecoef of the spline
*              approximation g are multiplied by the square matrix dtau(-1/2),
*              i.e., pc2->ecoef[
*
*
* REFERENCES : Any book on general numerical analysis or numerical
*              linear algebra.
*              
*
* USE        :
*
*-
* CALLS      :   
*
* WRITTEN BY : Vibeke Skytt, SI, 05.92, on the basis of a routine
*              written by Tom Lyche and Knut Moerken, 12.85.
*
*********************************************************************
*/
{
   int ki,kj,kr;
   int kjh;
   int kk = pc1->ik;
   int km = pc1->in;
   int kn = pc2->in;
   int kj1,kj2;
   int kstop;
   double tkindv = (double)1.0/(double)kk;
   double tw;
   double thelp;
   double *st = pc1->et;
   double *sd = pc1->ecoef;
   double *stau = pc2->et;
   double *sc = pc2->ecoef;
   double *stemp = SISL_NULL;
   
   /* Allocate scratch for local array.  */
   
   if ((stemp = newarray(idim,DOUBLE)) == SISL_NULL) goto err101;
  
   /* Multiply the part of ec pointed to by kstart and kstop by the
      corresponding parts of the square matrix dtau(-1/2).   */
   
   for (kstop=kn-irend, ki=ilend; ki<kstop; ki++)
     {
	tw = sqrt((double)kk/(stau[ki+kk] - stau[ki]));
	for (kj=0; kj<idim; kj++)
	  sc[ki*idim+kj] *= tw;
     }
   
   /* Initiate arrays to zero.  */
   
   memzero(stemp,idim,DOUBLE);
   memzero(emxerr,idim,DOUBLE);
   memzero(el2err,idim,DOUBLE);
		  
   /* Express the approximating spline as a spline on et by multiplying
      ec by ea and then calculate the error in the spline approximation. */
   
   for (ki=0; ki<km; ki++)
     {
	memzero(stemp,idim,DOUBLE);
	
	/* Express the approximation as a spline on et.  */
	
	kj1 = nstart[ki];
	kj2 = nstop[ki];
	for (kjh=kk+kj1-kj2-1, kj=kj1; kj<=kj2; kjh++, kj++)
	  {
	     thelp = ea[ki*kk+kjh];
	     for (kr=0; kr<idim; kr++)
	       stemp[kr] += thelp*sc[kj*idim+kr];
	  }
	
	/* Calculate the maxerror and the weighted L2-error of the
	   approximation. */
	
	thelp = (st[ki+kk] - st[ki])*tkindv;
	for (kr=0; kr<idim; kr++)
	  {
	     stemp[kr] = fabs(stemp[kr] - sd[ki*idim+kr]);
	     el2err[kr] += thelp*stemp[kr]*stemp[kr];
	     if (stemp[kr] > emxerr[kr]) emxerr[kr] = stemp[kr];
	  }
     }
   for (kr=0; kr<idim; kr++)
     el2err[kr] = sqrt(el2err[kr]);
	
   
   *jstat = 0;
   goto out;
   
   /* Error in space allocation.  */
   
   err101: *jstat = -101;
   goto out;
   
   out:
      /* Free scratch used for local array.  */
      
      if (stemp != SISL_NULL) freearray(stemp);
	  
      return;
}
   
