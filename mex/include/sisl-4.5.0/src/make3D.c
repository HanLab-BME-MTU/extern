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
 * $Id: make3D.c,v 1.3 2001-03-19 16:13:07 afr Exp $
 *
 */


#define MAKE3D

#include "sislP.h"


#if defined(SISLNEEDPROTOTYPES)
void
    make3D (SISLSurf *ps, SISLSurf **rsnew, int *jstat)
#else
void 
   make3D (ps, rsnew, jstat)
     SISLSurf *ps;
     SISLSurf **rsnew;
     int *jstat;
#endif
/*
********************************************************************
*
*********************************************************************
*
* PURPOSE    : To convert a 1D surface to a 3D representation.
*
*
*
* INPUT      : ps	- 1D Surface.
*
*
*
* OUTPUT     : rsnew	- The 3D surface.
*              jstat	- status messages
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : We want to change f(u,v) ----> (u,v,f(u,v)).
*              Using marsdens identity on u and v will do the job.
*              (In fact the u and v are translated to origin and scaled.)
*
*
* REFERENCES :
*
*-
*
* WRITTEN BY : Ulf J. Krystad, SI, 10.94.
* Revised by : 
*
**********************************************************************/
{
  
  int kk1,kk2,kn1,kn2;  /* Orders and numbers of vertices              */
  double *st1,*st2,*scoef; /* Knots and vertices of input surface      */
  double *s3coef=SISL_NULL;  /* 3-D coeff                                   */
  int kkm1,kkm2;        /* Orders minus 1                              */
  int kincre;           /* Number of doubles in first vertex direction */
  int ki,kj,kl,kstop;
  double tsum,*sp,*sq;
  //int kstat=0,kpos=0;
   
  if (!ps) goto errnull;
  if (ps->idim != 1) goto errdim;


  kk1   = ps -> ik1;
  kk2   = ps -> ik2;
  kn1   = ps -> in1;
  kn2   = ps -> in2;
  st1   = ps -> et1;
  st2   = ps -> et2;
  scoef = ps -> ecoef;

  /* Allocate array for 3-D representation of surface */
  
  if((s3coef = newarray(kn1*kn2*3,DOUBLE)) == SISL_NULL) goto err101;
  
  
  
  /* Make 3-D description of the surface */
  
  
  /* Make representation of coefficients from Marsdens identity for the
   * function f(t) = t, this will be used as the x-coordinate in the 3-D
   * representation */
  
  kkm1    = kk1 - 1;
  kincre  = 3*kn1;
  
  for (ki=0,kl=0,sp=s3coef ; ki<kn1 ; ki++,kl+=3,sp+=3)
    {
      tsum = (double)0.0;
      kstop = ki+kk1;
      for (kj=ki+1;kj<kstop;kj++)
        tsum +=st1[kj];
      
      tsum = tsum/kkm1;
      
      
      /* Copy x-coordinate to the other vertex rows */
      for (kj=0,sq=sp ; kj<kn2 ; kj++,sq+=kincre) *sq = tsum;
      
    }
  
  /* Make representation of coefficients from Marsdens identity for the
   * function f(t) = t, with the knot vector in second parameter direction
   * scaled to [0,tfak].  This will be used as the x-coordinate in the 3-D
   * representation */
  
  kkm2 = kk2 - 1;
  for (ki=0,sp=s3coef+1 ; ki< kn2 ; ki++)
    {
      tsum = (double)0.0;
      kstop = ki+kk2;
      for (kj=ki+1;kj<kstop;kj++)
        tsum +=st2[kj];
      
      tsum  = tsum/kkm2;
      
      /*  Copy to remaining y-coordinates in first vertex row */
      
      for (kj=0 ; kj<kn1 ; kj++,sp+=3) *sp = tsum;
      
    }
  
  /* Copy z-coordinates */
  
  for (kj=0,sp=s3coef+2,sq=scoef ; kj < kn2 ; kj++)
     for (ki=0 ; ki<kn1 ; ki++,sp+=3,sq++)
	*sp = *sq;
  
  /* Make 3-D surface */
  
  if(((*rsnew) = newSurf(kn1,kn2,kk1,kk2,st1,st2,s3coef,1,3,1)) == SISL_NULL) goto err101;
   
  goto out;

  /* ____________________________________________________________________ */
  /* Error in alloc */
err101:
   *jstat = -101;
  s6err ("make3D", *jstat, 0);
  goto out;

  /* Null pointer */
errnull:
   *jstat = -200;
  s6err ("make3D", *jstat, 0);
  goto out;

    /* Dimension NE 1 */
errdim:
   *jstat = -201;
  s6err ("make3D", *jstat, 0);
  goto out;

out:
   if (s3coef) freearray(s3coef);
     
}
