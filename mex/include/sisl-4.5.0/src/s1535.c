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

#define S1535

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void s1535(double points[],double der10[],double der01[],double der11[],
	   int im1,int im2,int idim,double par1[],
	   double par2[],int con1,int con2,int con3,int con4, 
	   int order1, int order2, int iopen1, int iopen2,
	   SISLSurf **rsurf,int *jstat)
#else
void s1535(points,der10, der01,der11,im1,im2,idim,par1, par2,
	   con1,con2,con3,con4,order1, order2,iopen1,iopen2,rsurf,jstat)
     double points[];
     double der10[];
     double der01[];
     double der11[];
     int im1;
     int im2;
     int idim;
     double par1[];
     double par2[];
     int con1;
     int con2;
     int con3;
     int con4;
     int order1;
     int order2;
     int iopen1;
     int iopen2;
     SISLSurf **rsurf;
     int *jstat;
#endif
/*
************************************************************************
*
* PURPOSE: To compute a B-spline tensor surface interpolating a set
*          of points.
*
* INPUT:
*          points - Array of dimension idim*im1*im2 containing
*                   the positions of the nodes (using the same ordering
*                   as ecoef in the SISLSurf structure).
*
*          der10  - Array of dimension idim*im1*im2 containing the first
*                   derivatives in the first parameter direction.
*
*          der01  - Array of dimension idim*im1*im2 containing the first
*                   derivatives in the second parameter direction.
*
*          der11  - Array of dimension idim*im1*im2 containing the cross
*                   derivatives (the twists).
*
*          im1    - The number of interpolation points in the
*                   first parameter direction.
*
*          im2    - The number of interpolation points in the
*                   second parameter direction.
*
*          idim   - Dimension of the space we are working in.
*
*          par1   - Parametrization in first parameter direction.
*                   For closed curves, one additional parameter value
*                   must be spesified. The last entry contains
*                   the parametrization of the repeted start point.
*                   (if the endpoint is equal to the startpoint of
*                   the interpolation the lenght of the array could
*                   be equal to im1 also in the closed case).
*
*          par2   - Parametrization in second parameter direction.
*                   For closed curves, one additional parameter value
*                   must be spesified. The last entry contains
*                   the parametrization of the repeted start point.
*                   (if the endpoint is equal to the startpoint of
*                   the interpolation the lenght of the array could
*                   be equal to im2 also in the closed case).
*
*
*                          ^ Second par. direction 
*                          |     
*                          |    (2.)
*                          |-----------|
*                          |           |
*                     (3.) |           | (4.) 
*                          |           |
*                          |           |
*                          |-----------|-> First par. direction
*                               (1.)
*
*          con1      - Additional condition along edge 1:
*                           = 0: No additional condition.
*                           = 1: Zero curvature.
*
*          con2      - Additional condition along edge 2:
*                           = 0: No additional condition.
*                           = 1: Zero curvature.
*
*          con3      - Additional condition along edge 3:
*                           = 0: No additional condition.
*                           = 1: Zero curvature.
*
*          con4      - Additional condition along edge 4:
*                           = 0: No additional condition.
*                           = 1: Zero curvature.
*
*          order1    - Order of surface in first parameter direction.
*
*          order2    - Order of surface in second parameter direction.
*
*          iopen1    - Open/close parameter in first parameter direction.
*                      =  1 : open surface.
*                      =  0 : closed, non-periodic surface.
*                      = -1 : periodic surface
*
*          iopen2    - Open/close parameter in second parameter direction.
*                      =  1 : open surface.
*                      =  0 : closed, non-periodic surface.
*                      = -1 : periodic surface
*
*
* Output:
*          rsurf - Pointer to the surf produced
*          jstat  - Status variable
*                    < 0 - Error.
*
* Method:
*     The interpolation is accomplished by using a one dimensional
*     routine for spline interpolation called several times. 
*     First, the datapoints
*     are considered to be idim*im1 dimentional and so on...
*
*
* REFERENCES :
*
* CALLS      : s1357, s6chpar.
*
* WRITTEN BY : Christophe Rene Birkeland, SINTEF, June 1993.
* CHANGED BY : Vibeke Skytt, SINTEF, 0394. Introduced iopen1 and iopen2.
*
*********************************************************************
*/                                                               
{
  int i, j, k, len;   /* Loop control parameter                      */
  int kpek1, kpek2, kpek3;
  int idimm1, newindim;
  int maxim;          /* Max (im1, im2)                              */
  int kstat=0;        /* Status variable                             */
  int kpos=0;         /* Position of error                           */
  int newin1, newin2; /* Number of vertices along par. dir. 1 & 2    */
  int numpt;              /* Needed in call to s1357                 */
  double start=0.;        /* Needed in call to s1357                 */
  double end;
  int *typept=SISL_NULL;       /* Array needed for call to s1357          */
  double *pointpar=SISL_NULL;  /* Array needed for call to s1357          */
  double *coeffpos=SISL_NULL;  /* Array needed for call to s1357          */
  double *coeffder=SISL_NULL;  /* Array needed for call to s1357          */
  double *coeffposder=SISL_NULL; /* Array needed for call to s1357          */
  double *newpoint=SISL_NULL;
  double *newder=SISL_NULL;
  SISLCurve *curve1a=SISL_NULL, *curve1b=SISL_NULL;
  SISLCurve *curve2=SISL_NULL;

  /* Allocate and initialize necessary arrays for call to s1357 */

  maxim = 2 * MAX(im1, im2);
  if((typept = newarray(maxim, INT))==SISL_NULL) goto err101;
  for(i=0; i<maxim; i+=2)
    {
      typept[i] = 1;
      typept[i+1] = 4;
    }
    
  idimm1 = idim*im1;
  len = 2*im2*idimm1;
  if((newpoint = newarray(len, DOUBLE))==SISL_NULL) goto err101;
  if((newder = newarray(len, DOUBLE))==SISL_NULL) goto err101;
  for(i=0, kpek1=0, kpek2=0, kpek3 = idimm1; 
      i<im2; 
      i++, kpek1+=(2*idimm1), kpek3+=(2*idimm1), kpek2+=idimm1)
    {
      for(k=0; k<idimm1; k++)
	{
	  newpoint[kpek1+k] = points[kpek2+k];
	  newpoint[kpek3+k] = der01[kpek2+k];
	  newder[kpek1+k] = der10[kpek2+k];
	  newder[kpek3+k] = der11[kpek2+k];
	}
    }

  /* INTERPOLATION in SECOND direction : position. */

  s1357(newpoint, im2*2, idimm1, typept, par2, con1, con2, iopen2, order2, 
	start, &end, &curve1a, &pointpar, &numpt, &kstat);
  if(kstat < 0) goto error;
  if(pointpar != SISL_NULL) 
    {
      freearray(pointpar);
      pointpar = SISL_NULL;
    }
  newin2 = curve1a->in;

  /* INTERPOLATION in SECOND direction : derivative. */

  s1357(newder, im2*2, idimm1, typept, par2, con1, con2, iopen2, order2, 
	start, &end, &curve1b, &pointpar, &numpt, &kstat);
  if(kstat < 0) goto error;
  if(pointpar != SISL_NULL) 
    {
      freearray(pointpar);
      pointpar = SISL_NULL;
    }
  if(curve1b->in != newin2) goto err116;

  /* Transpose results, store new coefficients in 
   * arrays coeffpos and coeffder */

  newindim = newin2 * idim;
  if((coeffpos = newarray(im1*newindim, DOUBLE)) == SISL_NULL)
    goto err101;
  if((coeffder = newarray(im1*newindim, DOUBLE)) == SISL_NULL)
    goto err101;
  s6chpar(curve1a->ecoef, im1, newin2, idim, coeffpos);
  s6chpar(curve1b->ecoef, im1, newin2, idim, coeffder);

  if((coeffposder = newarray(2*im1*newindim, DOUBLE)) == SISL_NULL)
    goto err101;
  for(j=0, kpek1=0, kpek2=newindim, kpek3=0; 
      j<im1; 
      j++, kpek1+=(2*newindim), kpek2+=(2*newindim), kpek3+=newindim)
    {
      for(k=0; k<newindim; k++)
	{
	  coeffposder[kpek1+k] = coeffpos[kpek3+k];
	  coeffposder[kpek2+k] = coeffder[kpek3+k];
	}
    }

  /* Interpolation in FIRST parameter direction */

  s1357(coeffposder, 2*im1, idim*newin2, typept, par1, con3, con4, iopen1, 
	order1, start, &end, &curve2, &pointpar, &numpt, &kstat);
  if(kstat < 0) goto error;
  if(pointpar != SISL_NULL) 
    {
      freearray(pointpar);
      pointpar = SISL_NULL;
    }
  newin1 = curve2->in;

  /* Transpose back coefficients */

  if((coeffposder=increasearray(coeffposder, idim*newin1*newin2, DOUBLE)) 
     == SISL_NULL)  goto err101;
  s6chpar(curve2->ecoef, newin2, newin1, idim, coeffposder);

  /* Create instance of surface */

  if (((*rsurf) = newSurf(newin1, newin2, order1, order2, curve2->et,
		     curve1a->et, coeffposder, 1, idim, 1)) == SISL_NULL) 
     goto err101;
  
  /* Set periodicity flag. */
  
  (*rsurf)->cuopen_1 = curve2->cuopen;
  (*rsurf)->cuopen_2 = curve1a->cuopen;

  /* Success */
  
  *jstat = 0;
  goto out;  
  
  /* Allocation error. */

  err101: 
    *jstat = -101;
    s6err("s1535",*jstat,kpos);
    goto out;
      
  /* Error. */

  err116: 
    *jstat = -116;
    s6err("s1535",*jstat,kpos);
    goto out;
      
  /* Error in lower level routine. */

  error:  *jstat =kstat;
    s6err("s1535",*jstat,kpos);
    goto out;
  
  out:
    /* Free arrays */
  
     if (typept != SISL_NULL) freearray(typept);
    if(newpoint != SISL_NULL) freearray(newpoint);
    if(newder != SISL_NULL) freearray(newder);
    if(coeffpos != SISL_NULL) freearray(coeffpos);
    if(coeffder != SISL_NULL) freearray(coeffder);
    if(coeffposder != SISL_NULL) freearray(coeffposder);

    /* Free local SISL-curve objects */
   
    if(curve1a != SISL_NULL) freeCurve(curve1a);
    if(curve1b != SISL_NULL) freeCurve(curve1b);
    if(curve2 != SISL_NULL) freeCurve(curve2);
  
    return;
}
