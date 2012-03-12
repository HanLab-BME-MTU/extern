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
 *
 *
 */


#define S1516

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void s1516(double ep[],double epar[],int im,int idim,
           double **ev,int *jstat)
#else
void s1516(ep,epar,im,idim,ev,jstat)
     double ep[];
     double epar[];
     int    im;
     int    idim;
     double **ev;
     int    *jstat;
#endif
/*
************************************************************************
*
* Purpose:   To estimate the first derivative at each point in a sequence.
*
* Input:
*          ep     - Array containing the point in sequence
*                   (x,y,..,x,y,..), length idim * im.
*          epar   - Parametrization array. The array should be increasing
*                   in value.
*          im     - Number of point and derivatives
*          idim   - The dimension of the space the points and derivatives
*                   lie in
* Output:
*          ev     - Pointer to array containing the derivatives in sequence
*                   (x,y,..,x,y,..), length idim * im.
*          jstat  - Status variable
*                    < 0 - Error.
* Method:
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Michael Floater, SI 1993-10
*
*********************************************************************
*/
{
  int ki,kj;          /* Loop variables                              */
  int kk;             /* Polynomial order                            */
  int kpos=0;         /* Position of error                           */
  int kstat=0;        /* Status variable                             */
  double *gpar;
  int kcnsta;
  int kcnend;
  int iopen;
  int iorder;
  int ileft;
  double *ntype;
  SISLCurve *qc;
  int knbpar;
  double *evtemp;
  double cendpar;
  double *eder;




  /* Check input */

  if (idim < 1 || im < 2) goto err102;


  /* Allocate array for derivatives */

  evtemp    = newarray(idim*im,DOUBLE);
  if (evtemp == SISL_NULL) goto err101;

  ntype    = newarray(im,DOUBLE);
  if (ntype == SISL_NULL) goto err101;

  for(ki=0; ki<im; ki++)
  {
      ntype[ki] = 1.0;
  }

  eder    = newarray(2 * idim,DOUBLE);
  if (eder == SISL_NULL) goto err101;



  kcnsta = 1;
  kcnend = 1;
  iopen = 1;
  iorder = 4;

  s1358(ep, im, idim, ntype, epar, kcnsta, kcnend, iopen, iorder,
        epar[0],&cendpar, &qc, &gpar, &knbpar, &kstat);
     if(kstat < 0) goto error;

  for(ki=0,kk=0; ki<im; ki++,kk+=idim)
  {
      s1221(qc,1,epar[ki],&ileft,eder,&kstat);
      if(kstat < 0) goto error;

      for(kj=0; kj<idim; kj++)
      {
          evtemp[kk+kj] = eder[idim+kj];
      }
  }


  /* Calculation completed */

  /* Set result. */

  (*ev) = evtemp;

  *jstat = 0;
  goto out;



  /* Error in space allocation */

 err101: *jstat = -101;
  s6err("s1516",*jstat,kpos);
  goto out;


  /* Error in input. */

 err102: *jstat = -102;
  s6err("s1516",*jstat,kpos);
  goto out;

  /* Error in lower level routine. */

 error:  *jstat =kstat;
  s6err("s1516",*jstat,kpos);
  goto out;

 out:
  if (ntype != SISL_NULL) freearray(ntype);
  if (gpar != SISL_NULL) freearray(gpar);
  if (eder != SISL_NULL) freearray(eder);

  return;
}
