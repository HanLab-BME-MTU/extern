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
 * $Id: s1880.c,v 1.2 2001-03-19 15:58:55 afr Exp $
 *
 */


#define S1880

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1880(int ipar1,int ipar2,int *jpt,SISLIntpt **vpoint,int *jlist,SISLIntlist **vlist,
	   int *jpar,double **gpar1,double **gpar2,int *jcrv,SISLIntcurve ***wcrv,int *jstat)
#else
void s1880(ipar1,ipar2,jpt,vpoint,jlist,vlist,jpar,gpar1,gpar2,jcrv,
           wcrv,jstat)
     int      ipar1;
     int      ipar2;
     int      *jpt;
     SISLIntpt    **vpoint;
     int      *jlist;
     SISLIntlist  **vlist;
     int      *jpar;
     double   **gpar1;
     double   **gpar2;
     int      *jcrv;
     SISLIntcurve ***wcrv;
     int      *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Transform intersection points and curves from internal
*              format in the recursive part of intersection routines
*              to output format.
*
*
*
* INPUT      : ipar1  - Number of parameter directions of first object.
*              ipar2  - Number of parameter directions of second object.
*              vpoint - Array containing intersection points on the
*                       internal format.
*              vlist  - Array representing intersection curves on the
*                       internal format.
*
* INPUT/OUTPUT : jpt    - Number of intersection point in the vpoint array.
*                jlist  - Number of lists representing intersection curves
*                         in the array vlist.
*
*
*
*
* OUTPUT     : jpar   - Number of single intersection points.
*              gpar1  - Parameter values of the single intersection points
*                       in the parameter area of the first object.
*              gpar2  - Parameter values of the single intersection points
*                       in the parameter area of the second object.
*              jcrv   - Number of intersection curves.
*              wcrv   - Array containing description of intersection curves.
*              jstat  - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     :
*
*
* REFERENCES :
*
*-
* CALLS      : newIntcurve - Create a new instance of Intcurve.
*              freeIntpt   - Free space occupied by intersection point.
*
* WRITTEN BY : Vibeke Skytt, SI, 88-05.
*
*********************************************************************
*/
{
  int kpos = 0;         /* Position of error.                          */
  int ki,kj,kk;         /* Counters.                                   */
  int kpoint;           /* Number of points in an intersection list.   */
  int klst;             /* Kind of intersection list. (See SISLIntlist).   */
  int ktype;            /* Kind of intersection curve. (See SISLIntcurve). */
  int kpt;              /* Used to find number of single intersection points.*/
  double *spar1,*spar2; /* Values of points belonging to an intersection
			   curve in the parameter area of the objects
			   involved in the intersection.               */
  double *stpar1,*stpar2,*stpar3; /* Pointers used to travers arrays 
				     containing parameter values.      */ 
  SISLIntcurve **ucrv;      /* Pointer used to traverse *wcrv array.    */
  SISLIntlist **ulst;       /* Pointer used to traverse vlist array.     */
  SISLIntpt *qpt;           /* Pointer to an intersection point.         */
  SISLIntpt **upt;          /* Pointer used to travers vpoint array.     */
  
  /* Initiate output arrays.  */
  
  *gpar1 = *gpar2 = SISL_NULL;  *wcrv = SISL_NULL;
  
  /* Allocate space for intersection curve array.  */
  
  *jcrv = *jlist;
  *wcrv = newarray(*jlist,SISLIntcurve*);
  if ((*jcrv) > 0 && *wcrv == SISL_NULL) goto err101;
  
  /* Transfer curve-information from vlist array to wcrv array. */
  
  ucrv = *wcrv;
  ulst = vlist;
  kpt = 0;
  for (ki=0; ki<(*jlist); ki++)
    {
      qpt = (*ulst) -> pfirst;
      
      /* Allocate space for arrays containing parameter vlaues of points 
	 in intersection curves.                                          */
      
      kpoint = (*ulst) -> inumb;
      if (kpoint == 0) goto err137;
      spar1 = newarray(ipar1*kpoint,double);
      spar2 = newarray(ipar2*kpoint,double);
      if ((ipar1 > 0 && spar1 == SISL_NULL) ||
	  (ipar2 > 0 && spar2 == SISL_NULL)) goto err101;
      
      /* Collect parameter values of the points in this intersection list
	 and distribute values to the objects in the intersection.         */
      
      kj = 0;
      stpar1 = spar1;   
      stpar2 = spar2;
      while (qpt != SISL_NULL && qpt -> ipar != -1)
	{
	  stpar3 = qpt -> epar;
	  for (kk=0; kk<ipar1; kk++) *(stpar1++) = *(stpar3++);
	  for (kk=0; kk<ipar2; kk++) *(stpar2++) = *(stpar3++);
	  qpt -> ipar = -1;
	  qpt = qpt -> pcurve;
	  kj++;
	}
      
      /* Find kind of intersection curve.  */
      
      klst = (*ulst) -> itype;
      if (klst == 0) ktype = 4;
      else if (klst == 1) ktype = 2;
      else if (klst == 2) ktype = 5;
      else if (klst == 3) ktype = 6;
      else if (klst == 4) ktype = 7;
      else if (klst == 5) ktype = 8;
      else goto err146;             
      
      /* Create new intersection curve.  */
      
      *ucrv = newIntcurve(kj,ipar1,ipar2,spar1,spar2,ktype);
      if (*ucrv == SISL_NULL) goto err101;
      
      kpt += kj;
      ucrv++;
      ulst++;
    }                  
  
  /* Find number of single intersection points.  */
  
  kpt = *jpt - kpt;
  
  /* Create arrays to keep parameter values of intersection points.  */
  
  *gpar1 = newarray(ipar1*kpt,double);
  *gpar2 = newarray(ipar2*kpt,double);
  if ((ipar1*kpt > 0 && *gpar1 == SISL_NULL) 
      || (ipar2*kpt > 0 && *gpar2 == SISL_NULL)) goto err101;
  
  /* Copy parameters of single intersection points into output-arrays. */
  
  kj = 0;
  upt = vpoint; 
  stpar1 = *gpar1;
  stpar2 = *gpar2;
  for (ki=0; ki<(*jpt); ki++)
    {
      qpt = *upt;     
      if (qpt != SISL_NULL)
	{  
	  if (qpt -> ipar != -1)
	    {
	      kj++;
	      stpar3 = qpt -> epar;
	      for (kk=0; kk<ipar1; kk++) *(stpar1++) = *(stpar3++);
	      for (kk=0; kk<ipar2; kk++) *(stpar2++) = *(stpar3++);
	    }     
	  
	  /* Free space occupied by current intersection point.  */
	  
	  freeIntpt(qpt);
	}
      
      upt++;
    }
  
  *jpar = kj;
  
  /* Adjust output arrays to correct length.  */
  
  if (kj*ipar1 > 0)
    {
      if ((*gpar1 = increasearray(*gpar1,kj*ipar1,double)) == SISL_NULL) goto err101;
    }
  else 
    {
      if (*gpar1 != SISL_NULL) freearray(*gpar1);
      *gpar1 = SISL_NULL;
    }
  if (kj*ipar2 > 0)
    {
      if ((*gpar2 = increasearray(*gpar2,kj*ipar2,double)) == SISL_NULL) goto err101;
    }
  else 
    {
      if (*gpar2 != SISL_NULL) freearray(*gpar2);
      *gpar2 = SISL_NULL;
    }
  
  /* Intersections copied to output format.  */
  
  *jpt = 0;
  *jstat = 0;
  goto out;
  
  /* Error in space allocation.  */
  
 err101: *jstat = -101;
  s6err("s1880",*jstat,kpos);
  goto out;
  
  /* Error in data-strucuture. Expected intersection point not found. */
  
 err137: *jstat = -137;
  s6err("s1880",*jstat,kpos);
  goto out;
  
  /* Unknown kind of intersection type.  */
  
 err146: *jstat = -146;
  s6err("s1880",*jstat,kpos);
  goto out;
  
 out: return;
}                                   

