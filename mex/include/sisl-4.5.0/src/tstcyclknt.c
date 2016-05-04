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
 * $Id: tstcyclknt.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define TEST_CYCLIC_KNOTS

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
  test_cyclic_knots(double et[],int in,int ik,int *jstat)
#else
void test_cyclic_knots(et,in,ik,jstat)
     double et[];
     int    in;
     int    ik;
     int    *jstat;
#endif
     /*
      *********************************************************************
      *                                                                   
      * PURPOSE    : To test if a knot vector is cyclic.
      *
      *             
      *
      * INPUT      : et     - Knot vector
      *              in     - Number of knots
      *              ik     - Polynomial order
      *
      * OUTPUT     : 
      *              jstat  - status messages  
      *                                        = 2 : Cyclic full freedom.
      *                                        = 1 : Cyclic NOT full freedom.
      *                                        = 0 : NOT cyclic.
      *                                        < 0 : Error.
      *
      * METHOD     : 1. Check multiplicity of start and end parameter value
      *              2. Check if knots before start parameter value and after
      *                 end parameter value are repeated in a cyclic way.
      *              3. Check that number of basis function not repeated
      *                 are at least ik.
      *
      *
      * REFERENCES :
      *
      *-                                      
      * CALLS      : 
      *
      * WRITTEN BY : Tor Dokken, SI, Oslo, Norway, April 1992
      *
      *********************************************************************
      */
{
  int    kleft;          /* Pointer into knot interval             */
  int    kmult1;         /* Multiplicity of start parameter value  */
  int    kmult2;         /* Multiplicity of end  parameter value   */
  int    ki;             /* Control variable in loop               */
  int    kpos = 1;       /* Position of error                      */
  int    kant;           /* Number of knots before start parameter value */
  int    kcyclic;        /* Flag telling if cyclic basis           */
  int    kstat;          /* Local status variable                  */
  
  double tperiode;       /* Periode of basis                       */
  
  /* Find multiplicity of et[ik-1] and et[in] */
  
  kleft = ik-1;
  
  kmult1 = s6knotmult(et,ik,in,&kleft,et[ik-1],&kstat);
  if(kstat<0) goto error;
  
  kleft = in;
  
  kmult2 = s6knotmult(et,ik,in,&kleft,et[in],&kstat);
  if(kstat<0) goto error;
  
  if (kmult1 != kmult2 || kmult1 == ik) goto noncyclic;
  
  kant = ik - kmult1;
  tperiode = et[in] -et[ik-1];  
  
  /* Test that the first kant knots are repetitions of the knots in-kant,...,in-1 */
  
  for (ki=0, kcyclic=1; ki<kant ; ki++)
    if (DNEQUAL((et[ki]+tperiode),et[in-kant+ki])) kcyclic = 0;
  
  /* Test that the last kant knots are repetions of knots ik,..,ik+kant-1 */
  
  for (ki=0; ki<kant ; ki++)
    if (DNEQUAL((et[ik+ki]+tperiode),et[in+kmult1+ki])) kcyclic = 0;
  
  if (kcyclic == 0) goto noncyclic;
  
  /* The basis should have at least kant+ik degrees of freedom to allow for
     a proper cyclic curve with ik degrees of real freadom since kant vertices
     are repeated */
  
  if (in<ik+kant) goto missing_freedom; 
  
  /* Cyclic with enough degrees of freedom */
  
  *jstat = 2;
  goto out;
  
  /* Cyclic with less than ik+kant degrees of freedom */
 missing_freedom:
  
  *jstat = 1;
  goto out;
  
  /* Noncyclic basis */
 noncyclic:
  
  *jstat = 0;
  goto out;
  
 error:
  *jstat = kstat;
  s6err("test_cyclic_knots",*jstat,kpos);
  goto out;
  
 out:
  
  return;
  
}


