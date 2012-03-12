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
 * $Id: s6decomp.c,v 1.2 2001-03-19 15:59:01 afr Exp $
 *
 */


#define S6DECOMP

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s6decomp(double ea[],double gx[],double eb1[],double eb2[],double eb3[],int *jstat)
#else
void s6decomp(ea,gx,eb1,eb2,eb3,jstat)
     double ea[];
     double gx[];
     double eb1[];
     double eb2[];
     double eb3[];
     int    *jstat;
#endif
/*
***********************************************************************
*
************************************************************************
*
*   PURPOSE : In dimention tree we change basis from a euclides bases
*             to eb1, eb2, eb3.
*
*
*   INPUT   : ea   - The orginale vector.
*             eb1  - the first bases vector.
*             eb2  - the first bases vector.
*             eb3  - the first bases vector.
*
*
*   
*   OUTPUT  : gx    - The new vector
*             jstat - Status variable.
*                       < 0 : error
*                       = 0 : ok
*                       > 0 : warning
*
*
*   METHOD  : 
*
*
*   REFERENCES : 
*-
*   CALLS      :
*
*   WRITTEN BY : Arne Laksaa, SI, 89-07.
*
************************************************************************
*/
{
  int kstat =0;       /* Local status variable.    */
  int ki;             /* Counter.                  */
  int n1[3];          /* Array for use in lufac.   */
  double sc[9],se[3]; /* Matrix and help vector.   */
  
  
  /* Copy new bases into local matrix.  */
  
  memcopy(sc,eb1,3,double);
  memcopy(sc+3,eb2,3,double);
  memcopy(sc+6,eb3,3,double);
  
  
  s6lufacp(sc,n1,3,&kstat);
  if (kstat < 0) goto error;
  else if (kstat > 0) goto warn1;
  
  for (ki=0; ki<3; ki++)
    {                      
      se[0] = se[1] = se[2] = DZERO;
      se[ki] = (double)1;
      
      s6lusolp(sc,se,n1,3,&kstat);
      if (kstat < 0) goto error;
      else if (kstat > 0) goto warn1;
      
      gx[ki] = s6scpr(ea,se,3);
    }
  
  /* Change of bases performed.  */
  
  *jstat = 0;
  goto out;

/* Singular equation system.  */

warn1 : *jstat = 1;
        goto out;

/* Error in subrutines.  */

error: *jstat = kstat;
        s6err("s6decomp",*jstat,0);
        goto out;

 out: ;
}

                                    

                                        
