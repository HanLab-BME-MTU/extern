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


#define S1951

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
   s1951(double etau[], double ecoef[], int in, int ik, int idim, 
	 int ilend, int irend, int incont, double efac[])
#else
void s1951(etau, ecoef, in, ik, idim, ilend, irend, incont, efac)
   double etau[];
   double ecoef[];
   int    in;
   int    ik;
   int    idim;
   int    ilend;
   int    irend;
   int    incont;
   double efac[];
#endif     
/*
*********************************************************************
* 
* PURPOSE    : Multiply the coefficients by dtau(-1/2) and express 
*              the incont last coefficients as a weighted sum
*              of the incont first coeffecients. The weights are given
*              in efac.
* 
* 
* INPUT      : ecoef  - Coefficients of spline curve.
*              in     - Number of coefficients.p in
*              idim   - Dimension of geometry space.
*              incont - Number of continuity conditions, i.e. number of
*                       coefficients at the end to be expressed by
*                       coefficients at the start.
*              efac   - Factors.
*              
*
* 
* OUTPUT     : ecoef  - Coefficients of spline curve.
*             
* 
* METHOD     : 
*
*
* REFERENCES : 
*              
*
* USE        :
*
*-
* CALLS      :   
*
* WRITTEN BY : Vibeke Skytt,  SINTEF Oslo, 01.95.
*
*********************************************************************
*/
{
   int ki, kj, kr;   /* Counters.  */
   int kstop;
   double tw;

   /* Multiply the part of ec pointed to by kstart and kstop by the
      corresponding parts of the square matrix dtau(-1/2).   */
   
   for (kstop=in-MAX(incont,irend), ki=ilend; ki<kstop; ki++)
     {
	tw = sqrt((double)ik/(etau[ki+ik] - etau[ki]));
	for (kj=0; kj<idim; kj++)
	  ecoef[ki*idim+kj] *= tw;
     }  
   
   /* Express the incont last coefficients by the incont first ones given
      the factors stored in efac. See s1947. */
   
   for (ki=0; ki<incont; ki++)
   {
      for (kr=0; kr<idim; kr++)
      {
	 ecoef[(in-ki-1)*idim+kr] = DZERO;
	 for (kj=0; kj<=ki; kj++)
	    ecoef[(in-ki-1)*idim+kr] += ecoef[kj*idim+kr]*efac[ki*incont+kj];
      }
   }
   
}
   
