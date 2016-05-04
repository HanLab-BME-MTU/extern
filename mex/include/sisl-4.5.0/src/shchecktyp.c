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
 * $Id: shchecktyp.c,v 1.2 2001-03-19 16:06:04 afr Exp $
 *
 */


#define SHCHECKTYPE

#include "sislP.h"                            


#if defined (SISLNEEDPROTOTYPES)
int
    shchecktype(SISLObject *pobj,double *parval)
#else
int shchecktype(pobj,parval)

     SISLObject *pobj;
     double *parval;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : To classify (in the extremal sense) a onedimensional 
*              object's behaviour at a given parameter value.
*
*
*
* INPUT      : pobj     - Pointer to a surface or curve object.
*              parval   - Parameter value(s).
*
*
* OUTPUT     : shchecktype :
*                            0 - No extremal point
*                            1 - Maximum point
*                            2 - Minimum point
*                            3 - Saddle point (used in surface only)
*                            4 - Extremum, test inconclusive.
*                          < 0 - error
*
*
* METHOD     : Test are performed on the first and second derivatives.
*
*
* REFERENCES :
*
*-
* CALLS      : s1221      - Evaluate curve.
*              s1421      - Evaluate surface.
*
* WRITTEN BY :Ulf J. Krystad, January 1991.
*
*********************************************************************
*/                                     
{
   
   int kstat;              /* Local status variable.            */
   int kleft1 = 0;         /* Knot navigator                    */
   int kleft2 = 0;         /* Knot navigator                    */
   int kder   = 2;         /* Flag, compute 2. derivative.      */
   int kdim   = 1;         /* Dimension is one !                */
   double sval[9];         /* Storing uppdated parametervalues. */
   double sval1[9];        /* Storing uppdated parametervalues. */
   double *snorm = sval+6; /* Dummy normal pointer              */
   double tmax;        	   /* Size of derivatives               */
   double tdet;            /* Size of Hessian determinant       */
   int mult = 0;           /* Knot multiplicity                 */
   double ttol = 1000000.0*REL_COMP_RES;
   /* --------------------------------------------------------- */
   
   if (pobj == SISL_NULL ||
       (pobj->iobj != SISLCURVE && pobj->iobj != SISLSURFACE))
      return -1;
   
   if (pobj->iobj == SISLCURVE )
   {
      /* Curve case */
      if (pobj->c1->idim != kdim) return -1;
      
      /* Get function values */
      s1221(pobj->o1->c1,kder,parval[0],&kleft1,sval,&kstat);
      if (kstat < 0) return -2;
      
      mult = s6knotmult(pobj->o1->c1->et,pobj->o1->c1->ik,
		     pobj->o1->c1->in,&kleft1,parval[0],&kstat);
      if (kstat < 0) return -2;
      
      if (mult >= pobj->o1->c1->ik - 1)
      {	 
	 /* Get left side function values */
	 s1227(pobj->o1->c1,kder,parval[0],&kleft1,sval1,&kstat);
	 if (kstat < 0) return -2;
	 /* Test function values */
	 if (sval[1] < -ttol && sval1[1] >ttol)       return 1;
	 else if (sval[1] > ttol && sval1[1] < -ttol) return 2;
	 else                                       return 4;
      }
      else
      {
	 /* Test if first derivative is zero */
	 if (fabs(sval[1]) > ttol) return 0;
	 
	 /* Test if max, min or inconclusive point */
	 if      (sval[2] < -ttol) return 1;
	 else if (sval[2] >  ttol) return 2;
	 else                      return 4;
      }
   }
   else
   {
      /* Surface case */
      if (pobj->s1->idim != kdim) return -1;
      
      /* Get function values */
      
      s1421(pobj->o1->s1,kder,parval,&kleft1,&kleft2,
	    sval,snorm,&kstat);
      if (kstat < 0) return -2;
      
      /* Test function values */
      
      /* Test if first derivative is zero */
      tmax = sqrt(sval[1]*sval[1] + sval[2]*sval[2]);
      if (tmax > ttol) return 0;
      
      /* Test if max, min ,saddle or inconclusive point */
      tdet = (sval[3]*sval[5] - sval[4]*sval[4]);
      if      (tdet < -ttol) return 3;
      else if (tdet <  ttol) return 4;
      else if (sval[3] < DZERO) return 1;
      else return 2;
      
   }
   
}	 
   
