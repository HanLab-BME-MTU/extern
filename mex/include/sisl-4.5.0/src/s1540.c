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
 * $Id:
 *
 */
#define S1540


#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1540  ( double   et[],
	 int      ik,
	 int      in,
	 double   ax[],
	 int      im,
	 int      ider,
	 double   ebder[],
	 int      ileft[],
	 int     *jstat )
#else
void
s1540 ( et, ik, in, ax, im, ider, ebder, ileft, jstat )
  double   et[];
  int      ik;
  int      in;
  double   ax[];
  int      im;
  int      ider;
  double   ebder[];
  int      ileft[];
  int     *jstat;
#endif
/*
*********************************************************************
*
* PURPOSE:	Given a (polynomial) spline basis
*		(with k-regular knot vector)
*		and an array of parameters ax[0, ..., im-1]
*		calculate the basis values (and first derivatives if ider = 1)
*		and fill them in the array ebder[].
*
* INPUT:	ik     -  spline order
*		in     -  number of spline functions (only for checking ax).
*		ax     -  array of parameters
*		im     -  length of the array ax[]
*
* OUTPUT:	ebder  -  array containing the basis values
*		          ider = 0:  B(ax[0   ],i),...,B(ax[0   ],i+ik-1)
*		                     B(ax[1   ],i),...,B(ax[1   ],i+ik-1)
*		                     B(ax[im-1],i),...,B(ax[im-1],i+ik-1)
*
*		          ider = 1:  B (ax[0   ],i),...,B (ax[0   ],i+ik-1)
*		                     B'(ax[0   ],i),...,B'(ax[0   ],i+ik-1)
*		                     B (ax[1   ],i),...,B (ax[1   ],i+ik-1)
*		                     B'(ax[1   ],i),...,B'(ax[1   ],i+ik-1)
*						...............
*		                     B (ax[im-1],i),...,B (ax[im-1],i+ik-1)
*		                     B'(ax[im-1],i),...,B'(ax[im-1],i+ik-1)
*
*		Ie. if we have positions and derivatives, then
*		for each ax[i] first the positions come and then the
*		derivatives.
*              jstat    - Status messages
*
*                         = 0 : Ok.
*                         < 0 : Error.
*
* IMPORTANT:	ebder is filled differently by the corresponding SISL routine
*		(s1504):
*			 B (ax[0   ],i     ), B'(ax[0   ],i     ),
*			 B (ax[0   ],i+1   ), B'(ax[0   ],i+1   ),
*			        ...............
*		         B (ax[0   ],i+ik-1), B'(ax[0   ],i+ik-1),
*			 B (ax[1   ],i     ), B'(ax[1   ],i     ),
*			        ...............
*			 B (ax[im-1],i+ik-1), B'(ax[im-1],i+ik-1),
*
*
* WRITTEN BY: 	Geir Westgaard, SINTEF, Oslo, November 1999
*
*********************************************************************
*/
{
  int kstat = 0;      /* Local status variable.                          */
  int kpos  = 0;      /* The position of error.                          */
  int i, k;           /* Control variables in for loops and for stepping
                         through arrays.                                 */
  int size;           /* (ider+1) * ik.                                  */
  double *eder = SISL_NULL;/* B-spline evaluationas at a single value.        */
  double tmpeder[10]; /* meaning: tmpeder[(ider+1)*ik]
			 and assuming: ider <= 1 and ik <= 5             */



   /* Check the input. */

    if ( ider < 0 || ider > 1        ) goto err10;
    if ( ik < 2 || ik > 5            ) goto err10;
    if ( im  < 0                     ) goto err10;


   /* Set local variables. */

   size = (ider + 1)*ik;

   eder = ebder;


   if ( ider == 0 )
   {
      for( k = 0; k < im; k++, eder += ik )
      {
	 s1220( et, ik, in, ileft + k, ax[k], ider, eder, &kstat );

	 if ( kstat < 0 ) goto error;
      }
   }
   else
   {

      for( k = 0; k < im; k++, eder += size )
      {
	 s1220( et, ik, in, ileft + k, ax[k], ider, tmpeder, &kstat );

	 if ( kstat < 0 ) goto error;

	 for ( i = 0; i < ik; ++i )
	 {
	    eder[i     ] = tmpeder[2*i    ];
	    eder[i+ik  ] = tmpeder[2*i + 1];
	 }
      }
   }




  /* Successful computations.  */

   *jstat = 0;
   goto out;

/* Error in input. */
err10: *jstat = -10;
s6err( "s1540", *jstat, kpos );
goto out;

/* Error in lower level routine.  */
error:  *jstat = kstat;
s6err( "s1540", *jstat, kpos );
goto out;




out: return;

}
