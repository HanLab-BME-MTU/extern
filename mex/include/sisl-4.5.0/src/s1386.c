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
 * $Id: s1386.c,v 1.5 2005-02-28 09:04:48 afr Exp $
 *
 */


#define S1386

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1386(SISLSurf *ps,int ider1,int ider2,SISLSurf **rsnew,int *jstat)
#else
void s1386(ps,ider1,ider2,rsnew,jstat)
     SISLSurf 	*ps;
     int  	ider1;
     int  	ider2;
     SISLSurf 	**rsnew;
     int  	*jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : To express the (ider1,ider2)-th derivative of a B-pline
*              surface as a B-spline surface with order (ider1,ider2)
*              less than the original B-spline surface.
*
*
*
* INPUT      : ps       - Surface to make a derivative from.
*              ider1    - The derivative to be produced in the first
*                         parameter direction 0 <= ider1 < p2 ->ik1
*              ider2    - The derivative to be produced in the second
*                         parameter direction 0 <= ider2 < p2 ->ik2
*
*
*
* OUTPUT     : rsnew    - The result of the (ider1,ider2) differentiation of ps.
*              jstat    - status messages
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : For non-rational surfaces, we treate the surface as a
*                curve in both parameter directions and utilize the
*                curve differentiation method twice.
*              For nurbs, we find the new orders and dimensions, and then
*                interpolate a set of points on the derivative surface.
*
*
* REFERENCES : Larry L. Schumaker, Spline Functions:Basic Theory.  p.195.
*              Carl de Boor, A Practial Guide to Spline.           p.139.
*
*-
* CALLS      : newCurve  - Allocate space for a new curve-object.
*              freeCurve - Free space occupied by given curve-object.
*              s1705.c   - Remove unnecesery knots and vertices from a curve.
*
* WRITTEN BY : Arne Laksaa, SI, 88-06.
* REVISED BY : Johannes Kaasa, SI, May 92 (Introduced NURBS. Remark that
*              it is then impossible to make the derivative in a higher
*              dimensional space. We must treat the surface as a collection
*              of seperate curves in each direction, this will of course make
*              the function less efficient).
* REVISED BY : Christophe Birkeland, SI, July 1992 (memcopy-statement line 276
*              is changed, call s6chpar line 283: kn1 & kn2 interchanged)
* REVISED BY : Christophe Birkeland, SI, August 1992
*              The NURBS algorithm is changed completely. First we find the new
*              orders and dimensions of the surface, then we estimate points
*              on the derivative surface, and at last we interpolate the
*              resulting derivative surface.
* REVISED BY : Johannes Kaasa, SI, Aug. 92 (Established correct order and interior
*              knot multiplicity in the NURBS case. In addition I have used s1891
*              to solve the interpolation matrix in both parameter directions).
* Revised by : Paal Fugelli, SINTEF, Oslo, Norway, Oct. 1994.  Changed icopy
*              value from 0 to 2 in call to newSurf for rationals and free'ed 'rat',
*              'par1', 'par2', 'der1' and 'der2' at exit to remove memory
*              leakage problems.  Initialized 'rsnew'.  Increased the size
*              of 'scoef' for non-rationals by 'kn2*ider1*ps->idim', since the
*              first call to s1720() will increase "in" by ider1.
*              NOTE: closed and periodic case isn't handled correctly by this
*              routine and it will have to be handled separatly for NURBS and
*              polynomial case.
*
**********************************************************************/
{
  SISLCurve *qc1 = SISL_NULL;      /* Temporary curve                     */
  SISLCurve *qc2 = SISL_NULL;      /* Temporary curve                     */
  SISLCurve *qc3 = SISL_NULL;      /* Temporary curve                     */
  SISLCurve *qc4 = SISL_NULL;      /* Temporary curve                     */

  int kk1;                    /* Order in first parameter direction  */
  int kk2;                    /* Order in second parameter direction */
  int kn1;                    /* NumberOrder in first parameter direction */
  int kn2;                    /* Order in second parameter direction */
  int kdim;                   /* Dimension used in temporary calc    */
  int kstat;                  /* Local parameter value               */
  int kpos=0;                 /* Position of error                   */
  double *st1 = SISL_NULL;         /* Pointer to knot vector              */
  double *st2 = SISL_NULL;         /* Pointer to knot vector              */
  double *scoef = SISL_NULL;       /* Pointer to coefficients             */

  /* NURBS variables: */

  SISLSurf *rat=SISL_NULL;         /* The denominator surface.            */
  double *ratcoef=SISL_NULL;       /* The vertices of rat.                */
  int ki, kj, kl, km;         /* Index in for loop.                  */
  int rdim;                   /* Rational dimension.                 */

  double eps;                 /* Knot equality resolution.               */
  int multadd;                /* Added multiplicity of interior knots.   */
  int mult;                   /* Original multiplicity of the knots.     */
  int oldprev;                /* Previous index in the original knots.   */
  int oldcurr;                /* Current index in the original knots.    */
  int newcurr;                /* Current index in the new knots.         */
  double denom;               /* The denominator.                        */
  double limit;               /* Used as a loop control parameter        */
  int left1 = 0;              /* Interval indicator.                     */
  int left2 = 0;              /* Interval indicator.                     */
  double par[2];              /* Parameter values used for interpolation */
  double *par1 = SISL_NULL;        /* Parameter values used for interpolation */
  double *par2 = SISL_NULL;        /* Parameter values used for interpolation */
  int *der1 = SISL_NULL;           /* The derivative indicators (= 0).        */
  int *der2 = SISL_NULL;           /* The derivative indicators (= 0).        */
  double *deriv = SISL_NULL;       /* The derivates returned by s1221.        */
  double *tau = SISL_NULL;         /* Interpolation points.                   */
  int iopen = TRUE;           /* Open flag.                              */
  int in;                     /* Number of vertices in interpolation.    */
  int inlr = 0, inrc = 0;
  int indx1, indx2;

  /* Check that we have a surface to differentiable. */

  *rsnew = SISL_NULL;  /* Must be valid incase of error exit. */

  if (!ps) goto err150;

  /* Check that the number of derivation is legal. */

  if ( ider1 < 0  ||  ider2 < 0 )  goto err156;

  if (ps->ikind == 2 || ps->ikind == 4)
  {

     /* NURBS surface. */

     kdim=ps->idim;
     rdim = kdim + 1;

     /* Make the denominator surface. */

     ratcoef = newarray(ps->in1 * ps->in2, DOUBLE);
     limit = ps->in1*ps->in2;
     for (kj=0; kj < limit ; kj++)
	ratcoef[kj] = ps->rcoef[(kj + 1)*rdim - 1];
     rat = newSurf(ps->in1,ps->in2, ps->ik1,ps->ik2, ps->et1,
		   ps->et2,ratcoef, 1, 1, 1);
     if (ratcoef != SISL_NULL) freearray(ratcoef);



     /* Make resolution for testing of knot equality. */

     eps = fabs(ps->et1[ps->in1] - ps->et1[ps->ik1 - 1])*REL_PAR_RES;

     /* Make the new knot vector st1 in first direction. */

     kn1 = 0;
     kk1 = (ider1 + ider2 + 1)*ps->ik1 - (ider1 + ider2);
     multadd = (ider1+ider2)*ps->ik1 - ider2;

     st1 = newarray((2 + ps->in1 - ps->ik1)*kk1, DOUBLE);

     for (newcurr=0; newcurr<kk1; newcurr++)
     {
	st1[newcurr] = ps->et1[ps->ik1 - 1];
        kn1++;
     }

     oldcurr = ps->ik1;
     oldprev = oldcurr;
     limit = ps->et1[ps->in1] - eps;
     while (ps->et1[oldcurr] < limit)
     {
	mult = 0;
	while ((ps->et1[oldcurr] - ps->et1[oldprev]) < eps)
	{
	   oldcurr++;
	   mult++;
	}
	mult += multadd;
	if (mult > kk1) mult = kk1;
	for (kj=0; kj<mult; kj++)
	{
	   st1[newcurr + kj] = ps->et1[oldprev];
	   kn1++;
	}
	newcurr += mult;
	oldprev = oldcurr;
     }
     for (kj=0; kj<kk1; kj++)
	st1[newcurr + kj] = ps->et1[ps->in1];

     /* Resize new knot vector st1 */

     st1 = increasearray(st1,kn1+kk1,DOUBLE);
     if (st1 == SISL_NULL) goto err101;

     /* Calculate parameter values and derivate indicators. */

     s1890(st1, kk1, kn1, &par1, &der1, &kstat);
     if (kstat < 0) goto error;



     /* Knot vector in second parameter direction */

     /* Make resolution for testing of knot equality in second parameter
      * direction. */

     eps = fabs(ps->et2[ps->in2] - ps->et2[ps->ik2 - 1])*REL_PAR_RES;

     /* Make the new knot vector st2. */

     kn2 = 0;
     kk2 = (ider1 + ider2 + 1)*ps->ik2 - (ider1 + ider2);
     multadd = (ider1 + ider2)*ps->ik2 - ider1;

     st2 = newarray((2 + ps->in2 - ps->ik2)*kk2, DOUBLE);

     for (newcurr=0; newcurr<kk2; newcurr++)
     {
	st2[newcurr] = ps->et2[ps->ik2 - 1];
        kn2++;
     }

     oldcurr = ps->ik2;
     oldprev = oldcurr;
     limit = ps->et2[ps->in2]-eps;
     while (ps->et2[oldcurr] < limit)
     {
	mult = 0;
	while ((ps->et2[oldcurr] - ps->et2[oldprev]) < eps)
	{
	   oldcurr++;
	   mult++;
	}
	mult += multadd;
	if (mult > kk2) mult = kk2;
	for (kj=0; kj<mult; kj++)
	{
	   st2[newcurr + kj] = ps->et2[oldprev];
	   kn2++;
	}
	newcurr += mult;
	oldprev = oldcurr;
     }

     for (kj=0; kj<kk2; kj++)
	st2[newcurr + kj] = ps->et2[ps->in2];

     /* Resize new knot vector st2 */

     st2 = increasearray(st2,kn2+kk2,DOUBLE);
     if (st2 == SISL_NULL) goto err101;

     /* Calculate parameter values and derivate indicators. */

     s1890(st2, kk2, kn2, &par2, &der2, &kstat);
     if (kstat < 0) goto error;



     /* ------------------------------- */
     /* Calculate interpolation points. */
     /* ------------------------------- */

     deriv = newarray((ider1 + 1)*(ider2 + 1)*kdim, DOUBLE);
     tau = newarray(rdim*kn1*kn2, DOUBLE);
     if (tau == SISL_NULL) goto err101;

     km = 0;
     for (kj=0 ; kj<kn1 ; kj++)
	for (kl=0; kl<kn2; kl++)
	{
	   par[0]=par1[kj];
	   par[1]=par2[kl];
	   s1424(rat, 0, 0, par, &left1, &left2, &denom, &kstat);
	   if (kstat < 0) goto error;
	   denom = pow(denom, (ider1 + ider2 + 1));
	   s1424(ps, ider1, ider2, par, &left1, &left2, deriv, &kstat);
	   if (kstat < 0) goto error;
	   for (ki=0; ki<kdim; ki++)
	      tau[km++] =
		 deriv[(ider2*(ider1+1)+ider1)*kdim + ki]*denom;
	   tau[km++] = denom;
	}

     /* Solve the interpolation equation in the second parameter direction. */

     s1891(par2, tau, rdim, kn2, kn1, der2, iopen, st2, &scoef, &in, kk2,
          inlr, inrc, &kstat);
     if (kstat < 0) goto error;
     if (in != kn2) goto error;

     /* Transpose scoef and put in tau. */

     km = 0;
     for (kj=0; kj<kn2; kj++)
       {
	 indx1 = kj*rdim;
         for (kl=0; kl<kn1; kl++)
	   {
	     indx2 = kl*kn2*rdim;
             for (ki=0; ki<rdim; ki++)
                tau[km++] = scoef[indx2 + indx1 + ki];
	   }
       }

     if (scoef != SISL_NULL) freearray(scoef);

     /* Solve the interpolation equation in the first parameter direction. */

     s1891(par1, tau, rdim, kn1, kn2, der1, iopen, st1, &scoef, &in, kk1,
          inlr, inrc, &kstat);
     if (kstat < 0) goto error;
     if (in != kn1) goto error;

     *rsnew = newSurf(kn1,kn2,kk1,kk2,st1,st2,scoef,ps->ikind,3,2);
  }
  else
  {
     /* Not NURBS. */

     kk1 = ps -> ik1;
     kk2 = ps -> ik2;
     kn1 = ps -> in1;
     kn2 = ps -> in2;

     /* Create curve representing the surface a a curve in the second parameter
        direction, copy input arrays */

     kdim = (ps->in1)*(ps->idim);

     qc1 = newCurve(ps->in2,ps->ik2,ps->et2,ps->ecoef,1,kdim,1);
     if (qc1 == SISL_NULL) goto err101;

     /* Make the derivative in the second parameter direction */

     s1720(qc1,ider2,&qc2,&kstat);
     if (kstat<0) goto error;

     /* Remember new knot vector in second parameter direction */

     kk2 = qc2 -> ik;
     kn2 = qc2 -> in;
     st2 = newarray(kk2+kn2,DOUBLE);
     if (st2 == SISL_NULL) goto err101;

     memcopy(st2,qc2->et,kk2+kn2,DOUBLE);

     /* Allocate space for turned parameter directions */

     scoef = newarray((kn1*kn2 + kn2*ider1)*(ps->idim),DOUBLE);
     if (scoef == SISL_NULL) goto err101;

     /* Turn parameter directions */

     s6chpar(qc2->ecoef,kn1,kn2,(ps->idim),scoef);


     /* Represent the surface as curve using the first knot vector */

     kdim = kn2*(ps->idim);

     qc3 = newCurve(ps->in1,ps->ik1,ps->et1,scoef,1,kdim,1);
     if (qc3 == SISL_NULL) goto err101;

     /* Make the derivative in the first parameter direction */

     s1720(qc3,ider1,&qc4,&kstat);
     if (kstat<0) goto error;


     /* Remember new knot vector in first parameter direction */

     kk1 = qc4 -> ik;
     kn1 = qc4 -> in;
     st1 = newarray(kk1+kn1,DOUBLE);
     if (st1 == SISL_NULL) goto err101;

     memcopy(st1,qc4->et,kk1+kn1,DOUBLE);


     /* Turn parameter directions of coefficients to match surface */

     s6chpar(qc4->ecoef,kn2,kn1,(ps->idim),scoef);

     /* Create surface object containing the differentiated of the surface */

     *rsnew = newSurf(kn1,kn2,kk1,kk2,st1,st2,scoef,ps->ikind,(ps->idim),1);
     if (*rsnew == SISL_NULL) goto err101;

  }

  *jstat = 0;
  goto out;

  /* Error in space allocation */

 err101: *jstat = -101;
  s6err("s1386",*jstat,kpos);
  goto out;

  /* Error. No surface to differentiate.  */

 err150:
  *jstat = -150;
  s6err("s1386",*jstat,kpos);
  goto out;

  /* Error. Illegal number of derivatives.  */

 err156:
  *jstat = -156;
  s6err("s1386",*jstat,kpos);
  goto out;

  /* Error. Error in lower level function. */

 error:  *jstat = kstat;
  s6err("s1386",*jstat,kpos);
  goto out;


  /* Free local used memory. */

 out:
  if (qc1 != SISL_NULL) freeCurve(qc1);
  if (qc2 != SISL_NULL) freeCurve(qc2);
  if (qc3 != SISL_NULL) freeCurve(qc3);
  if (qc4 != SISL_NULL) freeCurve(qc4);
  if (ps->ikind != 2 && ps->ikind != 4)
  {
     if (st1 != SISL_NULL) freearray(st1);
     if (st2 != SISL_NULL) freearray(st2);
     if (scoef != SISL_NULL) freearray(scoef);
  }
  else
  {
    if (par1)  freearray(par1);
    if (par2)  freearray(par2);
    if (der1)  freearray(der1);
    if (der2)  freearray(der2);

    if (rat) freeSurf(rat);

    if (tau != SISL_NULL) freearray(tau);

    if (deriv != SISL_NULL) freearray(deriv);
  }

  return;
}
