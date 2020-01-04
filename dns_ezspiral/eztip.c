/* ------------------------------------------------------------------------- 
 * eztip.c -- Spiral Tip Routines 
 *
 * Copyright (C) 1992, 1994, 1997, 1998 Dwight Barkley                 
 *
 * RCS Information
 * ---------------------------
 * $Revision: 3.2.1.1 $
 * $Date: 2007/05/07 09:50:42 $
 * ------------------------------------------------------------------------- */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ezspiral.h"

/*---------------------------------------------------------------------------*
 * These subroutines compute the location of the spiral tips. The tip (x,y) is
 * defined implicitly by f1(x,y) = f2(x,y) = 0. The functions f1 and f2 are
 * defined at the end of this file and can in principle be almost any
 * functions. Currently they are such as to define the tip as the intersection
 * of a pair of u,v contours. Newton iterations are used to find the root
 * (x,y). Derivatives are determined by finite differencing.  The u and v
 * fields are approximated (interpolated) with polynomials of order
 * N_ORDER. **THE CODE ASSUMES N_ORDER IS ODD.** N_ORDER=3 gives accurate
 * results.
 *
 * Note that each time Find_tips () is called, the entire grid is searched for
 * possible tips. This is a bit inefficient if there is only one tip. However
 * timing tests show that it takes little time and allows for finding multiple
 * tips.
 *---------------------------------------------------------------------------*/

#define U_CONT  0.5                    /* U and */
#define V_CONT (0.5*a_param-b_param)   /* V contour values defining the tip */
/* #define V_CONT 0.4 */  /* changed Bjorn Sandstede: uncomment for FHN */

#define NORDER 3         /* interpolation order for fields (MUST BE ODD)     */
#define DX     0.05      /* finite diff used for computing derivatives       */
#define DY     0.05      /* finite diff used for computing derivatives       */
#define TOL2   1.e-10    /* tol (squared) for newton convergence             */
#define TOOBIG 10.       /* max size of eps2 in newton (to prevent blowup)   */
#define MXNWT  8         /* max number of newton iterations                  */
#define MAX_TIPS 100000  /* max number of tips                               */

#define PRNT_WARN(verbose,level,message) if(verbose>=level) printf(message);

/* 
 * Private functions 
 * ----------------- */
static int   Newton      (Real *px, Real *py);
static Real  Eval_func   (Real (*f)(int i, int j), Real x, Real y, int new);
static Real  hh          (Real x_local, int i); 
static Real  f1          (int i, int j);
static Real  f2          (int i, int j);

/* ========================================================================= */

void Tip_ini (void)
{
  tip_x =(Real *) malloc((unsigned)(MAX_TIPS)*sizeof(Real));
  tip_y =(Real *) malloc((unsigned)(MAX_TIPS)*sizeof(Real));
  tip_x[0] = 0.;
  tip_y[0] = 0.;
  ntips = 0;
}
/* ========================================================================= */

void Find_tips (void) 
{
  Real x,y;
  int i,j,count=0;
  unsigned int index;

  /* Search the grid (except close to edges) */
  for(i=2;i<NX-1;i++) {
    for(j=2;j<NY-1;j++) {
      
      index = ((f1(i  ,j  ) >= 0.)    ) |               /* index for u field */
	      ((f1(i+1,j  ) >= 0.) <<1) |
	      ((f1(i  ,j+1) >= 0.) <<2) |
	      ((f1(i+1,j+1) >= 0.) <<3) ;
      
      if ((index != 0) && (index != 15)) {        /* non-trivial index for u */

	index = ((f2(i  ,j  ) >= 0.)    ) |             /* index for v field */
	        ((f2(i+1,j  ) >= 0.) <<1) |
	        ((f2(i  ,j+1) >= 0.) <<2) |
	        ((f2(i+1,j+1) >= 0.) <<3) ;
      
	if ((index != 0) && (index != 15)) {      /* non-trivial index for v */

	  /* Both u and v have non-trivial index. Possible tip. Take initial
             guess to be center of this square */

	  x = (Real)i + 0.5;    
	  y = (Real)j + 0.5;

	  /* Use Newton iteration to find tip.  Do not do this if above guess
             is too close to the last computed tip at this time step */

	  if( (ntips==0) || (count==0) ||
	      ((x-tip_x[ntips-1])*(x-tip_x[ntips-1]) + 
	       (y-tip_y[ntips-1])*(y-tip_y[ntips-1])) > 4 ) {

	    if( Newton(&x, &y) ) {
	      tip_x[ntips] = x;
	      tip_y[ntips] = y;
	      ntips++;
	      count++;
	      if(write_tip) Write_tip_data (grid_h*(x-1.), grid_h*(y-1.));
	    }
	  }
	}
      }
    }
  }
  if(verbose>=2) printf("number of tips in one field= %d\n", count);
}
/* ========================================================================= */

static int Newton (Real *px, Real *py) 
{
  Real x=*px, y=*py, f1_val, f2_val, df1dx, df1dy, df2dx, df2dy, jac, nrm2;
  int nwt=0, new=1, old=0;    /* new & old determine whether new
			       * interpolating polynomials are
			       * generated */

  while(nwt<MXNWT) {

    f1_val = Eval_func(f1,x,y,new);           /* compute f1 at current (x,y) */
    f2_val = Eval_func(f2,x,y,old);           /* compute f2 at current (x,y) */
    nrm2 = f1_val*f1_val + f2_val*f2_val;               /* norm of functions */

    if(nrm2>TOOBIG) return(0);

    if(nrm2<TOL2) {                                              /* Success! */
      *px = x;  
      *py = y;
      return(1);
    }

    /* Compute new x and y by taking one Newton step */
    df1dx = (Eval_func(f1,x+DX,y,new) - f1_val)/DX;
    df2dx = (Eval_func(f2,x+DX,y,old) - f2_val)/DX;
    df1dy = (Eval_func(f1,x,y+DY,new) - f1_val)/DY;
    df2dy = (Eval_func(f2,x,y+DY,old) - f2_val)/DY;
    jac   = df1dx*df2dy - df1dy*df2dx;
    x -= ( df2dy*f1_val - df1dy*f2_val)/jac;
    y -= (-df2dx*f1_val + df1dx*f2_val)/jac;
    nwt++;
  }

  /* Too many Newton iterations */
  PRNT_WARN(verbose,2,"Warning: nwt > MXNWT in Newton()\n");
  return(0);
}
/* ========================================================================= */

Real Eval_func (Real (*f)(int i, int j), Real x, Real y, int new)
{
  /* Uses Lagrangian interpolation to find value of function f at point
   * (x,y).  If new=1, then this is a new (x,y) point and new interpolating
   * polynomials are computed.  */

  static Real x_local, y_local, hx[NORDER+1], hy[NORDER+1]; 
  static int ix_0, iy_0;
  Real sum=0.;
  int i,j;

  if(new) { 
    /* Generate new polynomials */
    ix_0 = (int)x - ((NORDER+1)/2 -1);
    iy_0 = (int)y - ((NORDER+1)/2 -1);
    x_local = x - ix_0;
    y_local = y - iy_0;
    for(i=0;i<(NORDER+1);i++) {
      hx[i] = hh(x_local,i);
      hy[i] = hh(y_local,i);
    }
  }

  /* Test for out of range */
  if(ix_0<1||ix_0>NX-NORDER||iy_0<1||iy_0>NY-NORDER) {
    PRNT_WARN(verbose,2,"Warning: (x,y) out of range in Eval_func() \n");
    return(TOOBIG+1.);
  }

  /* Interpolate function */
  for(i=0;i<(NORDER+1);i++) {
    for(j=0;j<(NORDER+1);j++) {
      sum += hx[i]*hy[j]*(*f)(ix_0+i,iy_0+j);
    }  
  }
  return(sum);
}
/* ========================================================================= */

Real hh (Real x_local, int i)
{
  /* Computes the ith Lagrange interpolating polynomial at point x_local.  */

  Real product=1.;
  int j;

  for(j=0;j<(NORDER+1);j++) 
    if(j!=i) product *= (x_local-(Real)j)/(Real)(i-j);

  return(product);
} 
/* ========================================================================= */


/* Define Functions f1 and f2 
 *
 * Here f1 and f2 are such as to define the spiral 
 * tip as the intersection of the contours: 
 * u=U_CONT and v=V_CONT 
 * ----------------------------------------------- */

Real f1 (int i, int j) 
{
  return(U(i,j)-U_CONT);
}
/* ========================================================================= */

Real f2 (int i, int j) 
{
  return(V(i,j)-V_CONT);
}
/* ========================================================================= */

