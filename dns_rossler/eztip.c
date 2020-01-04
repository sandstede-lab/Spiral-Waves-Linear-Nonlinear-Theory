/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*                                                                           */
/*                                 EZTIP.C                                   */
/*                Spiral-Tip Finding Subroutines for EZ-SPIRAL               */
/*                                                                           */
/*                                Version 2                                  */
/*                                                                           */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* 
 * $Revision: 2.3 $
 * $Date: 1997/11/13 23:18:44 $
 */
/*---------------------------------------------------------------------------*/

#include "ezspiral.h"

/*---------------------------------------------------------------------------*/
/*                                                                           */
/* These subroutines compute the location of the spiral tip. The tip (x,y)   */
/* is defined implicitly by f1(x,y) = f2(x,y) = 0. The functions f1 and f2   */
/* are defined at the end of this file and can in principle be almost any    */
/* functions. Currently they are such as to define the tip as the            */
/* intersection of a pair of u,v contours. Newton iterations are used to     */
/* find the root (x,y). Derivatives are determined by finite differencing.   */
/* The u and v fields are approximated (interpolated) with polynomials of    */
/* order N_ORDER. **THE CODE ASSUMES N_ORDER IS ODD.**  N_ORDER=3 gives      */
/* accurate results. In general, there should be no need to change the       */
/* parameters defined below.                                                 */
/*                                                                           */
/*---------------------------------------------------------------------------*/
                      /*                                                     */
#define NORDER 3      /* interpolation order for (u,v) fields (MUST BE ODD)  */
#define DX     0.01   /* finite diff used for computing derivatives          */
#define DY     0.01   /* finite diff used for computing derivatives          */
#define TOL2   1.e-10 /* tol (squared) for convergence of newton iterations  */
#define TOOBIG 20.    /* max size of eps2 in newton (to prevent blowup)      */
#define MXNWT  50     /* max number of newton iterations                     */
                      /*                                                     */

/* What to do if eztip fails to find a tip.  Here I set the do nothing except
 * last_failed = 1 so the finder knows to re-search the grid next time.
 */
#define FAIL_RETURN {last_failed = 1; return;}

/* External variables */
extern precision **u, **v, **w, *tip_x, *tip_y,
                 a_param, b_param, c_param;
extern int       ndim, ntips, verbose;

/* Functions local to eztip.c */
static int        tip_ini     (precision *x, precision *y); 
static precision  eval_func   (precision (*f)(int i, int j), 
			       precision x, precision y, int new);
static precision  hh          (precision x_local, int i); 
static precision  f1          (int i, int j);
static precision  f2          (int i, int j);

/*---------------------------------------------------------------------------*/

void tip_finder()

     /*  Main tip finding routine  */
{
  precision x, y, f1_val, f2_val, df1dx, df1dy, df2dx, df2dy, jac, nrm2;
  int nwt=0, restart=0, new=1, old=0;
  static int last_failed = 0;
  /* new & old determine whether new interpolating polynomials are generated */
  
  /* Initialization */
  if(ntips==0 || last_failed) {                          /* no guess for tip */
    if(!tip_ini(&x,&y)) FAIL_RETURN        
  }
  else if(ntips==1) {                      /* use first tip as initial guess */
    x = tip_x[0];
    y = tip_y[0];
  }
  else {                                    /* extrapolate for initial guess */
    x = 2.*tip_x[ntips-1] - tip_x[ntips-2];
    y = 2.*tip_y[ntips-1] - tip_y[ntips-2];
  }

  /* Main Loop (Newton iterations for finding the roots of f1 and f2) */
  while(nwt<MXNWT) {
    f1_val = eval_func(f1,x,y,new);           /* compute f1 at current (x,y) */
    f2_val = eval_func(f2,x,y,old);           /* compute f2 at current (x,y) */
    nrm2 = f1_val*f1_val + f2_val*f2_val;               /* norm of functions */
    if(nrm2<TOL2) {                                               
      /* Success! copy (x,y) into tip_x, tip_y, increment ntips, and return. */
      tip_x[ntips] = x;  
      tip_y[ntips] = y;
      ntips++;
      return;
    }
    if(nrm2>TOOBIG) {
      /* Newton iterations blowing up */
      PRNT_WARN(verbose,2,"Warning: nrm2 > TOOBIG in tip_finder() \n");
      if(!restart && tip_ini(&x,&y)) { 
	/* Try one restart */
	restart = 1;
	nwt = 0;
      }
      else {
	/* Already tried one restart, or else tip_ini() failed, so give up */
	PRNT_WARN(verbose,2,"Could not find tip. giving up \n");
	FAIL_RETURN
      }
    }
    /* Compute new x and y by taking one Newton step */
    df1dx = (eval_func(f1,x+DX,y,new) - f1_val)/DX;
    df2dx = (eval_func(f2,x+DX,y,old) - f2_val)/DX;
    df1dy = (eval_func(f1,x,y+DY,new) - f1_val)/DY;
    df2dy = (eval_func(f2,x,y+DY,old) - f2_val)/DY;
    jac   = df1dx*df2dy - df1dy*df2dx;
    x -= ( df2dy*f1_val - df1dy*f2_val)/jac;
    y -= (-df2dx*f1_val + df1dx*f2_val)/jac;
    nwt++;
  }
  /* Too many Newton iterations */
  PRNT_WARN(verbose,2,"Warning: nwt > MXNWT in tip_finder()\n");
  FAIL_RETURN
}
/*---------------------------------------------------------------------------*/

int tip_ini(precision *x, precision *y) 
     /* 
      *  Finds tip by brute force search of the grid. 
      */
{
  precision sign_1,sign_2,sign_3,sign_4; 
  int i,j;

  for(i=10;i<=ndim-10;i++) {
    for(j=10;j<=ndim-10;j++) {
      sign_1 = f1(i,j)*f1(i+1,j);
      sign_2 = f1(i,j)*f1(i,j+1);
      sign_3 = f2(i,j)*f2(i+1,j);
      sign_4 = f2(i,j)*f2(i,j+1);
      if(   ( (sign_1 < 0.) && (sign_3 < 0.) )
	 || ( (sign_1 < 0.) && (sign_4 < 0.) ) 
	 || ( (sign_2 < 0.) && (sign_3 < 0.) ) 
	 || ( (sign_2 < 0.) && (sign_4 < 0.) ) ) {
	/* found tip */
	*x = 0.5 + (precision)i;
	*y = 0.5 + (precision)j;
	return(1);
      }
    }
  }
  PRNT_WARN(verbose,2,"Warning: no tip found by tip_ini() \n");
  return(0);
}
/*---------------------------------------------------------------------------*/

precision eval_func(precision (*f)(int i, int j), precision x, precision y, 
		    int new)

/*
 *  Uses Lagrangian interpolation to find value of function f at point (x,y).
 *  If new=1, then this is a new (x,y) point and new interpolating polynomials 
 *  are computed.
 */
{
  static precision x_local, y_local, hx[NORDER+1], hy[NORDER+1]; 
  static int ix_0, iy_0;
  precision sum=0.;
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
  if(ix_0<1||ix_0>ndim-NORDER||iy_0<1||iy_0>ndim-NORDER) {
    PRNT_WARN(verbose,2,"Warning: (x,y) out of range in eval_func() \n");
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
/*---------------------------------------------------------------------------*/

precision hh(precision x_local, int i)

     /* 
      *  Computes the ith Lagrange interpolating polynomial at point x_local. 
      */
{
  precision product=1.;
  int j;

  for(j=0;j<(NORDER+1);j++) 
    if(j!=i) product *= (x_local-(precision)j)/(precision)(i-j);

  return(product);
} 
/*---------------------------------------------------------------------------*/
/*                                                                           */
/*                      Define Functions f1 and f2                           */
/*                                                                           */
/*---------------------------------------------------------------------------*/

/*  f1 and f2 define spiral tip as the intersection of contours  */

#define UCONT 0.0
#define WCONT ((c_param-sqrt(c_param*c_param-4.0*a_param*b_param))/2.0/a_param)

precision f1(int i, int j) 
{  return( (precision)(u[i][j]-UCONT) );  }


precision f2(int i, int j) 
{  return( (precision)(w[i][j]-WCONT) );  }

/*---------------------------------------------------------------------------*/
