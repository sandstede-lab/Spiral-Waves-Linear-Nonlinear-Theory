/* ------------------------------------------------------------------------- 
 * ezstep.c -- Time Stepping Routines 
 *
 * Copyright (C) 1992, 1994, 1997, 1998 Dwight Barkley                 
 *
 * RCS Information
 * ---------------------------
 * $Revision: 3.2.1.1 $
 * $Date: 2007/05/07 09:50:42 $
 * ------------------------------------------------------------------------- */

#include "ezspiral.h"
#include "ezstep.h"

/* -------------------------------------------------------------------------
 * I use many pre-processor macros to make the time-stepping flexible.
 * See ezstep.h for the definitions.
 * ------------------------------------------------------------------------- */

static void  Average_periodic_directions (void);
static void  Impose_boundary_conditions  (void);
static int   s1, s2; 

const Real zero   = 0.;
const Real half   = 0.5;
const Real third  = 1./3.;
const Real one    = 1.;
const Real two    = 2.;
const Real four   = 4.;
const Real twenty = 20.;

/* ========================================================================= */

void Step (void)
{
  Real     u_thresh;
  register int i, j; 
  register int index, index_1, index_2; 
  
  /* Interchange s1 and s2 */
  s1 = 1-s1;
  s2 = 1-s2;
  
  Average_periodic_directions();

  /* ------------------------------------- 
   * Main Loop (almost all work done here) 
   * ------------------------------------- */

  for(j=1;j<=NY;j++) {
    index = INDEX(0,j); 
    index_1 = index + s1*field_size;
    index_2 = index + s2*field_size;
    for(i=1;i<=NX;i++) {
      index++; index_1++; index_2++;

#if SPLIT        /* (split=1 is best for large time steps, else use split=0) */
      u[index] += U_DIFFUSION(index_1);
      v[index] += V_DIFFUSION(index_1);
      ZERO_USED_SUMS(index_1);
#endif

      if(u[index]<delta) {  
	u[index] = zero;
	v[index] = V_KINETICS(zero,v[index]);
#if !SPLIT
	u[index] += U_DIFFUSION(index_1);
	v[index] += V_DIFFUSION(index_1);
#endif
      }
      else { 
	u_thresh = U_THRESHOLD(v[index]);
	v[index] = V_KINETICS(u[index],v[index]);
	u[index] = U_KINETICS(u[index],u_thresh);
#if !SPLIT
	u[index] += U_DIFFUSION(index_1);
	v[index] += V_DIFFUSION(index_1);
#endif
	ADD_TO_U_SUM(index,index_2);
      }
      ADD_TO_V_SUM(index,index_2);
      ZERO_USED_SUMS(index_1);
    }
  }
  Impose_boundary_conditions();
}    
/* ========================================================================= */

void Step_ini (void)
{
  int i, j;
  int index, index_1, index_2;

  /* Set initial s1 and s2 */
  s1 = 0;
  s2 = 1;

  Average_periodic_directions();

  /* Initialize spatial sums. */

  for(j=1;j<=NY;j++) {
    for(i=1;i<=NX;i++) {
      index = INDEX(i,j);
      index_1 = index + s1*field_size;
      index_2 = index + s2*field_size;
      ADD_TO_U_SUM(index,index_2);
      ADD_TO_V_SUM(index,index_2);
      ZERO_USED_SUMS(index_1);
    }
  }
  Impose_boundary_conditions();
}
/* ========================================================================= */

void Average_periodic_directions (void)
{
  /* ----------------------------------------------------------------------- 
   * In my implementation, if periodic boundary conditions are imposed in a
   * direction (x-direction, say), then the fields at grid points 1 and NX
   * should be the same.  This routine insures this by replacing these values
   * by their average.  This should be called by step_ini() to insure these
   * values are initially the same. They should thereafter remain the same
   * except for the effects of roundoff.  I have put a call to this routine
   * in step() and have put a counter here so that the averaging is done only
   * after a large number of time steps (500).  In my tests the fields differ
   * after this number of time steps by only a few times machine epsilon.
   * Note: if the averaging is done every time step then a numerical
   * instability occurs. 
   * ----------------------------------------------------------------------- */

  static int count=0;

  if(count++%500) return;

#if PBC_x
  {
    int j;
    for(j=1;j<=NY;j++) {
      U(NX,j) = U(1,j) = half*(U(NX,j)+U(1,j));
      V(NX,j) = V(1,j) = half*(V(NX,j)+V(1,j));
    }
  }
#endif
#if PBC_y
  {
    int i;
    for(i=1;i<=NX;i++) {
      U(i,NY) = U(i,1) = half*(U(i,NY)+U(i,1));
      V(i,NY) = V(i,1) = half*(V(i,NY)+V(i,1));
    }
  }
#endif
}
/* ========================================================================= */

void Impose_boundary_conditions (void)
{
  int i,j;

  /* First the values of u on the "real" part of the mesh are copied to
   * "fictitious" boundary points.  The way points are copied depend on
   * whether Neumann or periodic boundary conditions are being imposed.
   *
   * Note that the indices i, j do not range over the same values in each
   * case (x, y).  This is correct though a bit subtle. */

  /* Set fictitious points in x-direction */
  for(j=1;j<=NY;j++) { 
#if PBC_x
    U(0   ,j) = U(NX-1,j);  
    U(NX+1,j) = U(2   ,j);  
#else
    U(0   ,j) = U(2   ,j);  
    U(NX+1,j) = U(NX-1,j);  
#endif
  }

  /* Set fictitious points in y-direction */
  for(i=0;i<=NX+1;i++) { 
#if PBC_y
    U(i,0   ) = U(i,NY-1);  
    U(i,NY+1) = U(i,2   );  
#else
    U(i,0   ) = U(i,2   );  
    U(i,NY+1) = U(i,NY-1);  
#endif
  }

#if NINEPOINT

/* -----------------------------
 * Nine-point Laplacian formulas
 * ----------------------------- */

  for(i=1;i<=NX;i++) {
    Sigma_u(s2,i  ,1) += four*U(i,0);  
    Sigma_u(s2,i+1,1) +=      U(i,0);  
    Sigma_u(s2,i-1,1) +=      U(i,0); 

    Sigma_u(s2,i  ,NY) += four*U(i,NY+1); 
    Sigma_u(s2,i+1,NY) +=      U(i,NY+1);  
    Sigma_u(s2,i-1,NY) +=      U(i,NY+1); 
  } 

  for(j=1;j<=NY;j++) {
    Sigma_u(s2,1,j  ) += four*U(0,j); 
    Sigma_u(s2,1,j+1) +=      U(0,j);  
    Sigma_u(s2,1,j-1) +=      U(0,j); 

    Sigma_u(s2,NX,j  ) += four*U(NX+1,j); 
    Sigma_u(s2,NX,j+1) +=      U(NX+1,j);  
    Sigma_u(s2,NX,j-1) +=      U(NX+1,j); 
  }

  Sigma_u(s2,1 ,1 ) += U(  0 ,  0 );  
  Sigma_u(s2,NX,1 ) += U(NX+1,  0 );  
  Sigma_u(s2,1 ,NY) += U(  0 ,NY+1);  
  Sigma_u(s2,NX,NY) += U(NX+1,NY+1);  

#else

/* -----------------------------
 * Five-point Laplacian formulas 
 * ----------------------------- */

  for(i=1;i<=NX;i++) {
    Sigma_u(s2,i,1 ) += U(i,0   );  
    Sigma_u(s2,i,NY) += U(i,NY+1); 
  }

  for(j=1;j<=NY;j++) {
    Sigma_u(s2,1 ,j) += U(0   ,j); 
    Sigma_u(s2,NX,j) += U(NX+1,j); 
  }

#endif  /* end if Nine-point Laplacian formulas */


#if V_DIFF_ON  

  /* ------------------------------
   * v is diffusing.  Apply its BCS
   * ------------------------------ */

  /* Set fictitious points in x-direction */
  for(j=1;j<=NY;j++) { 
#if PBC_x
    V(0   ,j) = V(NX-1,j);  
    V(NX+1,j) = V(2   ,j);  
#else
    V(0   ,j) = V(2   ,j);  
    V(NX+1,j) = V(NX-1,j);  
#endif
  }

  /* Set fictitious points in y-direction */
  for(i=0;i<=NX+1;i++) { 
#if PBC_y
    V(i,0   ) = V(i,NY-1);  
    V(i,NY+1) = V(i,2   );  
#else
    V(i,0   ) = V(i,2   );  
    V(i,NY+1) = V(i,NY-1);  
#endif
  }

#if NINEPOINT

/* -----------------------------
 * Nine-point Laplacian formulas
 * ----------------------------- */

  for(i=1;i<=NX;i++) {
    Sigma_v(s2,i  ,1) += four*V(i,0);  
    Sigma_v(s2,i+1,1) +=      V(i,0);  
    Sigma_v(s2,i-1,1) +=      V(i,0); 

    Sigma_v(s2,i  ,NY) += four*V(i,NY+1); 
    Sigma_v(s2,i+1,NY) +=      V(i,NY+1);  
    Sigma_v(s2,i-1,NY) +=      V(i,NY+1); 
  } 

  for(j=1;j<=NY;j++) {
    Sigma_v(s2,1,j  ) += four*V(0,j); 
    Sigma_v(s2,1,j+1) +=      V(0,j);  
    Sigma_v(s2,1,j-1) +=      V(0,j); 

    Sigma_v(s2,NX,j  ) += four*V(NX+1,j); 
    Sigma_v(s2,NX,j+1) +=      V(NX+1,j);  
    Sigma_v(s2,NX,j-1) +=      V(NX+1,j); 
  }

  Sigma_v(s2,1 ,1 ) += V(  0 ,  0 );  
  Sigma_v(s2,NX,1 ) += V(NX+1,  0 );  
  Sigma_v(s2,1 ,NY) += V(  0 ,NY+1);  
  Sigma_v(s2,NX,NY) += V(NX+1,NY+1);  

#else

/* -----------------------------
 * Five-point Laplacian formulas 
 * ----------------------------- */

  for(i=1;i<=NX;i++) {
    Sigma_v(s2,i,1 ) += V(i,0   );  
    Sigma_v(s2,i,NY) += V(i,NY+1); 
  }

  for(j=1;j<=NY;j++) {
    Sigma_v(s2,1 ,j) += V(0   ,j); 
    Sigma_v(s2,NX,j) += V(NX+1,j); 
  }

#endif  /* end if Nine-point Laplacian formulas */

#endif  /* end if V_DIFF_ON */

}
/* ========================================================================= */

