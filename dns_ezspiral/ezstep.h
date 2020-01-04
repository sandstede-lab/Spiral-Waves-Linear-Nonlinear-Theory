/* ------------------------------------------------------------------------- 
 * ezstep.h -- Macros for EZSTEP
 *
 * Copyright (C) 1992, 1994, 1997, 1998 Dwight Barkley                 
 *
 * RCS Information
 * ---------------------------
 * $Revision: 3.2.1.1 $
 * $Date: 2007/05/07 09:50:42 $
 * ------------------------------------------------------------------------- */

#include "ezspiral.h"
#ifndef _EZSTEP_
#define _EZSTEP_

#include "math.h"  /* added Bjorn Sandstede for Karma */

/* ---------------------------------------------------------------------
 * Define the macros U_THRESHOLD(v) and G(u,v) according to the model
 * you wish to simulate. In principle these can be almost anything you
 * want.  Three examples are included.  
 * --------------------------------------------------------------------- */

#if 1  /* Standard Model */
  #define U_THRESHOLD(v)  ( one_o_a*(v) + b_o_a )
  #define G(u,v)          ( (u)-(v) )
#endif

#if 0  /* Simple Model for Turbulent Spirals */
  #define U_THRESHOLD(v)  ( one_o_a*(v) + b_o_a )
  #define G(u,v)          ( (u)*(u)*(u) - (v) )
#endif

#if 0  /* Bar and Eiswirth Model (Phys. Rev. E V. 48, p. 1635, 1993).
	* I do not include their case u>1 because it is not relevant. */
  #define U_THRESHOLD(v)  ( one_o_a*(v) + b_o_a )
  #define G(u,v)          (-(v) + ( ((u)<third) ? 0 : \
                                    (one-6.75*(u)*(1-(u))*(1-(u))) ) )  
#endif

#if 0  /* Karma Model */
  #define U_THRESHOLD(v)  ( 1.5414-(v)*(v)*(v)*(v) )
  #define G(u,v)          ( 4.0*((1+tanh(4.0*(u)-4.0))*0.5/(1-exp(-a_param))-(v)) )
#endif

/* ------------------------------------------------------------------------- 
 *          In general you should not need to change anything below          
 * ------------------------------------------------------------------------- */

/* I rely on the compiler to take care of most optimization.  The only place
 * where I give the compiler help is with the Laplacian formulas because
 * none of the compiles do the obvious thing.  */

/* ------------------------------------- 
 * Defining the kinetics macros:         
 * U_KINETICS(u,uth) and V_KINETICS(u,v) 
 * ------------------------------------- */

#if EXPLICIT   
     /* Explicit u-kinetics */
  #define U_KINETICS(u,uth) ( (u)+dt_o_eps*(u)*(one-(u))*((u)-(uth)) )

     /* added Bjorn Sandstede: FHN right-hand side */
/* #define U_KINETICS(u,uth) ( (u)+dt_o_eps*((u)*((u)-half)*(one-(u))-(uth)) ) */

     /* added Bjorn Sandstede: Karma right-hand side */
/* #define U_KINETICS(u,uth) ( (u)+dt_o_eps*(-(u)+uth*(one-tanh((u)-3.0))*(u)*(u)*half) ) */
   
#else          
     /* Implicit u-kinetics */
     /* The original (Physica 49D) scheme can be obtained by defining
      * F2m and F2p to be one */

  #define F1m(u,uth) ( dt_o_eps*(one-(u))*((u)-(uth)) )
  #define F1p(u,uth) ( dt_o_eps*(u)*((u)-(uth)) )
  #define F2m(u,uth) ( one + dt_o_2eps*((uth)-(u)*(u)) )
  #define F2p(u,uth) ( one + dt_o_2eps*(two*(u) -(uth)-(u)*(u)) )

  #define U_KINETICS(u,uth) (                                      \
          ( (u) < (uth) ) ?                                        \
          (u)/(one-F1m(u,uth)*F2m(u,uth) ) :                       \
          ((u)+F1p(u,uth)*F2p(u,uth))/(one+F1p(u,uth)*F2p(u,uth)) )

#endif

#define V_KINETICS(u,v) ( (v)+dt*G(u,v) )

/* ---------------------------------------- 
 * Defining the diffusion macros:           
 * U_DIFFUSION, V_DIFFUSION, ZERO_USED_SUMS 
 * ---------------------------------------- */

#if V_DIFF_ON        
     /* v is diffusing */
  #define U_DIFFUSION(index_1)     (dt_o_wh2   * sigma_u[index_1]) 
  #define V_DIFFUSION(index_1)     (dtDv_o_wh2 * sigma_v[index_1]) 
  #define ZERO_USED_SUMS(index_1)  sigma_u[index_1] = sigma_v[index_1] = zero;
#else
     /* v is not diffusing */
  #define U_DIFFUSION(index_1)     (dt_o_wh2 * sigma_u[index_1]) 
  #define V_DIFFUSION(index_1)     zero
  #define ZERO_USED_SUMS(index_1)  sigma_u[index_1] = zero;
#endif

/* -------------------------------- 
 * Defining the spatial-sum macros: 
 * ADD_TO_U_SUM, ADD_TO_V_SUM     
 * -------------------------------- */

#if NINEPOINT
/* 9-point Laplacian formula */
  #define ADD_TO_U_SUM(index,index_2) \
   {Real stupid_cc = u[index]; \
     sigma_u[(index_2)]  -= twenty*stupid_cc; \
     sigma_u[(index_2)+J_INC]  += four*stupid_cc; \
     sigma_u[(index_2)-J_INC]  += four*stupid_cc; \
     sigma_u[(index_2)+I_INC]  += four*stupid_cc; \
     sigma_u[(index_2)-I_INC]  += four*stupid_cc; \
     sigma_u[(index_2)+I_INC+J_INC] += stupid_cc; \
     sigma_u[(index_2)+I_INC-J_INC] += stupid_cc; \
     sigma_u[(index_2)-I_INC+J_INC] += stupid_cc; \
     sigma_u[(index_2)-I_INC-J_INC] += stupid_cc; \
   }
#else
/* 5-point Laplacian formula */
  #define ADD_TO_U_SUM(index,index_2) \
   {Real stupid_cc = u[index]; \
      sigma_u[(index_2)]  -= four*stupid_cc; \
      sigma_u[(index_2)+J_INC] += stupid_cc; \
      sigma_u[(index_2)-J_INC] += stupid_cc; \
      sigma_u[(index_2)+I_INC] += stupid_cc; \
      sigma_u[(index_2)-I_INC] += stupid_cc; \
   }
#endif         

#if !V_DIFF_ON   
/* v is not diffusing */
  #define ADD_TO_V_SUM(index,index_2) 
#else            
/* v is diffusing */
#if NINEPOINT
/* 9-point Laplacian formula */
  #define ADD_TO_V_SUM(index,index_2) \
   {Real stupid_cc = v[index]; \
     sigma_v[(index_2)]  -= twenty*stupid_cc; \
     sigma_v[(index_2)+J_INC]  += four*stupid_cc; \
     sigma_v[(index_2)-J_INC]  += four*stupid_cc; \
     sigma_v[(index_2)+I_INC]  += four*stupid_cc; \
     sigma_v[(index_2)-I_INC]  += four*stupid_cc; \
     sigma_v[(index_2)+I_INC+J_INC] += stupid_cc; \
     sigma_v[(index_2)+I_INC-J_INC] += stupid_cc; \
     sigma_v[(index_2)-I_INC+J_INC] += stupid_cc; \
     sigma_v[(index_2)-I_INC-J_INC] += stupid_cc; \
   }
#else          
/* 5-point Laplacian formula */
  #define ADD_TO_V_SUM(index,index_2) \
   {Real stupid_cc = v[index]; \
      sigma_v[(index_2)]  -= four*stupid_cc; \
      sigma_v[(index_2)+J_INC] += stupid_cc; \
      sigma_v[(index_2)-J_INC] += stupid_cc; \
      sigma_v[(index_2)+I_INC] += stupid_cc; \
      sigma_v[(index_2)-I_INC] += stupid_cc; \
  }
#endif 
#endif
#endif /*  _EZSTEP_  */
