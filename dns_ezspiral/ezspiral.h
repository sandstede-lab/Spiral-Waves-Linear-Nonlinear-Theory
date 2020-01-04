/* ------------------------------------------------------------------------- *
 * ezspiral.h -- header file for EZ-Spiral
 *
 * Copyright (C) 1992, 1994, 1997, 1998, 2007 Dwight Barkley                 
 *
 * RCS Information
 * ---------------------------
 * $Revision: 3.2.1.1 $
 * $Date: 2007/05/07 09:50:42 $
 * ------------------------------------------------------------------------- */

#ifndef _EZSPIRAL_
#define _EZSPIRAL_

/* ------------------------------------------------------------------------- */
/* These are the main parameters affecting the compilation of EZ-Spiral.     */
/* Define them according to your needs.                                      */
/*                                                                           */
typedef float  Real;      /* precision of Real variables (float or double)   */
                          /*                                                 */
#define V_DIFF_ON  1      /* if 1 then v-field diffuses                      */
#define EXPLICIT   1      /* if 1 then explicit u-kinetics, else implicit    */
#define SPLIT      1      /* if 1 then split diffusion and kinetics steps    */
#define NINEPOINT  1      /* if 1 then 9 pt Laplacian formulas, else 5 pt    */
#define PBC_x      0      /* if 1 then periodic BCs in x, else Neumann BCs   */
#define PBC_y      0      /* if 1 then periodic BCs in y, else Neumann BCs   */
#define GRAPHICS   1      /* if 1 then run with interactive graphics         */
/* ------------------------------------------------------------------------- */

/* 
 * I always define min, max, and make sure M_PI is defined 
 * ------------------------------------------------------- */
#define min(a,b)      ((a)>(b) ? (b) : (a))   
#define max(a,b)      ((a)>(b) ? (a) : (b))   
#ifndef M_PI
#define M_PI	      3.14159265358979323846
#endif

/* --------------------------------------------------- 
 * Global variables used throughout the EZ-Spiral Code 
 * (I use lots of them)
 * --------------------------------------------------- */

extern Real *fields, *u, *v,             /* arrays for concentration fields */
            *sigma_u, *sigma_v;          /* and spatial sums (Laplacians) */
extern Real *tip_x, *tip_y;              /* spiral tip arrays */
extern Real a_param, b_param, one_o_eps, /* model parameters */
            delta, grid_h, dt,           /* numerical parameters */
            dt_o_eps, dt_o_2eps,         /* useful parameter combinations */
            one_o_a, b_o_a,              /* useful parameter combinations */
            dt_o_wh2, dtDv_o_wh2;        /* useful parameter combinations */
extern int  nx, ny, nz,                  /* # of grid points per direction */
            field_size,                  /* array size for each field */
            nsteps,                      /* # time steps to take */
            istep,                       /* current time step */
            write_tip,                   /* write tip flag */
            simulating_resolution,       /* graphics parameter */
            verbose,                     /* verbosity level */
            ntips;                       /* number of spiral tips */

/* -------------------------------------------------------------------------
 * All index ranges throughout the code are expressed in terms of NX, and 
 * NY.  The code is generally more efficient if these are known numbers
 * at compile time, but then one must recompile for each change.  If the
 * values are known (eg specified on the compile line) then do nothing
 * here, otherwise define NX etc to be the *variables* nx, etc.
 * ------------------------------------------------------------------------- */

#ifndef NX
  #define NX  nx
#endif

#ifndef NY
  #define NY  ny
#endif

/* ------------------------------------------------------------------------- 
 * Memory for the chemical fields (u and v) and the spatial sums (sigma_u    
 * and sigma_v) is allocated in Allocate_memory() in ezspiral.c.  These are      
 * allocated as long (single dimensional) arrays.  Here macros are defined   
 * so that one can easily reference a value corresponding to a particular    
 * grid point, i.e. macros are defined to treat all arrays as                
 * multi-dimensional. If you don't like the way I do this, you should be     
 * able to change it easily by making modifications here and in              
 * AllocateMem().  Let me know if you find a significant improvement.        
 *                                                                           
 * INDEX(i,j) converts grid point (i,j) to the array index.              
 *                                                                           
 * U(i,j)          --  u field at grid point (i,j)
 * V(i,j)          --  v field at grid point (i,j)
 * Sigma_u(s,i,j)  --  Spatial-sum array for u: s=0 or 1 (see ezstep.c)    
 * Sigma_v(s,i,j)  --  Spatial-sum array for v: s=0 or 1 (see ezstep.c)    
 * Fields (f,i,j)  --  array of fields: f=0 for u or f=1 for v.            
 * ------------------------------------------------------------------------- */

#define J_INC       (NX+2)
#define I_INC        1
#define FIELD_SIZE ((NY+2)*(NX+2)) 

#define INDEX(i,j)  ((i)*I_INC + (j)*J_INC)

#define U(i,j)          u[INDEX(i,j)]
#define V(i,j)          v[INDEX(i,j)]
#define Sigma_u(s,i,j)  sigma_u[(s)*FIELD_SIZE + INDEX(i,j)]
#define Sigma_v(s,i,j)  sigma_v[(s)*FIELD_SIZE + INDEX(i,j)]
#define Fields(f,i,j)   fields [(f)*FIELD_SIZE + INDEX(i,j)]


/* ------------------------------------------- 
 * Prototypes for public functions defined in: 
 * ------------------------------------------- */

/* ezspiral.c 
 * ---------- */
void Write_tip_data (float x, float y);

/* ezstep3d.c 
 * ---------- */
void  Step      (void);
void  Step_ini  (void);

/* ezgraph3d.c 
 * ----------- */
void  Draw         (void);
void  Draw_ini     (int initial_field);
int   Event_check  (void);
void  QuitX        (void);

/* eztip.c 
 * ------- */
void  Find_tips  (void);
void  Tip_ini    (void);

#endif /*  _EZSPIRAL_  */

