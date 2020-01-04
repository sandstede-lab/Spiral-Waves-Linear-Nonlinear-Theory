/*
 *  Define the macros U_THRESHOLD(v) and G(u,v) according to the model you 
 *  wish to simulate. In principle these can be almost anything you want.  
 *  Three examples are included.
 * 
 * $Revision: 2.3 $
 * $Date: 1997/11/13 23:17:22 $
 *---------------------------------------------------------------------------*/

#if 1  /* Roessler */
#  define F(u,v,w)  ( -(v)-(w) )
#  define G(u,v,w)  ( (u)+a_param*(v) )
#  define H(u,v,w)  ( b_param+(w)*(u)-c_param*(w) )
#endif

/*---------------------------------------------------------------------------*/
/*          In general you should not need to change anything below          */
/*---------------------------------------------------------------------------*/

/*  I rely on the compiler (gcc) to take care of most optimization.  
 *  The only place where I give the compiler help is with the nine-point 
 *  Laplacian formulas because apparently gcc does not do the obvious 
 *  thing.  I do not unroll the main loops in step() nor do I write explicit   
 *  pointer operations because gcc seems to handle all the indirections very 
 *  efficiently.
 */

/*---------------------------------------------------------------------------*/

#define U_KINETICS(u,v,w) ( (u)+dt*F(u,v,w) )
#define V_KINETICS(u,v,w) ( (v)+dt*G(u,v,w) )
#define W_KINETICS(u,v,w) ( (w)+dt*H(u,v,w) )

/*---------------------------------------------------------------------------*/

#define U_DIFFUSION ( NDu * u_lap[k1][i][j] ) 
#define V_DIFFUSION ( NDv * v_lap[k1][i][j] )
#define W_DIFFUSION ( NDv * w_lap[k1][i][j] )


#if NINEPOINT
/* Nine-point Laplacian formulas */
#  define ADD_TO_U_LAPLACIAN           \
     {precision stupid_cc = four*u[i][j]; \
      u_lap[k2][i][j] -= twenty*u[i][j];  \
      u_lap[k2][i+1][j] += stupid_cc;  \
      u_lap[k2][i-1][j] += stupid_cc;  \
      u_lap[k2][i][j+1] += stupid_cc;  \
      u_lap[k2][i][j-1] += stupid_cc;  \
      u_lap[k2][i+1][j+1] += u[i][j];  \
      u_lap[k2][i-1][j+1] += u[i][j];  \
      u_lap[k2][i+1][j-1] += u[i][j];  \
      u_lap[k2][i-1][j-1] += u[i][j];  \
     }
#  define ADD_TO_V_LAPLACIAN            \
     {precision stupid_cc = four*v[i][j]; \
      v_lap[k2][i][j] -= twenty*v[i][j];  \
      v_lap[k2][i+1][j] += stupid_cc;   \
      v_lap[k2][i-1][j] += stupid_cc;   \
      v_lap[k2][i][j+1] += stupid_cc;   \
      v_lap[k2][i][j-1] += stupid_cc;   \
      v_lap[k2][i+1][j+1] += v[i][j];   \
      v_lap[k2][i-1][j+1] += v[i][j];   \
      v_lap[k2][i+1][j-1] += v[i][j];   \
      v_lap[k2][i-1][j-1] += v[i][j];   \
     }
#  define ADD_TO_W_LAPLACIAN            \
     {precision stupid_cc = four*w[i][j]; \
      w_lap[k2][i][j] -= twenty*w[i][j];  \
      w_lap[k2][i+1][j] += stupid_cc;   \
      w_lap[k2][i-1][j] += stupid_cc;   \
      w_lap[k2][i][j+1] += stupid_cc;   \
      w_lap[k2][i][j-1] += stupid_cc;   \
      w_lap[k2][i+1][j+1] += w[i][j];   \
      w_lap[k2][i-1][j+1] += w[i][j];   \
      w_lap[k2][i+1][j-1] += w[i][j];   \
      w_lap[k2][i-1][j-1] += w[i][j];   \
     }
#else
/* Five-point Laplacian formulas */
#  define ADD_TO_U_LAPLACIAN       \
     u_lap[k2][i][j] -= four*u[i][j]; \
     u_lap[k2][i+1][j] += u[i][j]; \
     u_lap[k2][i-1][j] += u[i][j]; \
     u_lap[k2][i][j+1] += u[i][j]; \
     u_lap[k2][i][j-1] += u[i][j];  
#  define ADD_TO_V_LAPLACIAN       \
     v_lap[k2][i][j] -= four*v[i][j]; \
     v_lap[k2][i+1][j] += v[i][j]; \
     v_lap[k2][i-1][j] += v[i][j]; \
     v_lap[k2][i][j+1] += v[i][j]; \
     v_lap[k2][i][j-1] += v[i][j]; 
#  define ADD_TO_W_LAPLACIAN       \
     w_lap[k2][i][j] -= four*w[i][j]; \
     w_lap[k2][i+1][j] += w[i][j]; \
     w_lap[k2][i-1][j] += w[i][j]; \
     w_lap[k2][i][j+1] += w[i][j]; \
     w_lap[k2][i][j-1] += w[i][j]; 
#endif

#define ZERO_USED_LAPLACIANS \
    u_lap[k1][i][j] = zero;    \
    v_lap[k1][i][j] = zero;    \
    w_lap[k1][i][j] = zero; 

/*---------------------------------------------------------------------------*/
