/*
 *  The parameters below affect the compilation of EZ-Spiral.
 *  Define them according to your needs.  
 */

/*---------------------------------------------------------------------------*/
#define precision float  /* Precision of arithmetic (float or double)        */
                         /*                                                  */
#define V_DIFF_ON 1      /* if 1 then v-field diffuses                       */
#define NINEPOINT 1      /* if 1 then 9 point Laplacian formulas             */
#define NBC       1      /* if 1 then Neumann BCs, if 0 periodic BCs         */
#define GRAPHICS  1      /* if 1 then compile X11 graphics                   */
/*---------------------------------------------------------------------------*/

#define min(a,b) ((a)>(b) ? (b) : (a))
#define max(a,b) ((a)>(b) ? (a) : (b))

#define PRNT_WARN(verbose,level,message) if(verbose>=level) printf(message);

/* Include Standard C Headers */
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

/* Prototypes */
void     step          ();
void     step_ini      ();
void     plot          ();
void     plot_ini      ();
int      keyboard_chk  ();
void     tip_finder    ();
void     save_image    (int num);
