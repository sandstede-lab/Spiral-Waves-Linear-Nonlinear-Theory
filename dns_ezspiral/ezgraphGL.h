/* ------------------------------------------------------------------------- */
/* ezgraphGL.h -- Macros for EZGRAPHGL
 *
 * Copyright (C) 1992, 1994, 1997, 1998 Dwight Barkley                 
 *
 * RCS Information
 * ---------------------------
 * $Revision: 3.2.1.1 $
 * $Date: 2007/05/07 10:00:16 $
 * ------------------------------------------------------------------------- */

#ifndef _EZGRAPH_
#define _EZGRAPH_

#define WINDOW_TITLE "EZ-Spiral"

#define WINX         0      /* Window location, in pixels, from screen left */
#define WINY         0      /* Window location, in pixels, from screen top. */
#define WM_CTRLS_POS 0      /* If WM_CTRL_POS is 0 then the window is placed
			     * at (WINX, WINY). If WM_CTRL_POS is 1 then WINX
			     * and WINY are ignored and the location is
			     * determined by the window manager. */
#define WINSIZE      400    /* Window is square of this size in pixels. */
#define PLOT_SIZE    1.2    /* This controls the size of the simulation
			     * volume in the view port: >1.0 for larger
			     * size, <1.0 for smaller size. */

#define BACKGROUND   1.0    /* Background color (R=G=B=BACKGROUND, so 0.0 gives
			       BLACK, 1.0 gives WHITE) */

#define START_PAUSED 0      /* If 1 then window is opened in paused mode
			     * showing initial condition. */

#define CLASSIC_COLORS 0    /* If 1 then use original ezspiral colors */

/* Here I plot the tip path as a line. For multiple tips you will probably
   want to use points. by defining TIP_PLOT_TYPE to be GL_POINTS  */

#define TIP_PLOT_TYPE GL_LINE_STRIP
#define TIP_WT  1.0
#define TIP_R   1.0
#define TIP_G   1.0
#define TIP_B   1.0


/* --------------------------------------------- 
 * You probably should not change anything below 
 * --------------------------------------------- */

#define PX(x) ((rect_h*((x)-1.))-half_width)
#define PY(y) ((rect_h*((y)-1.))-half_height)


#define TRUE             1   
#define FALSE            0

#define U_FIELD          0     /* These are used to determine which field */
#define V_FIELD          1     /* (if any) is being plotted */
#define NO_FIELD        -1

#define MODE_SIMULATING  1   
#define MODE_VIEWING     2   

#endif /*  _EZGRAPH_  */
