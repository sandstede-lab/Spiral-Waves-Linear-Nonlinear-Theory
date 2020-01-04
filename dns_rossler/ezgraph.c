/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*                                                                           */
/*                                 EZGRAPH.C                                 */
/*                   X11-Graphics Subroutines for EZ-SPIRAL                  */
/*                                                                           */
/*                                Version 2                                  */
/*                                                                           */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* 
 * $Revision: 2.2 $
 * $Date: 1994/06/27 23:31:12 $
 */
/*---------------------------------------------------------------------------*/

#include "ezspiral.h"

#if GRAPHICS   

#include "ezgraph.h"

/* External variables */
extern precision  **u, **v, **w, **u_lap[2],**v_lap[2], **w_lap[2],
                  *tip_x, *tip_y;
extern int        ndim, k2, ntips, u_plot, v_plot, verbose;

/*---------------------------------------------------------------------------*/

void plot()

     /* 
      *  Main plotting routine. 
      *  Filled rectangles of appropriate colors are drawn to pixmaps.
      *  The spiral tip path is also drawn to the pixmaps if desired.
      *  After all graphics are drawn, the pixmaps are copied to the 
      *  X-window.  This is similar to double buffering and gives better
      *  animation than drawing directly to the X-window.
      */
{
  int i,j;

  if(!(u_plot||v_plot)) return;  /* neither field is plotted so return */

  switch (plot_type) {
     case Color :  break;
     case Mono :   return;
  }

  if(u_plot) {
    for(i=1;i<=ndim;i++) {
      for(j=1;j<=ndim;j++) {
	setColor(u[i][j],'u');
	XFillRectangle(theDisplay, theUPixmap, theGC, 
			 PX(i), PY(j), irsize, irsize);
      }
    }
  }

  if(u_plot) {
    for(i=1;i<=ndim;i++) {
      for(j=1;j<=ndim;j++) {
	setColor(v[i][j],'v');
	XFillRectangle(theDisplay, theVPixmap, theGC, 
		       PX(i), PY(j), irsize, irsize);
      }
    }
  }

  if(v_plot) {
    for(i=1;i<=ndim;i++) {
      for(j=1;j<=ndim;j++) {
	setColor(w[i][j],'w');
	XFillRectangle(theDisplay, theWPixmap, theGC, 
		       PX(i), PY(j), irsize, irsize);
      }
    }
  }

  /* If tip path is to be plotted, then do it here. */

  if(ntips) plot_tip_path();

  /*  Copy pixmaps to the window so that they become visible.
   *  The u_origin and v_origin determine where in the window the
   *  pixmaps are copied.  See plot_ini().
   */

  if(u_plot) { XCopyArea(theDisplay, theUPixmap, theWindow, theGC, 
	                 0, 0, sqsize, sqsize, u_origin[0], u_origin[1]); }
  if(u_plot) { XCopyArea(theDisplay, theVPixmap, theWindow, theGC, 
	                 0, 0, sqsize, sqsize, v_origin[0], v_origin[1]); }
  if(v_plot) { XCopyArea(theDisplay, theWPixmap, theWindow, theGC, 
	                 0, 0, sqsize, sqsize, w_origin[0], w_origin[1]); }

  XFlush(theDisplay);
}

/*---------------------------------------------------------------------------*/

void save_image(int iter)

  /* Saves window image to file */
{
  static char filename[80];
  sprintf(filename, "spiral%4d.xpm", iter);
  XpmWriteFileFromPixmap(theDisplay,filename,theWPixmap,NULL,NULL);
}

/*---------------------------------------------------------------------------*/

int setColor (precision val, char cval)

     /* 
      *  Sets the foreground color.
      *  val is the value of the field of "type" cval (cval = 'u' or 'v').
      *  Returns 1 if cval == 'u' and the state is non-quiescent.
      */
{
  float scaled_u, scaled_v, scaled_w;
  int icol;

#define UMAX  5.0
#define UMIN -5.0
#define VMAX  5.0
#define VMIN -5.0

/*
#define WMAX  4.0
#define WMIN  0.0
*/
#define WMAX  4.0
#define WMIN  0.0

    scaled_u = (val-UMIN)/(UMAX-UMIN);
    scaled_u = max(0.0,  scaled_u);
    scaled_u = min(0.999,scaled_u);
    icol = (int)maxColors*scaled_u;
    XSetForeground(theDisplay, theGC, theColors[icol] );

  switch (cval) {
  case 'u' :
    scaled_u = (val-UMIN)/(UMAX-UMIN);
    scaled_u = max(0.0,  scaled_u);
    scaled_u = min(0.999,scaled_u);
    icol = (int)maxColors*scaled_u;
    XSetForeground(theDisplay, theGC, theColors[icol] );
    return(1);
  case 'v' :
    scaled_v = (val-VMIN)/(VMAX-VMIN);
    scaled_v = max(0.0,  scaled_v);
    scaled_v = min(0.999,scaled_v);
    icol = (int)maxColors*scaled_v;
    XSetForeground(theDisplay, theGC, theColors[icol] );
    return(1);
  case 'w' :
    scaled_w = (val-WMIN)/(WMAX-WMIN);
    scaled_w = max(0.0,  scaled_w);
    scaled_w = min(0.999,scaled_w);
    icol = (int)maxColors*scaled_w;
    XSetForeground(theDisplay, theGC, theColors[icol] );
    return(1);
  }

}
/*---------------------------------------------------------------------------*/

void plot_tip_path()

     /*  Plots the stored tip path  */
{
  XPoint *thePoints;
  int i;
  
  if(ntips<1) return;

  XSetForeground(theDisplay, theGC, theWhitePixel);
  thePoints=(XPoint *)malloc((unsigned)(ntips)*sizeof(XPoint));

  for(i=0;i<ntips;i++) {
    thePoints[i].x = (short) PX(tip_x[ntips-i-1]);
    thePoints[i].y = (short) PY(tip_y[ntips-i-1]);
  }

  if(u_plot) { XDrawLines(theDisplay, theUPixmap, theGC,
                          thePoints, ntips, CoordModeOrigin);  }
  if(u_plot) { XDrawLines(theDisplay, theVPixmap, theGC, 
	                  thePoints, ntips, CoordModeOrigin);  }
  if(v_plot) { XDrawLines(theDisplay, theWPixmap, theGC, 
	                  thePoints, ntips, CoordModeOrigin);  }

  free(thePoints);
}
/*---------------------------------------------------------------------------*/

void plot_ini()

     /* 
      *  Initializes the graphics. 
      */
{
  int width, height;

  if(!(u_plot||v_plot)) return;  /* neither field is plotted so return */

  /* Initialize X */
  initX();

  /*  Set up the graphics coordinates.
   *  sqsize is the length, in pixels, of one side of the square area in 
   *  which the u,v fields are plotted.  PLOTSIZE is the fraction of the
   *  full screen height for one of the square plots and thePHeight the 
   *  screen height in pixels.
   *  
   *  width and height are the width and height of the X-window in pixels.
   *  The height is always sqsize + 2*BORDER (to account for BORDER pixels
   *  above and BORDER pixels below the square).
   *  In the case of both fields plotted, the width is 2*sqsize (for the
   *  2 fields) + 3*BORDER (right BORDER, left BORDER and a BORDER 
   *  between the two plots).
   *  In the case of one field plotted, the width is sqsize + 2*BORDER 
   *  exactly as for the height.
   *
   *  The arrays u_origin and v_origin control where the u and v squares 
   *  appear within the X-window.  See also plot().  
   *  In the case of both fields, the upper-left corner of the u-field 
   *  is copied to the window BORDER pixels down and BORDER pixels to 
   *  the right of the upper-left corner of the window.  The upper-left 
   *  corner of the v-field is copied to the window BORDER pixels down and 
   *  sqsize+2*BORDER pixels to the right of the upper-left corner of the 
   *  window.  
   *  In the case of one field, the upper-left corner of that field
   *  is copied to the window BORDER pixels down and BORDER pixels to 
   *  the right of the upper-left corner of the window.  Both u_origin
   *  and v_origin are set but only one is used.
   */

  /*  sqsize = PLOTSIZE*thePHeight;  */
  sqsize = PLOTSIZE;

  if(u_plot && v_plot) {
    /*  (u,v,w) are plotted */
    height =   sqsize + 2*BORDER;
    width  = 3*sqsize + 4*BORDER;
    u_origin[0] = BORDER;
    u_origin[1] = BORDER;
    v_origin[0] = sqsize + 2*BORDER;
    v_origin[1] = BORDER;
    w_origin[0] = 2*sqsize + 3*BORDER;
    w_origin[1] = BORDER;
  }
  else {
    /*  Only one field plotted */
    height = width = sqsize + 2*BORDER;
    w_origin[0] = BORDER;
    w_origin[1] = BORDER;
  }

  /* rsize = the "exact" (i.e. floating point) size of the rectangles drawn 
   * at each grid point.  irsize is the integer equivalent, with the +1 
   * added to ensure that rectangles are overlapping.
   */
  rsize  = (float)sqsize/(float)ndim;
  irsize = (int)(rsize)+1;

  /* Open a window and call other necessary X routines */
  openWindow(WINX, WINY, width, height);

  /* If Color plot, then initialize colors */
  if(plot_type==Color) initColors();
}
/*---------------------------------------------------------------------------*/

void initX()

     /* 
      *  Initializes X.
      */
{
  char *theDisplayName = NULL;
  
  /* Open the display */
  if( (theDisplay = XOpenDisplay(NULL)) == NULL) {
    fprintf(stderr, 
	    "ERROR: Could not open a connection to X on display %s\n",
	    XDisplayName(theDisplayName));
    exit(1);
  }

  theScreen     = DefaultScreen(theDisplay);
  theDepth      = DefaultDepth (theDisplay, theScreen);
  thePWidth     = DisplayWidth (theDisplay, theScreen);
  thePHeight    = DisplayHeight(theDisplay, theScreen);
  theBlackPixel = BlackPixel   (theDisplay, theScreen);
  theWhitePixel = WhitePixel   (theDisplay, theScreen);
  theColormap   = DefaultColormap(theDisplay, theScreen);

  if(theDepth==1) {
    plot_type = Mono;
  }
  else {
    plot_type = Color;
  }

  /* Print useful information. I suggest printing this at least once */
  if(verbose>1) {

    if(theDepth==1) {
      printf("Color plane depth...........%d (monochrome)\n", theDepth);
    }
    else {
      printf("Color plane depth...........%d \n",             theDepth);
    }

    printf("Display Width...............%d \n", thePWidth);
    printf("Display Height..............%d \n", thePHeight);
    printf("The display %s\n\n", XDisplayName(theDisplayName));
  }
}
/*---------------------------------------------------------------------------*/

void openWindow(int winx, int winy, int width, int height)

     /* 
      *  Creates and maps the X window.
      *  Also creates the Graphics Context and Pixmaps.
      */
{
  XSizeHints    theSizeHints;
  XTextProperty theWindowName, theIconName;
  char          *Name="EZ-Spiral";
  unsigned long background;

  /* Set the background color (i.e. color of border around the square plots */
  background = BlackPixel(theDisplay,theScreen);
  if(BACKGROUND) background = WhitePixel(theDisplay,theScreen);

  /* Create the window */
  theWindow = XCreateSimpleWindow(theDisplay,
				  RootWindow(theDisplay, theScreen),
				  winx, winy, width, height, 0,
				  background, background);  

  /* Set standard window properties. 
   * theWindowName and theIconName are set to Name, which unless you
   * change it, it will be EZ-Spiral.
   * theSizeHints.flags are set such that if winx and winy both positive, 
   * they try to locate the window there, else let the window manager decide.
   */
  XStringListToTextProperty(&Name,1,&theWindowName);
  XStringListToTextProperty(&Name,1,&theIconName);

  if(winx>0 && winy>0) {
    theSizeHints.flags  = USPosition | PSize;
  } else {
    theSizeHints.flags  = PPosition | PSize;
  }

  XSetWMProperties(theDisplay, theWindow, &theWindowName, &theIconName,
		   NULL, 0, &theSizeHints, NULL, NULL);

  /* This selects which events the program will see.  This is used, e.g. to 
   * pass keyboard or mouse input to the program.  The event types are set 
   * in EV_MASK.  (See ezgraph.h)
   */
  XSelectInput(theDisplay, theWindow, EV_MASK);

  /* Create the Graphics Context. */
  theGC = XCreateGC(theDisplay, theWindow, (unsigned long) 0, NULL);

  if(theGC==0) {
    XDestroyWindow(theDisplay,theWindow);
    printf("Could not create the GC \n");
    exit(1);
  }

  XSetGraphicsExposures(theDisplay,theGC,0);

  /* Now that the Graphics Context has been created, setting and changing
   * most graphics characteristics (such as line width or line type) can 
   * easily be done with appropriate X routine calls.  I include an example 
   * of how to set the line width to 2 pixels.  LineSolid has the obvious
   * meaning.  CapRound and JoinRound are probably not important to you.
   *
   *  XSetLineAttributes(theDisplay, theGC, (unsigned int) 2, 
   *                     LineSolid, CapRound, JoinRound);
   */

  /* Map the window, i.e. make it visible */
  XMapWindow(theDisplay, theWindow); 
  XFlush(theDisplay);

  /* Create pixmaps of size sqsize by sqsize for the u/v plots */
  if(u_plot) {
    theUPixmap = XCreatePixmap(theDisplay,theWindow,
			       sqsize, sqsize, theDepth);
    if(theUPixmap==False) printf("couldn't open UPixmap\n");
  }

  if(u_plot) {
    theVPixmap = XCreatePixmap(theDisplay,theWindow,
			       sqsize, sqsize, theDepth);
    if(theVPixmap==False) printf("couldn't open VPixmap\n");
  }

  if(v_plot) {
    theWPixmap = XCreatePixmap(theDisplay,theWindow,
			       sqsize, sqsize, theDepth);
    if(theWPixmap==False) printf("couldn't open WPixmap\n");
  }
}
/*---------------------------------------------------------------------------*/

void initColors()

{
  XColor theRGBColor, theHardwareColor;
  int    theStatus;
  int    cv,i;

  /* Colors for the u-field: Red, Blue, (and black already defined) */
  
  if( XAllocNamedColor(theDisplay, theColormap, "Red",
		       &theRGBColor, &theHardwareColor) ) {
    theRedPixel = theHardwareColor.pixel; 
  }
  else {
    printf("Cannot allocate the Red Pixel \n");
  }
  
  if( XAllocNamedColor(theDisplay, theColormap, "Blue",
		       &theRGBColor, &theHardwareColor) ) {
    theBluePixel = theHardwareColor.pixel; 
  }
  else {
    printf("Cannot allocate the Blue Pixel \n");
  }
  
  /*  Colors for the v-field: Continuous range from Red to Blue.
   *  A few possibilities are included; you should modify according
   *  to taste.
   */
  
  for(i=0;i<maxColors;i++) {
    
    cv = (int)( i*65535./(maxColors-1.) );
#if 0   /* blue-red */
    theRGBColor.red   = theHardwareColor.red   = 65535 - cv;
    theRGBColor.green = theHardwareColor.green = cv/1.6; 
    theRGBColor.blue  = theHardwareColor.blue  = cv;
#endif

#if 0  /* Barkley's alternative color map */
    theRGBColor.red   = theHardwareColor.red   = 65535 - cv;
    theRGBColor.green = theHardwareColor.green = cv; 
    theRGBColor.blue  = theHardwareColor.blue  = cv;
#endif

#if 1  /* yellow-blue */
    theRGBColor.red   = theHardwareColor.red   = (65535 - cv)*0.95;
    theRGBColor.green = theHardwareColor.green = (65535 - cv)*0.8;
    theRGBColor.blue  = theHardwareColor.blue  = cv;
#endif
    
    theStatus = XAllocColor(theDisplay, theColormap, 
			    &theHardwareColor);
    if(theStatus) {
      theColors[i] = theHardwareColor.pixel;
      if(verbose>2) {
	printf("Pixel and RGBvalues (requested/received) = ");
	printf("%d, %d, %d, %d, %d, %d, %d \n", 
	       (int)theHardwareColor.pixel, 
	       (int)theRGBColor.red, 
	       (int)theRGBColor.green, 
	       (int)theRGBColor.blue,
	       (int)theHardwareColor.red, 
	       (int)theHardwareColor.green, 
	       (int)theHardwareColor.blue);
      }
    }
    else {
      printf("Failed ! Pixel and RGBvalues (requested) = %d, %d, %d, %d \n", 
	     (int)theHardwareColor.pixel, 
	     (int)theRGBColor.red, 
	     (int)theRGBColor.green, 
	     (int)theRGBColor.blue);
    }
  }
}
/*---------------------------------------------------------------------------*/

void quitX()

     /* 
      *  Free memory and close the theDisplay.
      */
{
  XFreeGC(theDisplay,theGC);
  if(u_plot) XFreePixmap(theDisplay,theUPixmap);
  if(u_plot) XFreePixmap(theDisplay,theVPixmap);
  if(v_plot) XFreePixmap(theDisplay,theWPixmap);
  XCloseDisplay(theDisplay);
}
/*---------------------------------------------------------------------------*/

int keyboard_chk()

     /* 
      *  Checks keyboard event queue and takes appropriate action. 
      */
{
  XEvent         theEvent;
  KeySym         theKeySym;
  int            theKeyBufferMaxLen = 64;
  char           theKeyBuffer[65];
  XComposeStatus theComposeStatus;
  int            n_up=0, n_rt=0;

  if(!(u_plot||v_plot)) return(0);   /* neither field is plotted so return 0 */

  while(XCheckWindowEvent(theDisplay, theWindow, EV_MASK, &theEvent) ) {
    /* If an event then ... */

    switch(theEvent.type) {

#if FOCUS_CK
      /*  This case is necessary for window managers which do not set 
       *  keyboard focus to theWindow when the pointer enters theWindow.
       *  See also ezgraph.h.
       */
    case EnterNotify:
      XSetInputFocus(theDisplay, theWindow, RevertToPointerRoot, CurrentTime); 
      break;
#endif
      
    case KeyPress:
      /* A KeyPress event has occurred.  Take appropriate action */

      XLookupString((XKeyEvent *)&theEvent,   theKeyBuffer, theKeyBufferMaxLen,
		    &theKeySym, &theComposeStatus);
      switch(theKeySym) {
      case XK_Up:     n_up++; break;                           /* move up    */
      case XK_Down:   n_up--; break;                           /* move down  */
      case XK_Right:  n_rt++; break;                           /* move right */
      case XK_Left:   n_rt--; break;                           /* move left  */
      case XK_Escape: quitX(); return(1);                           /* quit  */
      case XK_Q:      quitX(); return(1);                           /* quit  */
      case XK_q:      quitX(); return(1);                           /* quit  */
      case XK_P:                                                    /* pause */  
	XWindowEvent(theDisplay, theWindow, EV_MASK, &theEvent);
	break;      
      case XK_p:                                                    /* pause */  
	XWindowEvent(theDisplay, theWindow, EV_MASK, &theEvent);
	break;      
      }
      break;
    }
  }

  /* If a move has been requested, call mover() */
  if(n_up||n_rt) mover(n_up, n_rt);  
  return(0);
}
/*---------------------------------------------------------------------------*/

void mover (int n_up, int n_rt)

     /* 
      *  Moves the spiral(s) in response to a key press. 
      *  To save time, v_lap is not moved.  I have not implemented "wrap 
      *  around" moving in the case of periodic boundary conditions.
      */
{
  int n_move,i,j;
  
  if( n_up > 0 ) {                                                     /* up */
    n_move = n_up*MOVESTEP;
    for(j=ndim;j>n_move;j--) {
      for(i=1;i<=ndim;i++) {
	u[i][j] = u[i][j-n_move]; 
	v[i][j] = v[i][j-n_move];
	w[i][j] = w[i][j-n_move];
	u_lap[k2][i][j] = u_lap[k2][i][j-n_move]; 
	v_lap[k2][i][j] = v_lap[k2][i][j-n_move]; 
	w_lap[k2][i][j] = w_lap[k2][i][j-n_move]; 
      }
    }
  } 
  if( n_up < 0 ) {                                                   /* down */
    n_move = n_up*MOVESTEP;
    for(j=1;j<=ndim+n_move;j++) {
      for(i=1;i<=ndim;i++) {
	u[i][j] = u[i][j-n_move]; 
	v[i][j] = v[i][j-n_move];
	w[i][j] = w[i][j-n_move];
	u_lap[k2][i][j] = u_lap[k2][i][j-n_move]; 
	v_lap[k2][i][j] = v_lap[k2][i][j-n_move]; 
	w_lap[k2][i][j] = w_lap[k2][i][j-n_move]; 
      }
    }
  } 
  if( n_rt > 0 ) {                                                  /* right */
    n_move = n_rt*MOVESTEP;
    for(i=ndim;i>n_move;i--) {
      for(j=1;j<=ndim;j++) {
	u[i][j] = u[i-n_move][j]; 
	v[i][j] = v[i-n_move][j];
	w[i][j] = w[i-n_move][j];
	u_lap[k2][i][j] = u_lap[k2][i-n_move][j]; 
	v_lap[k2][i][j] = v_lap[k2][i-n_move][j]; 
	w_lap[k2][i][j] = w_lap[k2][i-n_move][j]; 
      }
    }
  }
  if( n_rt < 0 ) {                                                   /* left */
    n_move = n_rt*MOVESTEP;
    for(i=1;i<=ndim+n_move;i++) {
      for(j=1;j<=ndim;j++) {
	u[i][j] = u[i-n_move][j]; 
	v[i][j] = v[i-n_move][j];
	w[i][j] = w[i-n_move][j];
	u_lap[k2][i][j] = u_lap[k2][i-n_move][j]; 
	v_lap[k2][i][j] = v_lap[k2][i-n_move][j]; 
	w_lap[k2][i][j] = w_lap[k2][i-n_move][j]; 
      }
    }
  } 

  ntips = 0;
  plot();
}
/*---------------------------------------------------------------------------*/

#else   /* else #if GRAPHICS */

/* if NOT GRAPHICS, define these dummy routines. */

void plot_ini    () {}
void plot        () {}
int  keyboard_chk() {return(0);}

#endif

/*---------------------------------------------------------------------------*/

