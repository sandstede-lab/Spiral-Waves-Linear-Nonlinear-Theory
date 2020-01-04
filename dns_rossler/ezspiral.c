/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*                                                                           */
/*                                EZ-SPIRAL                                  */
/*                    A Code for Simulating Spiral Waves                     */
/*                                                                           */
/*           Copyright (C) 1992, 1994, 1997 Dwight Barkley                   */
/*                                          barkley@maths.warwick.ac.uk      */
/*                                                                           */
/*                                Version 2                                  */
/*                                                                           */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* 
 * $Revision: 2.4 $
 * $Date: 1997/11/13 23:22:00 $
 */
/*---------------------------------------------------------------------------*/

#include "ezspiral.h"

/* Global variables used throughout the EZ-Spiral Code */
precision  **u, **v, **w,                        /* u and v fields */
           **u_lap[2], **v_lap[2], **w_lap[2],   /* u and v Laplacian fields */
           *tip_x, *tip_y,              /* spiral tip arrays */
           a_param, b_param, c_param,   /* model parameters */
           length, Dv, diffconst,        /* model parameters */
           grid_h, dt, ts,              /* numerical parameters */
           NDu, NDv;                    /* useful parameter combinations */
int        ndim,                        /* # of grid points per direction */
           u_plot, v_plot,              /* flags indicating fields to plot */
           verbose,                     /* verbosity level */
           k1, k2, ntips;               /* useful integers */

/* Global variables and functions for this file only */
static int      nits, itip, iplot, iwrt, iwl, hist_x, hist_y;
static FILE     *history, *path, *wl;
static void     initialize();
static void     finish_up ();
static void     write_dat (int iter);
static void     wavelength (int iter);
static void     write_fc  (char *filename);
static void     read_ic   (char *filename);
static void     next_line (FILE *filep);

/*---------------------------------------------------------------------------*/

int main()

     /*  Main routine. Takes nits time steps of the model  */
{
  int iter;

  initialize();

  for(iter=0;iter<nits||nits<0;iter++) {
    if(itip && (iter%itip)==0)     tip_finder();
    if((iter%iplot)==0)            plot();
    if(iwrt && (iter%iwrt)==0)     write_dat(iter);
    if(keyboard_chk())             break;
    step();
  }

  finish_up();
  return(0);
}
/*---------------------------------------------------------------------------*/

void initialize()

     /*  Opens files, reads data, initializes graphics and time stepping  */
{
  double p_in;
  FILE *fp;
  int i;

  /* Read parameters from task.txt */
  if((fp=fopen("task.txt","r"))==NULL) {
    if(verbose) fprintf(stderr,"Cannot open task file: task.txt \n");
    exit(1);
  }
  else {
    fscanf(fp,"%lg",&p_in); next_line(fp); a_param=p_in;
    fscanf(fp,"%lg",&p_in); next_line(fp); b_param=p_in;
    fscanf(fp,"%lg",&p_in); next_line(fp); c_param=p_in;
    fscanf(fp,"%lg",&p_in); next_line(fp); length=p_in;
    fscanf(fp,"%lg",&p_in); next_line(fp); Dv=p_in;
                            next_line(fp);
    fscanf(fp,"%d", &ndim); next_line(fp); 
    fscanf(fp,"%lg",&p_in); next_line(fp); ts=p_in;
                            next_line(fp);
    fscanf(fp,"%d",&nits);  next_line(fp);
    fscanf(fp,"%d",&iplot); next_line(fp);
    fscanf(fp,"%d",&itip);  next_line(fp);
    fscanf(fp,"%d",&iwrt);  next_line(fp);
    fscanf(fp,"%d",&iwl);   next_line(fp);
    fscanf(fp,"%d,%d",&hist_x,&hist_y); next_line(fp);
                            next_line(fp);
    fscanf(fp,"%d",&u_plot);next_line(fp);
    fscanf(fp,"%d",&v_plot);next_line(fp);
    fscanf(fp,"%d",&verbose);next_line(fp);
    fclose(fp);
  }

#if NBC  
  /* Neumann  boundary conditions are used */
  grid_h   = length/(ndim-1);
#else    
  /* Periodic boundary conditions are used */
  grid_h   = length/ndim;
#endif

#define STABILITY_LIMIT (0.25*grid_h*grid_h)

  diffconst = 0.4;
  dt        = ts*STABILITY_LIMIT;
  NDu       = diffconst*dt/(grid_h*grid_h);
  NDv       = Dv*NDu;

  if(verbose) {
    printf("\nModel Parameters: \n");
    printf("  a   = %g\n", a_param);
    printf("  b   = %g\n", b_param);
    printf("  c   = %g\n", c_param);
    printf("  L   = %g\n", length);
    printf("  Dv  = %g\n", Dv);
    printf("  NDu = %g\n", NDu);
    printf("\nNumerical Parameters: ");
#if NINEPOINT
    printf("(using 9 point Laplacians)");
#endif
    printf("\n");
    printf("  ndim = %-10d ts     = %-10g\n", ndim, ts);
    printf("  dt   = %-10g grid_h = %-10g\n", dt, grid_h);
    printf("  dt/(grid_h*grid_h) = %-10g\n", dt/(grid_h*grid_h) );
    printf("\nTime Steps: \n");
    printf("  total                   = %d\n", nits);
    printf("  per plotting xwindow    = %d\n", abs(iplot));
    printf("  per locating tip        = %d\n", itip);
    printf("  per writing wavelength  = %d\n", iwl);
    printf("  per writing tip/history = %d\n", iwrt);
    if(iplot<0) {  printf("  per saving xwindow      = %d\n", abs(iplot)); }
    else        {  printf("  xwindow not saved\n"); }
    if(hist_x) printf("  history point          = (%d, %d)\n",hist_x,hist_y );
    printf("\n");
  }

#if NINEPOINT 
  /* There is a 6 in the denominator of the 9-point Laplacian formula */
  NDu      = NDu/6.;
  NDv      = NDv/6.;
#endif

  /* Perform some tests */
  if( V_DIFF_ON && Dv==0.) {
    fprintf(stderr,"***** V_DIFF_ON is 1 and Dv == 0. ******\n");
    exit(1);
  }
  if(!V_DIFF_ON && Dv!=0.) {
    fprintf(stderr,"***** V_DIFF_ON is 0 and Dv != 0. ******\n");
    exit(1);
  }
  if(hist_x<0 || hist_x>ndim || hist_y<0 || hist_y>ndim) {
    fprintf(stderr,"***** history point out of range ******\n");
    exit(1);
  }
  if(ts > 1.0 ) {
    fprintf(stderr,"***** ts > 1 (the diffusion stability limit) ******\n");
    exit(1);
  }

  /*  Allocate memory for u[ndim+2][ndim+2], v[ndim+2][ndim+2], 
   *  u_lap[2][ndim+2][ndim+2], and if V_DIFF_ON, for v_lap[2][ndim+2][ndim+2]. 
   *  u and v could be of size [ndim][ndim]; however, because u_lap and v_lap 
   *  MUST be of size [2][ndim+2][ndim+2], I think it best that u and v also
   *  be of size [ndim+2][ndim+2].  The memory for each field is allocated 
   *  with the (ndim+2)*(ndim+2) locations contiguous.
   */
  u = (precision **) malloc((unsigned)(ndim+2)*sizeof(precision*));
  v = (precision **) malloc((unsigned)(ndim+2)*sizeof(precision*));
  w = (precision **) malloc((unsigned)(ndim+2)*sizeof(precision*));
  u[0] = (precision *) calloc((ndim+2)*(ndim+2), sizeof(precision));
  v[0] = (precision *) calloc((ndim+2)*(ndim+2), sizeof(precision));
  w[0] = (precision *) calloc((ndim+2)*(ndim+2), sizeof(precision));
  for(i=1;i<=ndim+1;i++){
    u[i] = u[0]+i*(ndim+2);
    v[i] = v[0]+i*(ndim+2);
    w[i] = w[0]+i*(ndim+2);
  }

  u_lap[0] = (precision **) malloc((unsigned)(ndim+2)*sizeof(precision*));
  u_lap[1] = (precision **) malloc((unsigned)(ndim+2)*sizeof(precision*));
  u_lap[0][0] = (precision *) calloc((ndim+2)*(ndim+2), sizeof(precision));
  u_lap[1][0] = (precision *) calloc((ndim+2)*(ndim+2), sizeof(precision));
  for(i=1;i<=ndim+1;i++){
    u_lap[0][i] = u_lap[0][0]+i*(ndim+2);
    u_lap[1][i] = u_lap[1][0]+i*(ndim+2);
  }

  v_lap[0] = (precision **) malloc((unsigned)(ndim+2)*sizeof(precision*));
  v_lap[1] = (precision **) malloc((unsigned)(ndim+2)*sizeof(precision*));
  v_lap[0][0] = (precision *) calloc((ndim+2)*(ndim+2), sizeof(precision));
  v_lap[1][0] = (precision *) calloc((ndim+2)*(ndim+2), sizeof(precision));
  for(i=1;i<=ndim+1;i++){
    v_lap[0][i] = v_lap[0][0]+i*(ndim+2);
    v_lap[1][i] = v_lap[1][0]+i*(ndim+2);
  }

  w_lap[0] = (precision **) malloc((unsigned)(ndim+2)*sizeof(precision*));
  w_lap[1] = (precision **) malloc((unsigned)(ndim+2)*sizeof(precision*));
  w_lap[0][0] = (precision *) calloc((ndim+2)*(ndim+2), sizeof(precision));
  w_lap[1][0] = (precision *) calloc((ndim+2)*(ndim+2), sizeof(precision));
  for(i=1;i<=ndim+1;i++){
    w_lap[0][i] = w_lap[0][0]+i*(ndim+2);
    w_lap[1][i] = w_lap[1][0]+i*(ndim+2);
  }

  ntips = 0;
  if(itip) {  /* memory for tip path. */
    tip_x=(precision *)malloc((unsigned)(nits+1)*sizeof(precision));
    tip_y=(precision *)malloc((unsigned)(nits+1)*sizeof(precision));
  }

  if(iwrt) {  /* open history, tip-path, wavelength files. */
    if(hist_x!=0) history = fopen("history.txt", "w");
    if(itip  !=0) path    = fopen("path.txt",    "w");
  }

  /* Read initial conditions and initialize time stepping and graphics */

  read_ic("ic.txt"); 
  step_ini();
  plot_ini();
}
/*---------------------------------------------------------------------------*/

void finish_up()

     /*  Writes fc, last history, tip-path data, closes data files  */
{
  write_fc("fc.txt");
  if(iwrt) {
    /*  if((nits%iwrt)==0) write_dat(nits);  */
    if(hist_x!=0) fclose(history);
    if(itip  !=0) fclose(path);
    if(iwl   !=0) fclose(wl);
  }
}
/*---------------------------------------------------------------------------*/

void write_dat(int iter)

     /* Writes history, tip-path data             */
{ int i,j;

  if(verbose>1) printf("iter = %d", iter);

/*  if(hist_x!=0) {
    i = (ndim-1)/2;  
    for(j=1;j<=ndim;j++) {
      	fprintf(history,"%g %g %g\n",(float)u[i][j],(float)v[i][j],(float)w[i][j]);
    }
  }
*/

  if(hist_x!=0) {
  	i = 300;
  	for(j=1;j<=ndim;j++) {
    	fprintf(history,"%g\n",(float)w[i][j]);
  	}
  }

  if(ntips>0) {
    fprintf(path, "%.8f %.8f %.8f \n", dt*(precision)iter, 
	    grid_h*(tip_x[ntips-1]-1), grid_h*(tip_y[ntips-1]-1) );
    if(verbose>1) printf("; tip (x,y) = %.5f %.5f", 
	    grid_h*(tip_x[ntips-1]-1), grid_h*(tip_y[ntips-1]-1) );
  }
  if(verbose>1) printf("\n");
}

/*---------------------------------------------------------------------------*/

void wavelength(int num)

  /* Computes and saves amplitudes and wavelengths to file */
{
  int i,js;

  js = (int)(ndim-tip_y[ntips-1]);

  fprintf(wl,"# time = %g \n",num*dt);
  for(i=1;i<ndim;i++) { fprintf(wl,"%d %g \n",i,w[i][js]); }
  fprintf(wl," \n \n");
}

/*---------------------------------------------------------------------------*/

void write_fc (char *filename)

     /* Writes final conditions to a file  */
{
  FILE *fp;
  time_t tt1;
  int i,j;

  time(&tt1); 

  if(verbose) printf("Writing fc.txt: ... \n");
  fp = fopen(filename, "w");
  fprintf(fp,"Model Parameters: a, b, c, L, Dv = %g, %g, %g, %g, %g\n",
	  a_param, b_param, c_param, length, Dv);
  fprintf(fp,"Numerical Parameters: ndim, ts = %d, %g \n", ndim, ts);
  fprintf(fp,"File written: %s",ctime(&tt1));
  fprintf(fp,"Comments:\n");
  fprintf(fp,"Values of u and v follow\n");
  for(i=1;i<=ndim;i++) {
    for(j=1;j<=ndim;j++) {
      fprintf (fp,"%g %g %g\n",(float)u[i][j],(float)v[i][j],(float)w[i][j]);
    }
  }
  fclose(fp);
  if(verbose) printf("                done\n\n");
}
/*---------------------------------------------------------------------------*/

void read_ic (char *filename)
     
     /*  Reads ic.txt / Generates ic if ic.txt does not exist  */
{
  precision **u_tmp, **v_tmp, **w_tmp, ratio;
  double u_in, v_in, w_in;
  int ndim_ic, i_tmp, j_tmp, dummy, i, j;
  FILE *fp;

  if((fp=fopen(filename,"r"))!=NULL) { /* if file can be opened then read it */
    next_line(fp);
    while( (dummy=getc(fp)) != '='); fscanf (fp, "%d", &ndim_ic); 
    next_line(fp); 
    next_line(fp);
    next_line(fp);
    next_line(fp);
    if(verbose) printf("Reading ic.txt: ... \n");

    /* Allocate temporary memory (appropriate for in ndim_ic) */
    u_tmp =((precision **) malloc((unsigned)(ndim_ic)*sizeof(precision*)))-1;
    v_tmp =((precision **) malloc((unsigned)(ndim_ic)*sizeof(precision*)))-1;
    w_tmp =((precision **) malloc((unsigned)(ndim_ic)*sizeof(precision*)))-1;
    for(i=1;i<=ndim_ic;i++){
      u_tmp[i]=((precision *) malloc((unsigned)(ndim_ic)*sizeof(precision)))-1;
      v_tmp[i]=((precision *) malloc((unsigned)(ndim_ic)*sizeof(precision)))-1;
      w_tmp[i]=((precision *) malloc((unsigned)(ndim_ic)*sizeof(precision)))-1;
    }

    /* Read into temporary memory */
    for(i=1;i<=ndim_ic;i++) {
      for(j=1;j<=ndim_ic;j++) {
       	fscanf (fp, "%lg %lg %lg\n", &u_in, &v_in, &w_in);
       	u_tmp[i][j] = u_in; v_tmp[i][j] = v_in; w_tmp[i][j] = w_in;
      }
    }
    
    /* Copy into u and v */
    ratio = (ndim_ic-1.0)/(ndim-1.0);
    for(i=1;i<=ndim;i++) {
      i_tmp = 1+(int)(ratio*(i-1));
      for(j=1;j<=ndim;j++) {
	j_tmp = 1+(int)(ratio*(j-1));
	u[i][j] = u_tmp[i_tmp][j_tmp];
	v[i][j] = v_tmp[i_tmp][j_tmp];
	w[i][j] = w_tmp[i_tmp][j_tmp];
      }
    }
    
    /* Free temporary memory */
    for(i=ndim_ic;i>=1;i--) free((char*)(u_tmp[i]+1));
    for(i=ndim_ic;i>=1;i--) free((char*)(v_tmp[i]+1));
    for(i=ndim_ic;i>=1;i--) free((char*)(w_tmp[i]+1));
    free((char*) (u_tmp+1));
    free((char*) (v_tmp+1));
    free((char*) (w_tmp+1));

    if(verbose) printf("                done\n\n");
  }
  else {     
    /* If file cannot be opened then generate new initial condition */
    if(verbose) printf("No file ic.txt: generating initial condition ...\n");
    /*    for(i=1;i<=ndim;i++) {
      for(j=1;j<=ndim;j++) {
        u[i][j] = 0.0; v[i][j] = 0.0; w[i][j] = 0.0;
      }
      } */
    for(i=1;i<=ndim;i++) {
      for(j=1;j<=ndim;j++) {
        u[i][j] = -4.0*(ndim-j)/ndim+4.0*j/ndim;
      }
    }
    for(i=1;i<=ndim;i++) {
      for(j=1;j<=ndim;j++) {
        v[i][j] = -4.0*(ndim-i)/ndim+4.0*i/ndim;
      }
    }
    w_in = (c_param-sqrt(c_param*c_param-4.0*a_param*b_param))/2.0/a_param;
    for(i=1;i<=ndim;i++) {
      for(j=1;j<=ndim;j++) {
        w[i][j] = w_in;
      }
    }
    if(verbose) printf("                done\n\n");
  }

}
/*---------------------------------------------------------------------------*/

void next_line(FILE *filep)
     
     /*  Skips to next input line  */
{
  int dummy;
  while( (dummy=getc(filep)) != '\n');
}
/*---------------------------------------------------------------------------*/
