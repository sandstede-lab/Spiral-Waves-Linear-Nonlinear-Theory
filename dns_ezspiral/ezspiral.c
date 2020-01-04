/* ------------------------------------------------------------------------- *
 *                                                                           *
 *                                EZ-SPIRAL                                  *
 *                    A Code for Simulating Spiral Waves                     *
 *                                                                           *
 *       Copyright (C) 1992, 1994, 1997, 1998, 2007                          *
 *                                            Dwight Barkley                 *
 *                                            D.Barkley@warwick.ac.uk        *
 *                                                                           *
 *                                Version 3                                  *
 *                                                                           *
 * ------------------------------------------------------------------------- *
 * RCS Information
 * ---------------------------
 * $Revision: 3.2.1.1 $
 * $Date: 2007/05/07 09:50:42 $
 * ------------------------------------------------------------------------- */

/* changed Bjorn Sandstede: suffix .dat changed to .txt */

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include "ezspiral.h"

/* 
 * Global variables used throughout the EZ-Spiral Code (see ezspiral.h)
 * -------------------------------------------------------------------- */

Real *fields, *u, *v, *sigma_u, *sigma_v, *tip_x, *tip_y;
Real  a_param, b_param, one_o_eps, delta, grid_h, dt, 
      dt_o_eps, dt_o_2eps, one_o_a, b_o_a, dt_o_wh2, dtDv_o_wh2;
int   nx, ny, nz, field_size, nsteps, istep, write_tip, 
      graphic_resolution, verbose, ntips; 

/* Global variables for this file only
 * ----------------------------------- */

static  Real  length_x, length_y, Dv, ts;
static  int   plot_step, hist_step, hist_x, hist_y, binary_write;
static  FILE *history_file, *tip_file; 

/* Private functions 
 * ----------------- */

static void  Initialize      (void);
static void  Allocate_memory (void);
static void  Finish_up       (void);
static void  Write_history   (int wrt_step);
static void  Write_slice     (void);   /* added Bjorn Sandstede */
static void  Write_fc        (void);
static void  Read_ic         (void);
static void  Generate_ic     (void);

#define NEXT_LINE(fp) while(getc(fp)!='\n');  /* macro used to skip to end 
						 of input line */

/* ========================================================================= */

int main(void)
{
  Initialize();

  for(istep=0; istep<nsteps; istep++) {
    if( Event_check() )          break;
    if( (istep%plot_step) == 0 ) Draw();
/*    if( hist_step && (istep%hist_step)==0 ) Write_history(istep); */
	if( hist_step && (istep%hist_step)==0 ) Write_slice(); /* added Bjorn Sandstede */
    Step();
  }

  Finish_up();
  return 0;
}
/* ========================================================================= */

void Initialize (void)
{
  /* Read task file, open output files, and initialize graphics and time
   * stepping.  */

  double p_in;
  FILE *fp;
  int initial_field;

  /* ----------------------------- 
   * Read parameters from task.txt 
   * ----------------------------- */

  if((fp=fopen("task.txt","r"))==NULL) {
    if(verbose) fprintf(stderr,"Cannot open task file: task.txt \n");
    exit(1);
  }
  else {
    fscanf(fp,"%lg",&p_in);  NEXT_LINE(fp); a_param=p_in;
    fscanf(fp,"%lg",&p_in);  NEXT_LINE(fp); b_param=p_in;
    fscanf(fp,"%lg",&p_in);  NEXT_LINE(fp); one_o_eps=p_in;
    fscanf(fp,"%lg",&p_in);  NEXT_LINE(fp); length_x=p_in;
    fscanf(fp,"%lg",&p_in);  NEXT_LINE(fp); Dv=p_in;
                             NEXT_LINE(fp);
    fscanf(fp,"%d", &nx);    
    fscanf(fp,",%d",&ny);    NEXT_LINE(fp);
    fscanf(fp,"%lg",&p_in);  NEXT_LINE(fp); ts=p_in;
    fscanf(fp,"%lg",&p_in);  NEXT_LINE(fp); delta=p_in;
                             NEXT_LINE(fp);
    fscanf(fp,"%d",&nsteps); NEXT_LINE(fp);
    fscanf(fp,"%d",&plot_step);       NEXT_LINE(fp);
    fscanf(fp,"%d",&write_tip);  NEXT_LINE(fp);
    fscanf(fp,"%d",&hist_step);       NEXT_LINE(fp);
    fscanf(fp,"%d,%d",&hist_x,&hist_y);     NEXT_LINE(fp);
                                            NEXT_LINE(fp);
    fscanf(fp,"%d",&initial_field);         NEXT_LINE(fp);
    fscanf(fp,"%d",&graphic_resolution);    NEXT_LINE(fp); 
    fscanf(fp,"%d",&binary_write);          NEXT_LINE(fp);
                                            NEXT_LINE(fp);
    fscanf(fp,"%d",&verbose);

    fclose(fp);
  }

  /* ----------------------------------------------------------------- 
   * Process input parameters and write various things to screen 
   *
   * Note: I choose to read length_x and then use NY and NZ to
   * determine length_y and length_z.  You can easily change this, the
   * only requirement is that grid_h be the same for each direction.
   * ----------------------------------------------------------------- */

#if NINEPOINT 
  #define WEIGHT 6.
  #define STABILITY_LIMIT (3.*grid_h*grid_h/8.) 
#else
  #define WEIGHT 1.
  #define STABILITY_LIMIT (grid_h*grid_h/4.)
#endif

  grid_h    = length_x/(NX-1);
  length_y  = grid_h * (NY-1);

  dt         = ts*STABILITY_LIMIT;
  dt_o_wh2   = dt/(WEIGHT*grid_h*grid_h);
  dtDv_o_wh2 = Dv*dt_o_wh2;
  one_o_a    = 1.0/a_param; 
  b_o_a      = b_param/a_param; 
  dt_o_eps   = dt*one_o_eps; 
  dt_o_2eps  = dt_o_eps/2.; 

  if(verbose) {
    printf("\n\nModel Parameters: \n");
    printf("a     = %g\n", a_param);
    printf("b     = %g\n", b_param);
    printf("eps   = 1/%g = %g\n", one_o_eps, 1.0/one_o_eps);
    printf("L_x, L_y = %g, %g\n", length_x, length_y);
    printf("Dv    = %g\n", Dv);
    printf("\nNumerical Parameters: \n");
    printf("NX = %-6d NY = %-6d ts  = %-10g delta  = %-10g\n", 
	   NX, NY, ts, delta);
    printf("dt   = %-10g dt/eps = %-10g grid_h = %-10g\n", 
	   dt, dt_o_eps, grid_h);

#if NINEPOINT
    printf("\n-- Using 9 point Laplacians");
#else
    printf("\n-- Using 5 point Laplacians");
#endif

#if EXPLICIT
    printf("; explicit kinetics");
#else
    printf("; implicit kinetics");
#endif
#if SPLIT
    printf("; operator splitting --\n");
#else
    printf(" --\n");
#endif

    printf("\nNumber of time steps = %d\n", nsteps);
    printf("  time steps per plot (and/or tip computation) = %d\n", 
	   plot_step);
    if(write_tip) printf("  writing tip data\n");
    else printf("  not writing tip data\n");
    if(hist_step) printf("  writing history data\n");
    else printf("  not writing history data\n");
  }

  /* ------------------ 
   * Perform some tests 
   * ------------------ */

  if( V_DIFF_ON && Dv==0.) {
    fprintf(stderr,"***** V_DIFF_ON is 1 and Dv == 0. ******\n");
    exit(1);
  }
  if(!V_DIFF_ON && Dv!=0.) {
    fprintf(stderr,"***** V_DIFF_ON is 0 and Dv != 0. ******\n");
    exit(1);
  }
  if(hist_step && (
     hist_x<1 || hist_x>NX || 
     hist_y<1 || hist_y>NY ) ) {
    fprintf(stderr,"***** history point out of range ******\n");
    exit(1);
  }
  if(ts > 1.0 ) {
    fprintf(stderr,"***** ts > 1 (the diffusion stability limit) ******\n");
    exit(1);
  }
  if (plot_step == 0) {
    fprintf(stderr,"***** plot_step cannot be 0 ******\n");
    plot_step = nsteps+1;
  }


  /* ------------ 
   * Final things
   * ------------ */

  Allocate_memory();
  Read_ic(); 
  Step_ini();
  if(hist_step) history_file = fopen("history.txt", "w");
  if(write_tip) tip_file = fopen("tip.txt", "w"); 
  Draw_ini(initial_field);
  Tip_ini();
}
/* ========================================================================= */

void Allocate_memory (void)
{
  /* -----------------------------------------------------------------------
   * There are NX and NY actual grid points in each direction.
   * field_size = (NX+2)*(NY+2) because of the way boundaries are
   * treated.  See ezspiral.h for definition of NX, NY.
   *                                                                         
   * I allocate one big block of memory for everything (u, v, sigma_u,
   * and if needed sigma_v) and then set appropriate pointers.  The
   * handling of memory can be changed with appropriate changes to the
   * macros INDEX etc in ezspiral.h.
   * ----------------------------------------------------------------------- */

  field_size = (NY+2)*(NX+2);

#if V_DIFF_ON
  fields = (Real *) calloc(6*field_size, sizeof(Real));
#else
  fields = (Real *) calloc(4*field_size, sizeof(Real));
#endif

  if (fields != NULL) {
    /* --------------------------------------------------- 
     * We have one big block of memory.  Set some pointers 
     * --------------------------------------------------- */
    u       = fields;
    v       = fields +   field_size;
    sigma_u = fields + 2*field_size;
#if V_DIFF_ON
    sigma_v = fields + 4*field_size;
#endif
    if(verbose>2) printf("Memory allocated in one block.\n");
  }
  else {
    /* ------------------- 
     * Did not get memory. 
     * ------------------- */
    fprintf(stderr, "\n ****Error: Not enough memory available\n");
    exit(-1);
  }
}
/* ========================================================================= */

void Finish_up (void)
{
  QuitX();
  Write_fc(); 
  if(hist_step) fclose(history_file);
  if(write_tip) fclose(tip_file); 
}
/* ========================================================================= */

void Write_tip_data (float x, float y)
{
  /* This routine for outputting the tip data is called from Find_tips() in
   * eztip.c where the tip is computed.  */

  fprintf(tip_file, "%.5f %.5f %.5f\n", dt*(Real)istep, x, y);
}
/* ========================================================================= */

void Write_slice (void) /* added Bjorn Sandstede */
{ 
  /* This routine outputs a specific slice of the spiral profile */
   * eztip.c where the tip is computed.  */

  int i,j;
  i = 516;
  for(j=1;j<=NY;j++) {
	fprintf(history_file,"%g %g\n",(float)U(i,j),(float)V(i,j));
  }
}

/* ========================================================================= */

void Write_history (int wrt_step)
{
  fprintf(history_file, "%.5f %.5f %.5f\n", dt*(Real)wrt_step, 
	  U(hist_x,hist_y), V(hist_x,hist_y) ); 

  if(verbose>2) {
    printf("istep = %d", wrt_step);
    printf("; (u,v) = %.5f, %.5f", 
	   U(hist_x,hist_y), V(hist_x,hist_y) ); 
    printf("\n");
  }
}
/* ========================================================================= */

void Write_fc (void)
{
  /* Write final condition file */

  FILE *fp;
  time_t tt1;
  int i,j;

  time(&tt1); 
  fp = fopen("fc.txt", "w");

  fprintf(fp,"Model Parameters: a, b, 1/eps, Lx, Ly, Dv = ");
  fprintf(fp,"%g, %g, %g, %g, %g, %g\n", a_param, b_param, one_o_eps, 
	  length_x, length_y, Dv);
  fprintf(fp,"Numerical Parameters: NX, NY, ts, delta = ");
  fprintf(fp,"%d, %d, %g, %g \n", NX, NY, ts, delta);
  fprintf(fp,"File written: %s",ctime(&tt1));
  fprintf(fp,"Comments: 2D \n");
  fprintf(fp,"\n");
  fprintf(fp,"\n");
  fprintf(fp,"\n");
  fprintf(fp,"\n");
  fprintf(fp,"\n");

  if(binary_write) {
    /* Write binary data */
    fprintf(fp,"Binary values of u and v follow\n");
    for(j=1;j<=NY;j++) {
      fwrite(&U(1,j), sizeof(Real), NX, fp);
    }
    for(j=1;j<=NY;j++) {
      fwrite(&V(1,j), sizeof(Real), NX, fp);
    }
  }
  else {
    /* Write ascii data */
    fprintf(fp,"Ascii values of u and v follow\n");
    for(j=1;j<=NY;j++) {
      for(i=1;i<=NX;i++) {
	fprintf(fp, "%g %g\n", (float)U(i,j), (float)V(i,j) ); 
      }
    }
  }  
  fclose(fp); 
}
/* ========================================================================= */

void Read_ic (void)
{
  /* Reads initial condition file (ic.txt) if it exists, otherwise calls
   * Generate_ic() to generate new initial condition. */

  double u_in, v_in;
  Real   *u_tmp, *v_tmp, rx, ry;
  int    nx_ic, ny_ic, npts_ic2, index, i_tmp, j_tmp, i, j;
  char   f_type, dummy;
  FILE  *fp;

  if((fp=fopen("ic.txt","r"))==NULL) { 
    Generate_ic();
    return;
  }

  /* Read nx_ic etc following = sign on second line of file */
  NEXT_LINE(fp);                                     
  while( (dummy=getc(fp)) != '=');                   
  fscanf(fp, "%d, %d", &nx_ic, &ny_ic);  

  /* Skip to 10th line and read first character to determine type 
     B(inary) or A(scii) */
  for(i=0;i<8;i++) NEXT_LINE(fp); 
  f_type = getc(fp); NEXT_LINE(fp); 
  
  if ( (f_type !='B') && (f_type !='A') ) {
    if(verbose) 
      printf("\n ic.txt exists but of unrecognized type Binary or Ascii \n"); 
    exit(1);
  }

  if(verbose) printf("\nReading ic.txt with nx, ny = %d, %d... \n\n", 
		     nx_ic, ny_ic);

  npts_ic2 = nx_ic * ny_ic;  

  /* Allocate temporary memory and read from file */

  u_tmp =(Real *) malloc((unsigned)(npts_ic2)*sizeof(Real));
  v_tmp =(Real *) malloc((unsigned)(npts_ic2)*sizeof(Real));

  if(f_type =='B') {    
    /* Binary data file */
    fread(u_tmp, sizeof(Real), npts_ic2, fp);
    fread(v_tmp, sizeof(Real), npts_ic2, fp);
  }
  else {          
    /* Ascii data file */
    for(index=0;index<npts_ic2;index++) {
      fscanf (fp, "%lg %lg\n", &u_in, &v_in);
      u_tmp[index] = u_in; v_tmp[index] = v_in;
    }
  }

  /* Copy into u and v */
  rx = (nx_ic-1.0)/(NX-1.0);
  ry = (ny_ic-1.0)/(NY-1.0);

  for(j=1;j<=NY;j++) {
    j_tmp = nx_ic * (int)(ry*(j-1));
    for(i=1;i<=NX;i++) {
      i_tmp = (int)(rx*(i-1));
      U(i,j) = u_tmp[i_tmp + j_tmp];
      V(i,j) = v_tmp[i_tmp + j_tmp];
    }
  }

  free(u_tmp);
  free(v_tmp);  
}
/* ========================================================================= */

void Generate_ic (void)
{
  /* Generate new initial condition */

  int i, j;

  if(verbose) printf("\nGenerating initial condition \n\n");

  for(i=1;i<=NX;i++) {
    for(j=1;j<=NY;j++) {
      U(i,j) = V(i,j) = 0.; 
      if (i<NX/2) V(i,j) = a_param/2.;
      if(j>(10+NY/2)) U(i,j) = 1.0;  
    }
  }
}
/* ========================================================================= */


