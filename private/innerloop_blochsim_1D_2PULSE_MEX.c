/* [ff_list, gg_list, mm_list, mm_init_list] = ...              */
/* innerloop_blochsim_1D_2PULSE(zz_list, grad, dt, RF_EXC, ...  */
/* numTRs, E1, CRUSH_FACTOR, RF_REF, E1p);                      */
/*                                                              */
/* called by blochsim_1D_2PULSE (in the directory above)        */
/*                                                              */
/* units:                                                       */
/* -- zz_list and TH: cm                                        */
/* -- dt: sec                                                   */
/* -- grad: gauss/cm                                            */
/*                                                              */
/* (Mukund Balasubramanian, 2025)                               */

#include "mex.h"
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

void mexFunction(int nlhs,       mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  /* ************** Inputs *******************/
  double *zz_list, *grad, *dt, *RF_EXC, *numTRs, *E1;
  double *CRUSH_FACTOR, *RF_REF, *E1p;
  
  /* ************** Outputs *******************/
  double *ff_list, *gg_list, *mm_list, *mm_init_list;
  
  /* ************** Functions **************** */
  void blochsim_loop1D(/* Inputs */
                       int num_zz,
                       double *zz_list,
                       double grad,
                       double dt,
                       int num_tt,
                       double *RF_EXC,
                       int numTRs,
                       double E1,
                       double CRUSH_FACTOR,
                       double *RF_REF,
                       double E1p,
                       /* Outputs */
                       double *ff_list,
                       double *gg_list,
                       double *mm_list,
                       double *mm_init_list);
  
  /* ************** Others ******************* */
  int num_zz, num_tt;

  /* Inputs */
  zz_list = mxGetPr(prhs[0]);
  grad = mxGetPr(prhs[1]);
  dt = mxGetPr(prhs[2]);
  RF_EXC = mxGetPr(prhs[3]);
  numTRs = mxGetPr(prhs[4]);
  E1 = mxGetPr(prhs[5]);
  CRUSH_FACTOR = mxGetPr(prhs[6]);
  RF_REF = mxGetPr(prhs[7]);
  E1p = mxGetPr(prhs[8]);

  /* Derived from Inputs */
  num_zz = mxGetM(prhs[0]);
  num_tt = mxGetM(prhs[3]);

  /* Outputs */
  plhs[0] = mxCreateDoubleMatrix(num_zz,1,mxREAL);
  ff_list = mxGetPr(plhs[0]);

  plhs[1] = mxCreateDoubleMatrix(num_zz,1,mxREAL);
  gg_list = mxGetPr(plhs[1]);

  plhs[2] = mxCreateDoubleMatrix(num_zz,1,mxREAL);
  mm_list = mxGetPr(plhs[2]);
  
  plhs[3] = mxCreateDoubleMatrix(num_zz,1,mxREAL);
  mm_init_list = mxGetPr(plhs[3]);

  /* Call Function */
  blochsim_loop1D(num_zz, zz_list, grad[0], dt[0], num_tt, RF_EXC,
                  numTRs[0], E1[0], CRUSH_FACTOR[0], RF_REF, E1p[0],
                  ff_list, gg_list, mm_list, mm_init_list);

}

/* ------------------------------------------------------------ */
void blochsim_loop1D(/* Inputs */
                     int num_zz,
                     double *zz_list,
                     double grad,
                     double dt,
                     int num_tt,
                     double *RF_EXC,
                     int numTRs,
                     double E1,
                     double CRUSH_FACTOR,
                     double *RF_REF,
                     double E1p,
                     /* Outputs */
                     double *ff_list,
                     double *gg_list,
                     double *mm_list,
                     double *mm_init_list)
/* ------------------------------------------------------------ */
{

  int ii, jj, kk;
  double zz;
  double gamma = 2 * M_PI * 4258;
  double Delta_w;
  double Mrot[3][3];
  double fgm[3];
  double tmp[2];
  double phi;

  void RF_matrix(double M[3][3], double Delta_w, double RF, double dt);
  void RF_matrix_times_fgm(double M[3][3], double fgm[3]);

  /* initialize output arrays */
  for(kk = 0; kk < num_zz; kk++)
    {
      ff_list[kk] = 0.0;
      gg_list[kk] = 0.0;
      mm_list[kk] = 0.0;
      mm_init_list[kk] = 1.0;
    }

  /* loop through TRs */
  for(ii = 0; ii < numTRs; ii++)
    {

      /* loop through z locations */
      for(kk = 0; kk < num_zz; kk++)
        {
          /* get zz and fgm for this voxel */
          zz = zz_list[kk];
          fgm[0] = 0.0;
          fgm[1] = 0.0;
          fgm[2] = mm_init_list[kk];

          /* get the Delta_w from the slice-select gradient */
          Delta_w = gamma * grad * zz;

          /* excitation RF pulse */
          for(jj = 0; jj < num_tt; jj++)
            {
              /* Mrot = RF_matrix(Delta_w, RF_EXC(jj), dt) */
              RF_matrix(Mrot, Delta_w, RF_EXC[jj], dt); 

              /* fgm = Mrot * fgm */
              RF_matrix_times_fgm(Mrot, fgm);
            }

          /* rewind: T_EXC = dt * num_tt */
          tmp[0] = fgm[0];
          tmp[1] = fgm[1];
          phi = gamma * grad * zz * dt * num_tt / 2;
          fgm[0] = tmp[0]*cos(phi) - tmp[1]*sin(phi);
          fgm[1] = tmp[0]*sin(phi) + tmp[1]*cos(phi);

          /* transverse relaxation between EXC and REF */
          /* (ignore for now) */

          /* longitudinal relaxation between EXC and REF */
          fgm[2] = (1-E1) + fgm[2] * E1;

          /* crush left before REF */
          tmp[0] = fgm[0];
          tmp[1] = fgm[1];
          phi = gamma * grad * zz * dt * num_tt * CRUSH_FACTOR;
          fgm[0] = tmp[0]*cos(phi) - tmp[1]*sin(phi);
          fgm[1] = tmp[0]*sin(phi) + tmp[1]*cos(phi);

          /* refocusing RF pulse */
          for(jj = 0; jj < num_tt; jj++)
            {
              /* Mrot = RF_matrix(Delta_w, RF_REF(jj), dt) */
              RF_matrix(Mrot, Delta_w, RF_REF[jj], dt); 

              /* fgm = Mrot * fgm */
              RF_matrix_times_fgm(Mrot, fgm);
            }

          /* crush right after REF */
          tmp[0] = fgm[0];
          tmp[1] = fgm[1];
          phi = gamma * grad * zz * dt * num_tt * CRUSH_FACTOR;
          fgm[0] = tmp[0]*cos(phi) - tmp[1]*sin(phi);
          fgm[1] = tmp[0]*sin(phi) + tmp[1]*cos(phi);

          /* transverse relaxation between REF and TE */
          /* (ignore for now) */

          /* longitudinal relaxation between REF and TE */
          fgm[2] = (1-E1) + fgm[2] * E1;

          /* store magnetization at TE for the kkth zz location */
          ff_list[kk] = fgm[0];
          gg_list[kk] = fgm[1];
          mm_list[kk] = fgm[2];
  
          /* Mz at end of TR will be the new Mz_init */
          mm_init_list[kk] = (1-E1p) + fgm[2] * E1p;

        } /* matches for kk */

    } /* matches for ii */

}

/* ------------------------------------------------------------ */
void RF_matrix(double M[3][3], double Delta_w, double RF, double dt)
/* implements Appendix B of Mulkern and Williams (1993)         */
/* (w --> Delta_w, w1 --> RF, b --> beta)                       */
/*                                                              */
/* NOTE: RF is assumed to be real-valued, corresponding to      */
/* a pulse along x in the rotating frame                        */
/*                                                              */
/* (RF_matrix.m is the MATLAB version of this function)         */
/* ------------------------------------------------------------ */
{
  double beta;
  double Delta_w_over_beta;
  double RF_over_beta;

  beta = sqrt(Delta_w*Delta_w + RF*RF);

  /* should we be using fabs and EPS below? */
  if((Delta_w == 0) && (beta == 0))
    Delta_w_over_beta = 1.0;
  else
    Delta_w_over_beta = Delta_w / beta;

  /* should we be using fabs and EPS below? */
  if((RF == 0) && (beta == 0))
    RF_over_beta = 1.0;
  else
    RF_over_beta = RF / beta;

  M[0][0] = 1 - Delta_w_over_beta*Delta_w_over_beta * (1-cos(beta*dt)); 
  M[1][1] = cos(beta*dt);
  M[2][2] = 1 - RF_over_beta*RF_over_beta * (1-cos(beta*dt)); 

  M[0][1] = Delta_w_over_beta * sin(beta*dt);
  M[1][0] = -M[0][1];

  M[0][2] = RF_over_beta * Delta_w_over_beta * (1-cos(beta*dt)); 
  M[2][0] = +M[0][2];

  M[1][2] = RF_over_beta * sin(beta*dt);
  M[2][1] = -M[1][2];

}

/* ------------------------------------------------------------ */
void RF_matrix_times_fgm(double M[3][3], double fgm[3])
/* ------------------------------------------------------------ */
{
  double tmp[3];

  tmp[0] = fgm[0];
  tmp[1] = fgm[1];
  tmp[2] = fgm[2];

  fgm[0] = M[0][0]*tmp[0] + M[0][1]*tmp[1] + M[0][2]*tmp[2];
  fgm[1] = M[1][0]*tmp[0] + M[1][1]*tmp[1] + M[1][2]*tmp[2];
  fgm[2] = M[2][0]*tmp[0] + M[2][1]*tmp[1] + M[2][2]*tmp[2];
}
