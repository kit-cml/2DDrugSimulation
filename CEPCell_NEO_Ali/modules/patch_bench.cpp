#include "patch_bench.hpp"
#include "commons.hpp"
#include "../cellmodels/ohara_rudy_cipa_v1_2017.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <cvode/cvode.h>
#include <cvode/cvode_dense.h>
#include <nvector/nvector_serial.h>


void patch_bench(int argc, char **argv, param_t *p_param)
{
  patch_clamp *p_patch;

  double bcl,dt;
  int pace;
  int iout, imax, idx;

  FILE *fp_res;

  // SUNDIALs variables
  int retval;
  double tnext, tcurr;
  void *cvode_main;
  N_Vector states_vec;

  p_patch = new ohara_rudy_cipa_v1_2017();
  p_patch->initConsts();

  pace = 1;
  bcl = 2000.0;
  dt = 0.1;

  tcurr = 0.0;
  iout = 0;
  imax = ( pace * bcl ) / dt;

  // Set up tnext and CVode options
  tnext = dt;
  cvode_main = CVodeCreate( CV_BDF,CV_NEWTON );
  CVodeSetUserData( cvode_main, p_patch );
  states_vec = N_VMake_Serial( p_patch->states_size, p_patch->STATES );
  CVodeInit( cvode_main, rhs_fn, 0.0, states_vec );
  CVodeSetMaxStep( cvode_main, 0.5 );
  CVDense( cvode_main, p_patch->states_size );
  CVodeSStolerances( cvode_main, 1.0e-7, 1.0e-7 );

  fp_res = fopen("result.plt", "w");
  fprintf(fp_res, "%s %s %s\n", "TIME", "Curr", "Vm");

  while(iout < imax){
    fprintf(fp_res, "%lf %lf %lf\n", tcurr, p_patch->ALGEBRAIC[INa], p_patch->STATES[V]);
    retval = CVode( cvode_main, tnext, states_vec, &tcurr, CV_NORMAL  );
    if( retval == CV_SUCCESS ){
      iout++;
      tnext += dt;
    }
    else{
      break;
    }
  }

  fclose(fp_res);

  delete p_patch;
  printf("FINISHED!\n");

}
