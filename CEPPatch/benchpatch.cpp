#include "patches/tentusscher_noble_noble_panfilov_2004_a.hpp"
#include "patches/Ohara_Rudy_2011.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <cvode/cvode.h>
#include <cvode/cvode_dense.h>
#include <nvector/nvector_serial.h>

double glob_vm = 0.0;


int rhs_fn( realtype t, N_Vector y, N_Vector ydot, void *user_data );

int rhs_fn( realtype t, N_Vector y, N_Vector ydot, void *user_data )
{
  NV_Ith_S(y, V) = glob_vm;
  patch_clamp *data = (patch_clamp*)user_data;
  data->computeRates( t,
                      data->CONSTANTS,
                      N_VGetArrayPointer_Serial(ydot),
                      N_VGetArrayPointer_Serial(y),
                      data->ALGEBRAIC );
  return 0;
}


int main( int argc, char **argv )
{
  if(argc < 2){
    printf("Please provide the membrane potential input file!!\n");
    return 1;
  }

  // I/O variables
  char buff[255];
  char *pch;
  FILE *fp_curr;
  FILE *fp_vm;

  // Patch clamp variables
  patch_clamp *p_patch;
  double bcl,dt;

  // SUNDIALs variables
  int retval;
  double tnext, tcurr;
  void *cvode_main;
  N_Vector states_vec;

  p_patch = new Ohara_Rudy_2011();
  p_patch->initConsts();
  bcl = 1000.0;
  dt = 0.001;
  tcurr = 0.0;

  // Create CVODE solver
  cvode_main = CVodeCreate( CV_BDF,CV_NEWTON );
  // Give p_cell as CVode User Data
  CVodeSetUserData( cvode_main, p_patch );
  // Create the states vector based on the STATES array
  states_vec = N_VMake_Serial( p_patch->states_size, p_patch->STATES );
  // Initalize CVODE solver
  CVodeInit( cvode_main, rhs_fn, 0.0, states_vec );
  // Set up the future time (must be more than 0.0)
  tnext = dt;
  // Set up the maximum step size (dt)
  CVodeSetMaxStep( cvode_main, dt );
  // Set up the linear solver type
  CVDense( cvode_main, p_patch->states_size );
  // Set up the numerical error tolerances
  CVodeSStolerances( cvode_main, 1.0e-7, 1.0e-7 );

  fp_curr = fopen("curr_ORd_verapamil_CaK.plt", "w");
  fp_vm = fopen(argv[1], "r");
  fgets( buff, 255, fp_vm );

  fprintf(fp_curr, "%s %s\n", "TIME", "ICaK");

  while(fgets( buff, 255, fp_vm ) != NULL){
    pch = strtok(buff, ", ");
    pch = strtok(NULL, ", ");

    // Grab the Vm from input file
    glob_vm = strtod(pch, NULL);

    retval = CVode( cvode_main, tnext, states_vec, &tcurr, CV_NORMAL  );
    if( retval == CV_SUCCESS ){
      tnext += dt;
    }
    else{
      break;
    }

    fprintf(fp_curr, "%lf %lf\n", tcurr, p_patch->ALGEBRAIC[ICaK]);

  }

  fclose(fp_curr);
  fclose(fp_vm);

  delete p_patch;
  printf("FINISHED!\n");

  return 0;
}
