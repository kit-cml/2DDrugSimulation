#include "ep.hpp"
#include "rhs.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <cvode/cvode.h>
#include <cvode/cvode_dense.h>
#include <nvector/nvector_serial.h>

int EP( int argc, char **argv, param_t *p_param )
{
  /* CellML object pointer */
  CellML *p_cell;
  /* time related variables */
  int print_freq;
  int err_code;
  FILE *fp_ca;
  FILE *fp_vm;
  FILE *fp_ires;

  /* SUNDIALs variables */
  int retval;
  double tnext, tcurr;
  void *cvode_main;
  N_Vector states_vec;

  /* Initialize p_cell based on the ionmodel value */
  p_cell = init_cell( );
  if( p_cell == NULL ){
    return EXIT_FAILURE;
  }

  /* Initialize mutation object of p_cell if possible */
  p_cell->isMutated = 0;
  if( strcmp(p_param->variant, "ori") != 0 ){
    p_cell = init_mutation( p_param->variant, 0.5, p_cell);
    if( p_cell->mutation != NULL ){
      p_cell->isMutated = 1;
    }
  }


  /*
     Invoke the states' initial condition of the cell model
     References:
     https://www.cellml.org/getting-started/tutorials/usinggeneratedc
  */
  p_cell->initConsts();
  p_cell->isS1 = 1;
  p_cell->CONSTANTS[stim_period] = p_param->bcl_init;
  p_cell->CONSTANTS[stim_amplitude] *= p_param->stim_amp;

  // Create CVODE solver
  cvode_main = CVodeCreate( CV_BDF,CV_NEWTON );
  // Give p_cell as CVode User Data
  CVodeSetUserData( cvode_main, p_cell );
  // Create the states vector based on the STATES array
  states_vec = N_VMake_Serial( p_cell->states_size, p_cell->STATES );
  // Initalize CVODE solver
  CVodeInit( cvode_main, rhs_fn, 0.0, states_vec );
  // Set up the future time (must be more than 0.0)
  tnext = p_param->dt;
  // Set up the maximum step size (dt)
  CVodeSetMaxStep( cvode_main, p_param->dt );
  // Set up the linear solver type
  CVDense( cvode_main, p_cell->states_size );
  // Set up the numerical error tolerances
  CVodeSStolerances( cvode_main, 1.0e-7, 1.0e-7 );

  int iout = 0;
  int imax = ( p_param->num_pace1 * p_cell->CONSTANTS[stim_period] ) / p_param->dt;
  print_freq = (1. / p_param->dt) * p_param->dt_write;
  fp_ca = fopen( "ca_i_result.plt", "w" );
  fp_vm = fopen( "vmcheck.plt", "w" );
  fp_ires = fopen( "ires.plt", "w" );
  fprintf( fp_ca, "%s %s\n", "TIME", "Ca_i" );
  fprintf( fp_vm, "%s %s \n",
           "TIME",
           "Vm" );
#if defined ORUDY2011_STATIC || defined ORUDY2017_DYNAMIC
  fprintf( fp_ires, "%s %s %s %s %s %s %s %s %s\n",
           "TIME",
           "INaL",
           "ICaL",
           "Ito",
           "IKr",
           "IKs",
           "IK1",
           "INaCa",
           "INaK" );
#endif
  // begin simulation loop
  time_t begin = time(NULL);
  while(1){

    retval = CVode( cvode_main, tnext, states_vec, &tcurr, CV_NORMAL  );
    if( retval == CV_SUCCESS ){
      iout++;
      tnext += p_param->dt;
    }

    if(iout % print_freq == 0){
      fprintf( fp_vm, "%lf %lf\n", tcurr, p_cell->STATES[V]);
      fprintf( fp_ca, "%lf %lf\n",
               tcurr,
               p_cell->STATES[Ca_i] );
#if defined ORUDY2011_STATIC || defined ORUDY2017_DYNAMIC
      fprintf( fp_ires, "%lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
               tcurr,
               p_cell->ALGEBRAIC[i_NaL],
               p_cell->ALGEBRAIC[i_CaL],
               p_cell->ALGEBRAIC[i_to],
               p_cell->ALGEBRAIC[i_Kr],
               p_cell->ALGEBRAIC[i_Ks],
               p_cell->ALGEBRAIC[i_K1],
               p_cell->ALGEBRAIC[i_NaCa_i],
               p_cell->ALGEBRAIC[i_NaK] );
#endif

    }

   if (iout >= imax) break;
  }

  time_t end = time(NULL);
  // end simulation loop
  printf("Time elapsed for simulation is: %ld sec", end - begin);

  fclose( fp_vm );
  fclose( fp_ca );
  fclose( fp_ires );

  // Memory Cleanup
  N_VDestroy(states_vec);
  CVodeFree(&cvode_main);
  delete p_cell;

  return 0;
}
