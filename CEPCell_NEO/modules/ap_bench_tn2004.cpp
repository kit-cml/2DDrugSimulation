#include "ap_bench_tn2004.hpp"
#include "commons.hpp"
#include "globals.hpp"
#include "../cellmodels/tentusscher_noble_noble_panfilov_2004_b.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <cvode/cvode.h>
#include <cvode/cvode_dense.h>


void ap_bench_tn2004( int argc, char **argv, param_t *p_param )
{
  // I/O variables
  int print_freq;
  char buff[255];
  char *pch;
  FILE *fp_curr;
  FILE *fp_vm;
  FILE *fp_states;

  // Patch clamp variables
  patch_clamp *p_cell;

  // SUNDIALs variables
  int cvode_retval;
  double tnext, tcurr;
  void *cvode_mem;
  N_Vector states_vec;

  int i;

  switch (glob_var::A1656D_mode)
  {
    case 1:
      printf("Using A1656D mutation!\n");
      break;
    case 2:
      printf("Using A1656D mutation with MEXILETINE!\n");
      break;
    case 3:
      printf("Using A1656D mutation with FLECAINIDE!\n");
      break;
    case 4:
      printf("Using A1656D mutation with RANOLAZINE!\n");
      break;
    case 5:
      printf("Using A1656D WILDTYPE mutation!\n");
      break;
    default:
      printf("Not using any mutation!\n");
      break;
  }
  p_cell = new tentusscher_noble_noble_panfilov_2004_b();
  p_cell->initConsts();
  p_cell->CONSTANTS[stim_period] = p_param->bcl_init;

  char output_name[100] = "output.dat";
  if( is_file_exists(output_name) && p_param->is_using_output > 0 ){
    i = 0;
    printf("using output!!!\n");
    fp_states = fopen("output.dat", "r");
    while(fgets( buff, sizeof(buff), fp_states) != NULL){
      printf("STATES[%d] before: %lf\n", i, p_cell->STATES[i]);
      p_cell->STATES[i] = strtod(buff, NULL);
      printf("STATES[%d] after: %lf\n", i, p_cell->STATES[i]);
      i++;
    }
    fclose(fp_states);
  }

  // Create CVODE solver
  cvode_mem = CVodeCreate( CV_BDF,CV_NEWTON );
  // Give p_cell as CVode User Data
  CVodeSetUserData( cvode_mem, p_cell );
  // Create the states vector based on the STATES array
  states_vec = N_VMake_Serial( p_cell->states_size, p_cell->STATES );
  // Initalize CVODE solver
  CVodeInit( cvode_mem, rhs_fn, 0.0, states_vec );
  // Set up the maximum step size (dt)
  CVodeSetMaxStep( cvode_mem, p_param->dt );
  // Set up the linear solver type
  CVDense( cvode_mem, p_cell->states_size );
  // Set up the numerical error tolerances
  CVodeSStolerances( cvode_mem, 1.0e-6, 1.0e-6 );

  fp_curr = fopen("ires.plt", "w");
  fp_vm = fopen( "vmcheck.plt", "w");

  fprintf(fp_vm, "%s %s\n", "TIME", "Vm");
  fprintf(fp_curr, "%s %s %s %s %s %s %s\n",
           "TIME","INa","ICaL","Ito","IKr","IKs","IK1");

  int iout = 0;
  int imax = ( p_param->num_pace1 * p_cell->CONSTANTS[stim_period] ) / p_param->dt;
  tnext = p_param->dt;
  tcurr = 0.0;
  print_freq = (1. / p_param->dt) * p_param->dt_write;
  time_t begin = time(NULL);


  while( iout <= imax ){

    cvode_retval = CVode( cvode_mem, tnext, states_vec, &tcurr, CV_NORMAL  );
    if( cvode_retval == CV_SUCCESS ){
      iout++;
      tnext += p_param->dt;
    }
    else{
      printf("CVode calculation error at %lf msec!!!\n", tcurr);
      break;
    }

    if(iout % print_freq == 0){
      fprintf(fp_curr, "%lf %lf %lf %lf %lf %lf %lf\n",
               tcurr,p_cell->ALGEBRAIC[i_Na],p_cell->ALGEBRAIC[i_CaL],p_cell->ALGEBRAIC[i_to],p_cell->ALGEBRAIC[i_Kr],p_cell->ALGEBRAIC[i_Ks],p_cell->ALGEBRAIC[i_K1]);
      fprintf(fp_vm, "%lf %lf\n", tcurr, p_cell->STATES[V]);
    }

  }
  time_t end = time(NULL);
  // end simulation loop
  printf("Time elapsed for simulation is: %ld sec\n", end - begin);

  // save last states
  fp_states = fopen("output.dat", "w");
  for( i = 0; i < p_cell->states_size; i++ ){
    fprintf(fp_states, "%lf\n", p_cell->STATES[i]);
  }
  fclose(fp_states);

  // Memory Cleanup
  N_VDestroy(states_vec);
  CVodeFree(&cvode_mem);
  fclose(fp_curr);
  fclose(fp_vm);
  free(p_param);
  delete p_cell;

  printf("FINISHED!\n");

}
