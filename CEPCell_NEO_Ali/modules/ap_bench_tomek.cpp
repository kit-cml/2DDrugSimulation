#include "ap_bench_tomek.hpp"
#include "commons.hpp"
#include "../cellmodels/Tomek_model.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <cvode/cvode.h>
#include <cvode/cvode_dense.h>


void ap_bench_tomek( int argc, char **argv, param_t *p_param )
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

  int idx;

  p_cell = new Tomek_model();
  p_cell->initConsts(p_param->celltype);

  p_cell->CONSTANTS[BCL] = p_param->bcl_init;
  p_cell->CONSTANTS[duration] = p_param->stim_dur;

  char output_name[100] = "output.dat";
  if( p_param->is_using_output > 0 ){
    idx = 0;
    printf("using output!!!\n");
    fp_states = fopen("output.dat", "r");
    if(fp_states == NULL){
      fprintf(stderr, "Initial steady state file error!! Exiting...\n");
      return;
    }
    while(fgets( buff, sizeof(buff), fp_states) != NULL){
      printf("STATES[%d] before: %lf\n", idx, p_cell->STATES[idx]);
      p_cell->STATES[idx] = strtod(buff, NULL);
      printf("STATES[%d] after: %lf\n", idx, p_cell->STATES[idx]);
      idx++;
    }
    fclose(fp_states);
  }

  tnext = p_param->dt;
  tcurr = 0.0;
  print_freq = (1. / p_param->dt) * p_param->dt_write;

  // Create CVODE solver
  cvode_mem = CVodeCreate( CV_BDF,CV_NEWTON );
  // Give p_cell as CVode User Data
  CVodeSetUserData( cvode_mem, p_cell );
  // Create the states vector based on the STATES array
  states_vec = N_VMake_Serial( p_cell->states_size, p_cell->STATES );
  // Initalize CVODE solver
  CVodeInit( cvode_mem, rhs_fn, tcurr, states_vec );
  // Set up the maximum step size (dt)
  CVodeSetMaxStep( cvode_mem,  p_param->dt );
  // Set up the linear solver type
  CVDense( cvode_mem, p_cell->states_size );
  // Set up the numerical error tolerances
  CVodeSStolerances( cvode_mem, 1.0e-6, 1.0e-6 );

  fp_curr = fopen("ires.plt", "w");
  fp_vm = fopen( "vmcheck.plt", "w");

  fprintf(fp_vm, "%s %s\n", "TIME", "Vm");
  fprintf(fp_curr, "%s %s %s %s %s %s %s %s\n","TIME", "IStim","INa","ICaL","Ito","IKr","IKs","IK1");

  time_t begin = time(NULL);
  int iout = 0;
  int imax = (p_cell->CONSTANTS[BCL] * p_param->num_pace1) / p_param->dt;
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

    if(iout % print_freq ==0){
      fprintf(fp_curr, "%lf %lf %lf %lf %lf %lf %lf %lf\n",
             tcurr, p_cell->ALGEBRAIC[Istim],p_cell->ALGEBRAIC[INaL],p_cell->ALGEBRAIC[ICaL],p_cell->ALGEBRAIC[Ito],p_cell->ALGEBRAIC[IKr],p_cell->ALGEBRAIC[IKs],p_cell->ALGEBRAIC[IK1]);
      fprintf(fp_vm, "%lf %lf\n", tcurr, p_cell->STATES[V]);
    }
  }
  time_t end = time(NULL);

  // end simulation loop
  printf("Time elapsed for simulation is: %ld sec\n", end - begin);

  // save last states
  fp_states = fopen("output.dat", "w");
  for( idx = 0; idx < p_cell->states_size; idx++ ){
    fprintf(fp_states, "%.15lf\n", p_cell->STATES[idx]);
  }
  fclose(fp_states);

  // Memory Cleanup
  N_VDestroy(states_vec);
  CVodeFree(&cvode_mem);
  fclose(fp_curr);
  fclose(fp_vm);
  delete p_cell;
  free(p_param);

  printf("FINISHED!\n");

}
