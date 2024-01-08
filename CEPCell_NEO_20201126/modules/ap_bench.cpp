#include "ap_bench.hpp"
#include "commons.hpp"
#include "../cellmodels/tentusscher_noble_noble_panfilov_2004_b.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <cvode/cvode.h>
#include <cvode/cvode_dense.h>
#include <sys/stat.h>

int file_exists (char *filename);

void ap_bench( int argc, char **argv, param_t *p_param )
{
  // I/O variables
  int print_freq;
  char buff[255];
  char *pch;
  FILE *fp_curr;
  FILE *fp_curr10;
  FILE *fp_vm;
  FILE *fp_vm10;
  FILE *fp_states;

  // Patch clamp variables
  patch_clamp *p_cell;

  // SUNDIALs variables
  int retval;
  double tnext, tcurr;
  void *cvode_main;
  N_Vector states_vec;

  int i;
  p_cell = new tentusscher_noble_noble_panfilov_2004_b();
  p_cell->initConsts();
  p_cell->CONSTANTS[stim_period] = p_param->bcl_init;

  char output_name[100] = "output.dat";
  if( file_exists(output_name) && p_param->is_using_output > 0 ){
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
  cvode_main = CVodeCreate( CV_BDF,CV_NEWTON );
  // Give p_cell as CVode User Data
  CVodeSetUserData( cvode_main, p_cell );
  // Create the states vector based on the STATES array
  states_vec = N_VMake_Serial( p_cell->states_size, p_cell->STATES );
  // Initalize CVODE solver
  CVodeInit( cvode_main, rhs_fn, 0.0, states_vec );
  // Set up the maximum step size (dt)
  CVodeSetMaxStep( cvode_main, p_param->dt );
  // Set up the linear solver type
  CVDense( cvode_main, p_cell->states_size );
  // Set up the numerical error tolerances
  CVodeSStolerances( cvode_main, 1.0e-7, 1.0e-7 );

  fp_curr = fopen("ires.plt", "w");
  fp_curr10 = fopen("ires10.plt", "w");
  fp_vm = fopen( "vmcheck.plt", "w");
  fp_vm10 = fopen( "vmcheck10.plt", "w");

  fprintf(fp_vm, "%s %s\n", "TIME", "Vm");
  fprintf(fp_vm10, "%s %s\n", "TIME", "Vm");
  fprintf(fp_curr, "%s %s %s %s %s %s %s\n",
           "TIME","INa","ICaL","Ito","IKr","IKs","IK1");
  fprintf(fp_curr10, "%s %s %s %s %s %s %s\n",
           "TIME","INa","ICaL","Ito","IKr","IKs","IK1");

  int iout = 0;
  int imax = ( p_param->num_pace1 * p_cell->CONSTANTS[stim_period] ) / p_param->dt;
  tnext = p_param->dt;
  tcurr = 0.0;
  print_freq = (1. / p_param->dt) * p_param->dt_write;
  time_t begin = time(NULL);

  while( iout <= imax ){

    retval = CVode( cvode_main, tnext, states_vec, &tcurr, CV_NORMAL  );
    if( retval == CV_SUCCESS ){
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
    if( p_param->num_pace1 > 10 ){
      if( tcurr >= ( (p_param->num_pace1-10) * p_cell->CONSTANTS[stim_period] ) && iout % print_freq == 0 ){
        fprintf(fp_curr10, "%lf %lf %lf %lf %lf %lf %lf\n",
             tcurr,p_cell->ALGEBRAIC[i_Na],p_cell->ALGEBRAIC[i_CaL],p_cell->ALGEBRAIC[i_to],p_cell->ALGEBRAIC[i_Kr],p_cell->ALGEBRAIC[i_Ks],p_cell->ALGEBRAIC[i_K1]);
        fprintf(fp_vm10, "%lf %lf\n", tcurr, p_cell->STATES[V]);
      }
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
  CVodeFree(&cvode_main);
  fclose(fp_curr);
  fclose(fp_curr10);
  fclose(fp_vm);
  fclose(fp_vm10);
  delete p_cell;
  free(p_param);

  printf("FINISHED!\n");

}

int file_exists (char *filename) {
  struct stat   buffer;   
  return (stat (filename, &buffer) == 0);
}
