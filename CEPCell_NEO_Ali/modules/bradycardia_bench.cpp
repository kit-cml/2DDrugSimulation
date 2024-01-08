#include "bradycardia_bench.hpp"
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
#include <sys/stat.h>

void bradycardia_bench( int argc, char **argv, param_t *p_param )
{
  // I/O variables
  int print_freq;
  char buff[255];
  char *pch;
  FILE *fp_curr;
  FILE *fp_vm;
  FILE *fp_states;
  FILE *fp_res;

  // Patch clamp variables
  patch_clamp *p_cell;

  // SUNDIALs variables
  int retval;
  double tnext, tcurr;
  void *cvode_mem;
  N_Vector states_vec;

  // supporting variables
  double dvmdt_max, dvmdt_repol;
  double vm_peak, vm_valley, vm_repol90, vm_repol30;
  double t_depol, t_repol90, t_repol30, apd90;
  int pace;

  int i;
  p_cell = new tentusscher_noble_noble_panfilov_2004_b();
  p_cell->initConsts();
  p_cell->CONSTANTS[stim_period] = p_param->bcl_init;

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
  CVodeSStolerances( cvode_mem, 1.0e-7, 1.0e-7 );

  if( mympi::rank == 0 ){
    fp_curr = fopen("ires.plt", "w");
    fp_vm = fopen( "vmcheck.plt", "w");
    fp_res = fopen( "result.plt", "w");
  }

  fprintf(fp_vm, "%s %s\n", "TIME", "Vm");
  fprintf(fp_curr, "%s %s %s %s %s %s %s\n",
           "TIME","INa","ICaL","Ito","IKr","IKs","IK1");
  fprintf(fp_res, "%s %s %s %s\n", "PACE", "APD90", "DVMDT_MAX", "DVMDT_REPOL");

  int iout = 0;
  int imax = ( p_param->num_pace1 * p_cell->CONSTANTS[stim_period] ) / p_param->dt;
  tnext = p_param->dt;
  tcurr = 0.0;
  print_freq = (1. / p_param->dt) * p_param->dt_write;
  vm_valley = p_cell->STATES[V];
  vm_peak = vm_repol90 = dvmdt_max = 0 ;
  dvmdt_repol = -999.;

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
    default:
      printf("Using A1656D WILDTYPE mutation!\n");
      break;
  }
  
  time_t begin = time(NULL);
  while( iout <= imax ){

    retval = CVode( cvode_mem, tnext, states_vec, &tcurr, CV_NORMAL  );
    p_cell->computeRates(tnext, p_cell->CONSTANTS, p_cell->RATES, N_VGetArrayPointer_Serial(states_vec), p_cell->ALGEBRAIC);
    if( retval == CV_SUCCESS ){
      iout++;
      tnext += p_param->dt;
    }
    else{
      printf("CVode calculation error at %lf msec!!!\n", tcurr);
      break;
    }

    if( ((p_cell->CONSTANTS[stim_period]*pace)+p_cell->CONSTANTS[stim_start]+4.) < tcurr
        && tcurr < ((p_cell->CONSTANTS[stim_period]*pace)+p_cell->CONSTANTS[stim_start]+56.) ) {
        if( vm_peak < p_cell->STATES[V] ) vm_peak = p_cell->STATES[V];
    }
    else if( tcurr > ((p_cell->CONSTANTS[stim_period]*pace)+p_cell->CONSTANTS[stim_start]+56.) ){
      vm_repol90 = vm_peak - (0.9 * (vm_peak - vm_valley));
      vm_repol30 = vm_peak - (0.3 * (vm_peak - vm_valley));
    }

    if( vm_repol90 > p_cell->STATES[V] && p_cell->STATES[V] > vm_repol90-2 ){
      t_depol = ((p_cell->CONSTANTS[stim_period]*pace)+p_cell->CONSTANTS[stim_start]);
      t_repol90 = tcurr;
      apd90 = t_repol90 - t_depol;
    }

    // get maximum dvmdt
    if( dvmdt_max < p_cell->RATES[V] ) dvmdt_max = p_cell->RATES[V];

    // get the steepest dvmdt_repol around vm_repol30 and vm_repol90 if AP shape is eligible
    if( vm_repol30 >= p_cell->STATES[V] && p_cell->STATES[V] >= vm_repol90 && dvmdt_repol < p_cell->RATES[V] ){
      //printf( "time checking repol: %lf PeakVM: %lf Higher bound: %lf Lower bound: %lf detected value: %lf before value: %lf\n", tcurr, vm_peak, vm_repol30, vm_repol90, p_cell->RATES[V], dvmdt_repol );
      dvmdt_repol = p_cell->RATES[V];
    }

    // print result
    if(iout % print_freq == 0){
      fprintf(fp_curr, "%lf %lf %lf %lf %lf %lf %lf\n",
               tcurr,p_cell->ALGEBRAIC[i_Na],p_cell->ALGEBRAIC[i_CaL],p_cell->ALGEBRAIC[i_to],p_cell->ALGEBRAIC[i_Kr],p_cell->ALGEBRAIC[i_Ks],p_cell->ALGEBRAIC[i_K1]);
      fprintf(fp_vm, "%lf %lf %lf\n", tcurr, p_cell->STATES[V], p_cell->RATES[V]);
    }

    // executed when entering new pace
    if( iout % (int)( p_cell->CONSTANTS[stim_period] / p_param->dt ) == 0 ) {
      //printf("LOWER BOUND: %lf --- HIGHER BOUND: %lf\n", 
      //       ((p_cell->CONSTANTS[stim_period]*pace)+p_cell->CONSTANTS[stim_start]+4.),
      //       ((p_cell->CONSTANTS[stim_period]*pace)+p_cell->CONSTANTS[stim_start]+56.) );
      //printf("VM PEAK PACE %d: %lf -- at 90%%: %lf\n", pace, vm_peak, vm_repol90);
      //printf("REPOL: %lf --- DEPOL: %lf --- APD90: %lf\n", t_repol90, t_depol, apd90);
      fprintf( fp_res, "%d %.2lf %lf %lf\n", pace, apd90, dvmdt_max, dvmdt_repol );
      vm_peak = vm_repol90 = vm_repol30 = 0.;
      vm_valley = p_cell->STATES[V];
      pace++;
    }

  }

  time_t end = time(NULL);
  // end simulation loop
  printf("Time elapsed for simulation is: %ld mins\n", (end - begin)/60);


  // Memory Cleanup
  N_VDestroy(states_vec);
  CVodeFree(&cvode_mem);
  delete p_cell;
  free(p_param);

  if( mympi::rank == 0 ){
    fclose(fp_curr);
    fclose(fp_vm);
    fclose(fp_res);
  }
  printf("FINISHED!\n");
}
