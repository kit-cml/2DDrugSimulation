#include "apdr.hpp"
#include "rhs.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <cvode/cvode.h>
#include <cvode/cvode_dense.h>
#include <nvector/nvector_serial.h>

#if  defined TN2004ENDO || defined TN2004M || defined TN2004EPI || defined TN2006ENDO || defined TN2006M || defined TN2006EPI
int APDR( int argc, char **argv, param_t *p_param )
{
  /* CellML object pointer */
  CellML *p_cell;
  /* time related variables */
  double t_prev;
  double t_stim_start;
  double t_stim_end;
  double t_di;
  double t_depolarize;
  double t_repolarize;
  double bcl;
  double bcl_decrement; 
  double apd;
  /* volatge related variables */
  double v_prev;
  double v_top;
  double v_valley;
  double v_apd90;
  /* supporting variables */
  bool is_ap_increasing;
  bool is_action_potential;
  int depolarization_counter;
  int print_freq;
  int err_code;
  FILE *fp_apd90;
  FILE *fp_apdr;
  FILE *fp_vm;
  /* SUNDIALs variables */
  int retval;
  double t_out, t_curr;
  void *cvode_mem;
  N_Vector states_vec;
  
  /* Initialize p_cell based on the ionmodel value */
  p_cell = init_cell();
  if( p_cell == NULL ){
    return 20;
  }

  /* Initialize mutation object of p_cell if possible */
  p_cell->isMutated = 0;
  if( strcmp(p_param->variant, "ori") != 0 ){
    p_cell = init_mutation( p_param->variant, p_param->mutation_scaling, p_cell);
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

  cvode_mem = CVodeCreate( CV_BDF,CV_NEWTON );
  CVodeSetUserData( cvode_mem, p_cell );
  states_vec = N_VMake_Serial( p_cell->states_size, p_cell->STATES );
  CVodeInit( cvode_mem, rhs_fn, T0, states_vec );
  t_out = T1;
  CVodeSetMaxStep( cvode_mem, p_param->dt );
  CVDense( cvode_mem, p_cell->states_size );
  CVodeSStolerances( cvode_mem, 1.0e-7, 1.0e-7 );

  bcl = p_param->bcl_init;
  bcl_decrement = p_param->bcl_decrement;
  p_cell->CONSTANTS[stim_period] = p_param->bcl_init;
  p_cell->CONSTANTS[stim_amplitude] *= p_param->stim_amp;
  p_cell->CONSTANTS[stim_duration] = p_param->stim_dur;
  is_ap_increasing = false;
  is_action_potential = false;
  t_curr = 0.0;
  t_prev = 0.0;
  t_stim_start = 0.0;
  t_stim_end = t_stim_start + p_cell->CONSTANTS[stim_duration];
  apd = 0.0;
  depolarization_counter = 0;
  v_prev = 0.0;
  v_top = 0.0;
  v_valley = 0.0;
  v_apd90 = 0.0;
  t_di = 0.0;
  print_freq = (1. / p_param->dt) * 2;
 
  fp_apdr = fopen( "apdr.plt", "w" );
  fp_apd90 = fopen( "debug.plt", "w" );
  fp_vm = fopen( "vmcheck.plt", "w" );

  int iout = 0;

  /* Print the header of the output files */
  fprintf( fp_vm, "%s %s\n", "TIME", "Vm" );
  /* Extra flag for S1 stimulus */
  p_cell->isS1 = true;

  /* core simulation loop */
  for (int i = 0; ; i++) {
    v_prev = p_cell->STATES[0];
    /* 
      Control the flow of simulation.
      In case of APDR simulation, this part will also 
      decrease the bcl based on the bcl_decrement value.
    */
    if (bcl <= 300.0) bcl_decrement = 10;
    if (t_curr <= t_prev + bcl * p_param->num_pace1) {
      if (t_curr > t_stim_end) {
        printf( "%s = %d %s = %lf %s = %lf %s = %lf\n",
                "depolarization_counter", depolarization_counter,
                "t_stim_start", t_stim_start,
                "t_stim_end", t_stim_end,
                "bcl", bcl );
        t_stim_start += bcl;
        t_stim_end = t_stim_start + p_cell->CONSTANTS[stim_duration];
      }
    }
    else {
      bcl -= p_param->bcl_decrement;
      p_cell->CONSTANTS[stim_period] = bcl;
      t_prev = t_curr;
      if (bcl < p_param->bcl_end) break;
    }

    retval = CVode( cvode_mem, t_out, states_vec, &t_curr, CV_NORMAL  );
    if( retval == CV_SUCCESS ){
      iout++;
      t_out += p_param->dt;
    }
 
    /* find out whether the simulation can produce ACTION POTENTIAL or not */
    if (p_cell->STATES[0] > v_prev && !is_ap_increasing) {
      is_ap_increasing = true;
      v_valley = v_prev;
      is_action_potential = v_top - v_valley > p_param->v_apd90;
      if (is_action_potential) {
        //v_apd90 = v_top - (0.9 * (v_top - v_valley));
        fprintf( fp_apd90, "%lf %lf %lf %lf\n",
                 round(t_curr), v_top, v_prev, v_apd90 );
      }
    }
    else if (p_cell->STATES[0] < v_prev && is_ap_increasing) {
      is_ap_increasing = false;
      v_top = v_prev;
    }
    // repolarization, when the potential moves from HIGH to LOW
    if (v_prev > v_apd90 && p_cell->STATES[0] <= v_apd90) {
      t_repolarize = t_curr;
      apd = t_repolarize - t_depolarize;
    }
    // depolarization, when the potential moves from LOW to HIGH
    if (v_prev <= v_apd90 && p_cell->STATES[0] > v_apd90) {
      depolarization_counter++;
      t_depolarize = t_curr;
      if (depolarization_counter >= 2) {
        t_di = t_depolarize - t_repolarize;
        fprintf( fp_apdr, "%d %lf %lf %lf %lf\n",
                 depolarization_counter, bcl, apd,
                 t_di, t_curr );
     }
    }
    /* print the result */
    if (i % print_freq == 0 ){ 
      fprintf( fp_vm, "%lf %lf\n", t_curr, p_cell->STATES[V] );
    }
  }


  fclose( fp_vm );
  fclose( fp_apd90 );
  fclose( fp_apdr );

  N_VDestroy(states_vec);
  CVodeFree(&cvode_mem);
 
  /* generate APD vs DiastolicInterval and APD vs BCL result files */
  create_apdr_output(p_param->num_pace1, depolarization_counter + 1);

  delete p_cell;
  free( p_param );
  return 0;
}

#endif
