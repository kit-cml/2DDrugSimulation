#include "apdr_bench.hpp"
#include "commons.hpp"
#include "globals.hpp"
#include "../cellmodels/tentusscher_noble_noble_panfilov_2004_a.hpp"
#include "../cellmodels/tentusscher_noble_noble_panfilov_2004_b.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <map>

#include <cvode/cvode.h>
#include <cvode/cvode_dense.h>

void create_apdr_output(int num_pace1, int number_of_lines);

using std::multimap;

void apdr_bench( int argc, char **argv, param_t *p_param )
{
  bool is_replaced = false;

  // I/O variables
  int print_freq;
  char buffer[150];
  char *pch;
  FILE *fp_vm;
  FILE *fp_vm_all;
  FILE *fp_apdr;
  FILE *fp_debug;

  // Cell object variables
  patch_clamp *p_cell;

  // SUNDIALs variables
  int cvode_retval, iout, imax;
  bool cvode_firsttime;
  double tnext, tcurr;
  void *cvode_mem;
  N_Vector states_vec;

  // APDR-related variables
  bool is_alternant;
  double vm_peak, vm_valley, vm_repol, vm_valley_prev;
  double t_depol, t_depol_prev, t_repol, t_repol_prev, apd90, t_di, t_peak;
  double bcl;
  int pace,all_pace;
  multimap<double, double> vm_data;

  // map for saving the AP graph of 1 pace
  

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

  if( (int)floor(p_param->celltype) == 0 ){
    p_cell = new tentusscher_noble_noble_panfilov_2004_a();
    printf("Using Myocardial cell!\n");
  }
  else if( (int)floor(p_param->celltype) == 1 ){
    p_cell = new tentusscher_noble_noble_panfilov_2004_b();
    printf("Using Epicardial cell!\n");
  }
  p_cell->initConsts();
  cvode_mem = CVodeCreate( CV_BDF,CV_NEWTON );
  cvode_firsttime = true;
  bcl = p_param->bcl_init;

  fp_apdr = fopen( "result/apdr.plt", "w" );
  fp_vm_all = fopen( "result/vmcheck.plt", "w");
  fp_debug = fopen( "result/debug.plt","w" );

  tnext = p_param->dt;
  tcurr = 0.0;
  print_freq = (1. / p_param->dt) * p_param->dt_write;
  vm_valley = vm_peak = p_cell->STATES[V];
  vm_repol = 0;
  all_pace = 0;
  time_t begin = time(NULL);
 
 
  while( bcl >= p_param->bcl_end ){ // start BCL loop

    p_cell->CONSTANTS[stim_period] = bcl;
    printf("Current BCL: %lf --- Next BCL: %lf\n", bcl, bcl-p_param->bcl_decrement);
    snprintf(buffer, sizeof(buffer), "result/vmcheck_%.2lfBCL.plt", bcl);
    fp_vm = fopen( buffer, "w");

    // setup CVode setting
    tcurr = 0.0;
    tnext = p_param->dt;
    states_vec = N_VMake_Serial( p_cell->states_size, p_cell->STATES );
    if( cvode_firsttime == true ) {
      CVodeInit( cvode_mem, rhs_fn, tcurr, states_vec );
      CVodeSetUserData( cvode_mem, p_cell );
      CVodeSStolerances( cvode_mem, 1.0e-6, 1.0e-6 );
      CVodeSetMaxStep( cvode_mem, p_param->dt );
      CVDense( cvode_mem, p_cell->states_size );
      cvode_firsttime = false;
    }
    else CVodeReInit( cvode_mem, tcurr, states_vec );

    // repeated variable init
    pace = 0;

    iout = 0;
    imax = ( p_param->num_pace1 * p_cell->CONSTANTS[stim_period] ) / p_param->dt;
    t_depol = 0;
    t_repol = 0;
    while( iout <= imax ){ // begin pacing loop

      // solve ODE
      cvode_retval = CVode( cvode_mem, tnext, states_vec, &tcurr, CV_NORMAL  );
      // get the RATES array
      p_cell->computeRates(tnext, p_cell->CONSTANTS, p_cell->RATES, N_VGetArrayPointer_Serial(states_vec), p_cell->ALGEBRAIC);
      if( cvode_retval == CV_SUCCESS ){
        iout++;
        tnext += p_param->dt;
      }
      else{
        printf("CVode calculation error at BCL %lf in %.2lf msec!!!\n", p_cell->CONSTANTS[stim_period], tcurr);
        break;
      }

      vm_data.insert( std::pair<double, double> (tcurr, p_cell->STATES[V]) );

      // calculate vm_peak and vm_repol
      if( ((p_cell->CONSTANTS[stim_period]*pace)+p_cell->CONSTANTS[stim_start]+4.) < tcurr && tcurr < ((p_cell->CONSTANTS[stim_period]*pace)+p_cell->CONSTANTS[stim_start]+56.) ) {
          if( vm_peak < p_cell->STATES[V] ) {
            vm_peak = p_cell->STATES[V];
            t_peak = tcurr;
          }
      }
      else if( tcurr > ((p_cell->CONSTANTS[stim_period]*pace)+p_cell->CONSTANTS[stim_start]+56.) ){
        vm_repol = vm_peak - (0.9 * (vm_peak - vm_valley));
      }


      // repolarization phase
      if( tcurr > ((p_cell->CONSTANTS[stim_period]*pace)+p_cell->CONSTANTS[stim_start]+56.) && vm_repol >= p_cell->STATES[V] && p_cell->STATES[V] > vm_repol-2 ){
        t_depol = ((p_cell->CONSTANTS[stim_period]*pace)+p_cell->CONSTANTS[stim_start]);
        t_repol = tcurr;
        apd90 = t_repol - t_depol;
        is_replaced = true;
      }

      // depolarization phase
      if( iout % (int)( p_cell->CONSTANTS[stim_period] / p_param->dt ) == 0) {

        // recalibration in case vm_peak is less than 0
        if( vm_peak < 0 || (vm_peak >= 0 && vm_peak < 10) ){
          is_replaced = false;
          for(std::multimap<double, double>::iterator itrmap = vm_data.begin(); itrmap != vm_data.end() ; itrmap++ ){
            if( vm_peak < itrmap->second ) {
              vm_peak = itrmap->second;
              vm_repol = vm_peak - (0.9 * (vm_peak - vm_valley));
              t_peak = itrmap->first;
            }
          }

          for(std::multimap<double, double>::iterator itrmap = vm_data.begin(); itrmap != vm_data.end() ; itrmap++ ){
            if( itrmap->first > t_peak && vm_repol >= itrmap->second && itrmap->second > vm_repol-1 ){
              t_depol = ((p_cell->CONSTANTS[stim_period]*pace)+p_cell->CONSTANTS[stim_start]);
              t_repol = itrmap->first;
              apd90 = t_repol - t_depol;
              is_replaced = true;
            }
          }

        }
        // if still not replaced after recalibration, just use tcurr for the repol time
        if( is_replaced == false ){
          t_depol = ((p_cell->CONSTANTS[stim_period]*pace)+p_cell->CONSTANTS[stim_start]);
          t_repol = tcurr;
          apd90 = t_repol - t_depol;
        }

        vm_valley = p_cell->STATES[V];
        //fprintf( fp_debug, "%lf  %d %lf %lf %lf %lf %lf %lf %lf\n", bcl, pace, apd90, t_repol, t_depol, t_peak, vm_peak, vm_valley, vm_repol );
        t_depol_prev = t_depol;
        t_repol_prev = t_repol;
        vm_peak = -999;
        vm_repol = 0;
        pace++;
        all_pace++;
        is_replaced = false;
        vm_data.clear();
        if(pace > 2){
          t_di = (tcurr+p_cell->CONSTANTS[stim_start]) - t_repol;
          fprintf( fp_apdr, "%d %lf %lf %lf %lf\n", all_pace, bcl, apd90, t_di, tcurr );
        }
      }

      // print vm graph
      if(iout % print_freq == 0){
        fprintf(fp_vm, "%lf %lf\n", tcurr, p_cell->STATES[V]);
        fprintf(fp_vm_all, "%lf %lf\n", tcurr, p_cell->STATES[V]);
      }

    } // end pacing loop


    // clean memory
    N_VDestroy(states_vec);
    fclose(fp_vm);

    // decrease the BCL
    if (bcl <= 300.0) p_param->bcl_decrement = 10;
    bcl -= p_param->bcl_decrement;

  } // end BCL loop

  time_t end = time(NULL);
  // end simulation loop
  printf("Time elapsed for simulation is: %lf minutes\n", (double)((end - begin)/60.));

  /* generate APD vs DiastolicInterval and APD vs BCL result files */
  create_apdr_output(p_param->num_pace1, all_pace);
 
  // Memory Cleanup
  CVodeFree(&cvode_mem);
  delete p_cell;
  free(p_param);

  fclose( fp_debug );
  fclose( fp_vm_all );
  fclose( fp_apdr );

}
