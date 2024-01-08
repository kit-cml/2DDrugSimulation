#include "anm_bench_custom.hpp"
#include "commons.hpp"
#include "globals.hpp"
#include "../cellmodels/tentusscher_noble_noble_panfilov_2004_b.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <map>
#include <vector>

#include <cvode/cvode.h>
#include <cvode/cvode_dense.h>

using std::map;
using std::multimap;
using std::vector;

anm_result_t anm_bench_custom(int argc, char **argv, param_t *p_param, patch_clamp* p_cell, const char* scenario_name)
{
  bool is_replaced = false;

  // Struct to store result
  anm_result_t result;

  // I/O variables
  int print_freq;
  char buffer[150];
  char *pch;
  FILE *fp_vm;
  FILE *fp_apd;
  FILE *fp_apdr;
  FILE *fp_debug;

  // SUNDIALs variables
  int cvode_retval, iout, imax;
  bool cvode_firsttime;
  double tnext, tcurr;
  void *cvode_mem;
  N_Vector states_vec;

  // APDR-related variables
  bool is_local;
  double vm_peak, vm_valley, vm_repol, vm_valley_prev;
  double t_depol, t_depol_prev, t_repol, t_repol_prev, apd90, t_di, t_peak;
  double bcl, bcl_dec, aocl, anm_aocl, anm_alternant;
  int pace,all_pace;
  multimap<double, double> vm_data;
  multimap<double, double> vm_last20_data;
  map<double, double> anm_data;
  map<double, double> mean_apd_data;
  map<double, multimap<double, double> > bcl_vm_data;
  vector<double> apd90_vec;
  map<double, vector<double> > bcl_apd_data;

  
  if( p_cell == NULL ){
    p_cell = new tentusscher_noble_noble_panfilov_2004_b();
    p_cell->initConsts();
    p_cell->CONSTANTS[g_Ks] *= p_param->gks_scale;
    p_cell->CONSTANTS[g_Kr] *= p_param->gkr_scale;
    p_cell->CONSTANTS[g_K1] *= p_param->gk1_scale;
    p_cell->CONSTANTS[g_Na] *= p_param->gna_scale;
    p_cell->CONSTANTS[g_bna] *= p_param->gbna_scale;
    p_cell->CONSTANTS[g_CaL] *= p_param->gcal_scale;
    p_cell->CONSTANTS[g_bca] *= p_param->gbca_scale;
    p_cell->CONSTANTS[g_to] *= p_param->gto_scale;
    p_cell->CONSTANTS[g_pCa] *= p_param->gpca_scale;
    p_cell->CONSTANTS[g_pK] *= p_param->gpk_scale;
    is_local = true;
  }
  else{
    printf("Cell has been created!! Skipped initialization....\n");
    printf("Scenario: %s\n", scenario_name);
    is_local = false;
  }

  cvode_mem = CVodeCreate( CV_BDF,CV_NEWTON );
  cvode_firsttime = true;
  bcl = p_param->bcl_init;
  bcl_dec = p_param->bcl_decrement; 
  
  if(scenario_name == NULL || strlen(scenario_name) == 0) {
    fp_apdr = fopen( "result/apdr.plt", "w" );
    fp_debug = fopen( "result/debug.plt","w" );
  }

  if( scenario_name == NULL ) fprintf( fp_debug, "%s %s %s %s %s %s %s %s %s\n", "bcl", "pace", "apd90", "t_repol", "t_depol", "t_peak", "vm_peak", "vm_valley", "vm_repol" );

  tnext = p_param->dt;
  tcurr = 0.0;
  print_freq = (1. / p_param->dt) * p_param->dt_write;
  vm_valley = vm_valley_prev = vm_peak = p_cell->STATES[V];
  vm_repol = 0;
  all_pace = 0;
 
  while( bcl >= p_param->bcl_end ){ // start BCL loop

    p_cell->CONSTANTS[stim_period] = bcl;
    //printf("Current BCL: %lf --- Next BCL: %lf\n", bcl, bcl-bcl_dec);

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
      /*if(pace >= p_param->num_pace1-20)*/ vm_last20_data.insert( std::pair<double, double> (tcurr, p_cell->STATES[V]) );

      // calculate vm_peak and vm_repol
      if( ((p_cell->CONSTANTS[stim_period]*pace)+p_cell->CONSTANTS[stim_start]+20.) < tcurr && tcurr < ((p_cell->CONSTANTS[stim_period]*pace)+p_cell->CONSTANTS[stim_start]+50.) ) {
          if( vm_peak < p_cell->STATES[V] ) {
            vm_peak = p_cell->STATES[V];
            t_peak = tcurr;
          }
      }
      else if( tcurr > ((p_cell->CONSTANTS[stim_period]*pace)+p_cell->CONSTANTS[stim_start]+50.) ){
        vm_repol = vm_peak - (0.9 * (vm_peak - vm_valley));
      }


      // repolarization phase
      if( tcurr > ((p_cell->CONSTANTS[stim_period]*pace)+p_cell->CONSTANTS[stim_start]+40.) && vm_repol >= p_cell->STATES[V] && p_cell->STATES[V] > vm_repol-1 ){
        t_depol = ((p_cell->CONSTANTS[stim_period]*pace)+p_cell->CONSTANTS[stim_start]);
        t_repol = tcurr;
        apd90 = t_repol - t_depol;
        is_replaced = true;
      }

      // depolarization phase

      // executed at the beginning of the pace
      if( iout > 0 && iout % (int)( p_cell->CONSTANTS[stim_period] / p_param->dt ) == 0) {

        // recalibration in case vm_peak is less than 0 or around 0 to 9.99... mV
        if( vm_peak < 0 || (vm_peak >= 0 && vm_peak < 10) ){
          is_replaced = false;
          for(std::multimap<double, double>::iterator itrmap = vm_data.begin(); itrmap != vm_data.end() ; itrmap++ ){
            if( vm_peak < itrmap->second ) {
              vm_peak = itrmap->second;
              vm_repol = vm_peak - (0.9 * (vm_peak - vm_valley));
              t_peak = itrmap->first;
            }
          }

          // itrmap->first means time, itrmap->second means Vm(time)
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
        if(is_replaced == false)
        {
          t_depol = ((p_cell->CONSTANTS[stim_period]*pace)+p_cell->CONSTANTS[stim_start]);
          t_repol = tcurr;
          apd90 = t_repol - t_depol;
        }
        //if((int)bcl == 250) printf("RECALIBRATING VM_PEAK NEGATIVE: BCL %lf PACE %d VM_PEAK %lf VM_VALLEY %lf VM_REPOL %lf T_REPOL %lf T_DEPOL %lf T_CURR %lf APD90 %lf \n", bcl, pace, vm_peak, vm_valley, vm_repol, t_repol, t_depol, tcurr, apd90);

        vm_valley = p_cell->STATES[V];
        if( scenario_name == NULL ) fprintf( fp_debug, "%lf  %d %lf %lf %lf %lf %lf %lf %lf\n", bcl, pace, apd90, t_repol, t_depol, t_peak, vm_peak, vm_valley, vm_repol );
        t_depol_prev = t_depol;
        t_repol_prev = t_repol;
        vm_peak = -999;
        vm_repol = 0;
        pace++;
        all_pace++;
        vm_data.clear();
        is_replaced = false;
        apd90_vec.push_back( apd90 );
        if(pace > 2){
          t_di = (tcurr+p_cell->CONSTANTS[stim_start]) - t_repol;
          if( scenario_name == NULL ) fprintf( fp_apdr, "%d %lf %lf %lf %lf\n", all_pace, bcl, apd90, t_di, tcurr );
        }
      }

    } // end pacing loop

    // clean memory
    N_VDestroy(states_vec);

    // calculate ANM and mean APD for the corresponding BCL
    anm_data.insert( std::pair<double, double> (bcl, get_ANM(apd90_vec)) );
    mean_apd_data.insert( std::pair<double, double> (bcl, average(apd90_vec, apd90_vec.size()-10-1, apd90_vec.size()-1)) );

    // insert the last 20 paces for printing
    bcl_vm_data.insert( std::pair<double, multimap<double, double> > (bcl, vm_last20_data) );
    bcl_apd_data.insert( std::pair<double, vector<double> > (bcl, apd90_vec) );

    // decrease the BCL
    if (bcl <= 350.0) bcl_dec = 10.;
    bcl -= bcl_dec;
    apd90_vec.clear();
    vm_last20_data.clear();

  } // end BCL loop


  aocl = anm_aocl = anm_alternant = 0.0;
  // itrmap->first represent BCL value
  // itrmap->second represent ANM value
  if( scenario_name == NULL ) fprintf(fp_debug, "%s %s\n", "bcl", "ANM");
  for(std::map<double, double>::iterator itrmap = anm_data.begin(), itrmap_next = ++anm_data.begin(); itrmap_next != anm_data.end() ; itrmap++, itrmap_next++ ){
    if( itrmap->second > 0.05 && aocl < itrmap_next->first ) {
      printf("Current biggest one: %lf ANM: %lf Next one: %lf ANM: %lf\n", itrmap->first, itrmap->second, itrmap_next->first, itrmap_next->second);
      aocl = itrmap_next->first;
      anm_aocl = itrmap_next->second;
      anm_alternant = itrmap->second;
    }
    if( scenario_name == NULL ) fprintf(fp_debug, "%lf %lf\n", itrmap_next->first, itrmap_next->second);
  }

  result.aocl = aocl;
  result.mean_apd = mean_apd_data[aocl];
  result.anm_aocl = anm_aocl;
  result.anm_alternant = anm_alternant;

/*
  if( scenario_name == NULL || strlen(scenario_name) == 0 ) snprintf(buffer, sizeof(buffer), "result/vmcheck_%.2lfBCL.plt", aocl);
  else snprintf(buffer, sizeof(buffer), "result/vmcheck_%.2lfBCL_%s.plt", aocl, scenario_name);
  if( p_param->is_print_graph == 1 ){
    fp_vm = fopen( buffer, "w");
    for(std::multimap<double, double>::iterator itrmap = bcl_vm_data[aocl].begin(); itrmap != bcl_vm_data[aocl].end() ; itrmap++ ){
      fprintf(fp_vm, "%lf %lf\n", itrmap->first, itrmap->second);
    }
    fclose( fp_vm );
  }
*/

  double bcl_print_iter = p_param->bcl_print_start;
  while ( bcl_print_iter >= p_param->bcl_print_end ){
    snprintf(buffer, sizeof(buffer), "result/vmcheck_%.2lfBCL.plt", bcl_print_iter);
    fp_vm = fopen( buffer, "w");
    snprintf(buffer, sizeof(buffer), "result/apd_list_%.2lfBCL.plt", bcl_print_iter);
    fp_apd = fopen( buffer, "w");
    for(std::multimap<double, double>::iterator itrmap = bcl_vm_data[bcl_print_iter].begin(); itrmap != bcl_vm_data[bcl_print_iter].end() ; itrmap++ ){
      fprintf(fp_vm, "%lf %lf\n", itrmap->first, itrmap->second);
    }
    for(std::vector<double>::iterator itrvec = bcl_apd_data[bcl_print_iter].begin(); itrvec != bcl_apd_data[bcl_print_iter].end() ; itrvec++ ){
      fprintf(fp_apd, "%lf\n", *itrvec);
    }
    fprintf(fp_apd, "ANM: %lf meanAPD: %lf\n", anm_data[bcl_print_iter], mean_apd_data[bcl_print_iter]);
    fclose( fp_apd );
    fclose( fp_vm );
    if(bcl_print_iter > 350) bcl_print_iter -= p_param->bcl_decrement;
    else bcl_print_iter -= 10;
  }

  create_apdr_output(p_param->num_pace1, all_pace);

  // Memory Cleanup
  CVodeFree(&cvode_mem);
  if(is_local == true) delete p_cell;
  if( scenario_name == NULL )fclose( fp_debug );
  if( scenario_name == NULL )fclose( fp_apdr );

  return result;
}
