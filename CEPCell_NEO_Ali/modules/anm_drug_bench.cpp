#include "anm_drug_bench.hpp"
#include "commons.hpp"
#include "globals.hpp"
#include "../cellmodels/Ohara_Rudy_2011.hpp"

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

anm_result_t anm_drug_bench(int argc, char **argv, param_t *p_param, patch_clamp* p_cell, const int sample_id, const double concentration)
{
  bool is_replaced = false;

  // Struct to store result
  anm_result_t result;

  // I/O variables
  int print_freq;
  char buffer[150];
  char *pch;
  FILE *fp_vm;
  FILE *fp_apdr;
  FILE *fp_debug;

  // SUNDIALs variables
  int cvode_retval, iout, imax;
  bool cvode_firsttime;
  double tnext, tcurr;
  void *cvode_mem;
  N_Vector states_vec;

  // APDR-related variables
  double vm_peak, vm_valley, vm_repol, vm_valley_prev;
  double t_depol, t_depol_prev, t_repol, t_repol_prev, apd90, t_di, t_peak;
  double bcl, bcl_dec, aocl, anm_aocl, anm_alternant, cl_alternant;
  int pace,all_pace;
  multimap<double, double> vm_data;
  multimap<double, double> vm_last20_data;
  map<double, double> anm_data;
  map<double, double> mean_apd_data;
  map<double, multimap<double, double> > bcl_vm_data;
  vector<double> apd90_vec;


  cvode_mem = CVodeCreate( CV_BDF,CV_NEWTON );
  cvode_firsttime = true;
  bcl = p_param->bcl_init;
  bcl_dec = p_param->bcl_decrement;

  tnext = p_param->dt;
  tcurr = 0.0;
  print_freq = (1. / p_param->dt) * p_param->dt_write;
  vm_valley = vm_peak = p_cell->STATES[V];
  vm_repol = 0;
  all_pace = 0;

  while( bcl >= p_param->bcl_end ){ // start BCL loop

    p_cell->CONSTANTS[BCL] = bcl;
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
    imax = ( p_param->num_pace1 * p_cell->CONSTANTS[BCL] ) / p_param->dt;
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
        printf("CVode calculation error at BCL %lf in %.2lf msec!!!\n", p_cell->CONSTANTS[BCL], tcurr);
        break;
      }

      vm_data.insert( std::pair<double, double> (tcurr, p_cell->STATES[V]) );
      if(pace >= p_param->num_pace1-20) vm_last20_data.insert( std::pair<double, double> (tcurr, p_cell->STATES[V]) );

      // calculate vm_peak and vm_repol
      if( ((p_cell->CONSTANTS[BCL]*pace)+p_cell->CONSTANTS[stim_start]+20.) < tcurr && tcurr < ((p_cell->CONSTANTS[BCL]*pace)+p_cell->CONSTANTS[stim_start]+50.) ) {
          if( vm_peak < p_cell->STATES[V] ) {
            vm_peak = p_cell->STATES[V];
            t_peak = tcurr;
          }
      }
      else if( tcurr > ((p_cell->CONSTANTS[BCL]*pace)+p_cell->CONSTANTS[stim_start]+50.) ){
        vm_repol = vm_peak - (0.9 * (vm_peak - vm_valley));
      }

      // repolarization phase
      if( tcurr > ((p_cell->CONSTANTS[BCL]*pace)+p_cell->CONSTANTS[stim_start]+40.) && vm_repol >= p_cell->STATES[V] && p_cell->STATES[V] > vm_repol-1 ){
        t_depol = ((p_cell->CONSTANTS[BCL]*pace)+p_cell->CONSTANTS[stim_start]);
        t_repol = tcurr;
        apd90 = t_repol - t_depol;
        is_replaced = true;
      }

      // depolarization phase
      if( iout % (int)( p_cell->CONSTANTS[BCL] / p_param->dt ) == 0) {

        //if((int)floor(bcl) == 250) printf("BCL %lf PACE %d VM_PEAK %lf VM_VALLEY %lf VM_REPOL %lf T_REPOL %lf T_DEPOL %lf: APD: %lf\n", bcl, pace, vm_peak, vm_valley, vm_repol, t_repol, t_depol, apd90);
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
              t_depol = ((p_cell->CONSTANTS[BCL]*pace)+p_cell->CONSTANTS[stim_start]);
              t_repol = itrmap->first;
              apd90 = t_repol - t_depol;
              is_replaced = true;
            }
          }

        }
        // if still not replaced after recalibration, just use tcurr for the repol time
        if(is_replaced == false)
        {
          t_depol = ((p_cell->CONSTANTS[BCL]*pace)+p_cell->CONSTANTS[stim_start]);
          t_repol = tcurr;
          apd90 = t_repol - t_depol;
        }
        //if((int)bcl == 250) printf("RECALIBRATING VM_PEAK NEGATIVE: BCL %lf PACE %d VM_PEAK %lf VM_VALLEY %lf VM_REPOL %lf T_REPOL %lf T_DEPOL %lf T_CURR %lf APD90 %lf \n", bcl, pace, vm_peak, vm_valley, vm_repol, t_repol, t_depol, tcurr, apd90);

        vm_valley = p_cell->STATES[V];
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

    // decrease the BCL
    if (bcl <= 350.0) bcl_dec = 10.;
    bcl -= bcl_dec;
    apd90_vec.clear();
    vm_last20_data.clear();

  } // end BCL loop

  aocl = anm_aocl = anm_alternant = cl_alternant = result.mean_apd_alternant =  0.0;
  // itrmap->first represent BCL value
  // itrmap->second represent ANM value
  for(std::map<double, double>::iterator itrmap = anm_data.begin(), itrmap_next = ++anm_data.begin(); itrmap_next != anm_data.end() ; itrmap++, itrmap_next++ ){
    if( itrmap->second > 0.05 && aocl < itrmap_next->first ) {
      aocl = itrmap_next->first;
      anm_aocl = itrmap_next->second;
      cl_alternant = itrmap->first;
      anm_alternant = itrmap->second;
    }
  }

  result.aocl = aocl;
  result.mean_apd = mean_apd_data[aocl];
  result.anm_aocl = anm_aocl;
  result.anm_alternant = anm_alternant;
  result.cl_alternant = cl_alternant;
  result.mean_apd_alternant = mean_apd_data[cl_alternant];



  // Memory Cleanup
  CVodeFree(&cvode_mem);
  return result;
}
