#include "qnet_bench.hpp"
#include "commons.hpp"
#include "../cellmodels/Ohara_Rudy_2011.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <cvode/cvode.h>
#include <cvode/cvode_dense.h>

cipa_t qnet_all_bench(int argc, char **argv, param_t *p_param, patch_clamp* p_cell, const char* scenario_name)
{
  // Struct to store result
  cipa_t cipa_result, temp_result;

  // I/O variables
  int print_freq;
  char buffer[255];
  char *pch;
  FILE *fp_ires;
  FILE *fp_vm;
  FILE *fp_cai;

  bool is_eligible_AP, is_local;
  double qnet;
  double t_ca_peak, ca_amp50, cad50_curr, cad50_prev, ca_amp90, cad90_curr, cad90_prev;
  double t_depol, vm_repol30, vm_repol90, vm_repol50;
  
  int pace;

  // SUNDIALs variables
  int cvode_retval;
  double tnext, tcurr;
  void *cvode_mem;
  N_Vector states_vec;

  int idx;

  if( p_cell == NULL ){
    p_cell = new Ohara_Rudy_2011();
    p_cell->initConsts();
    p_cell->CONSTANTS[BCL] = p_param->bcl_init;
    is_local = true;
  }
  else{
    printf("Cell has been created!! Skipped initialization....\n");
    printf("Scenario: %s\n", scenario_name);
    is_local = false;
  }

  tnext = p_param->dt;
  tcurr = 0.0;
  print_freq = (1. / p_param->dt) * p_param->dt_write;

  // Create CVODE solver
  cvode_mem = CVodeCreate( CV_BDF,CV_NEWTON );
  // Create the states vector based on the STATES array
  states_vec = N_VMake_Serial( p_cell->states_size, p_cell->STATES );
  // Initalize CVODE solver
  CVodeInit( cvode_mem, rhs_fn, tcurr, states_vec );
  // Give p_cell as CVode User Data
  CVodeSetUserData( cvode_mem, p_cell );
  // Set up the maximum step size (dt)
  CVodeSetMaxStep( cvode_mem,  p_param->dt );
  // Set up the linear solver type
  CVDense( cvode_mem, p_cell->states_size );
  // Set up the numerical error tolerances
  CVodeSStolerances( cvode_mem, 1.0e-6, 1.0e-6 );

  if( scenario_name == NULL || strlen(scenario_name) == 0 ) snprintf(buffer, sizeof(buffer), "result/vmcheck.plt");
  else snprintf(buffer, sizeof(buffer), "result/vmcheck_%s.plt", scenario_name);
  fp_vm = fopen( buffer, "w");

  if( scenario_name == NULL || strlen(scenario_name) == 0 ) snprintf(buffer, sizeof(buffer), "result/vmcheck.plt");
  else snprintf(buffer, sizeof(buffer), "result/cai_%s.plt", scenario_name);
  fp_cai = fopen( buffer, "w");


  if( scenario_name == NULL || strlen(scenario_name) == 0 ) snprintf(buffer, sizeof(buffer), "result/ires.plt");
  else snprintf(buffer, sizeof(buffer), "result/ires_%s.plt", scenario_name);
  fp_ires = fopen( buffer, "w");


  fprintf(fp_vm, "%s %s\n", "TIME", "Vm");
  fprintf(fp_ires, "%s %s %s %s %s %s %s\n","TIME", "INaL","ICaL","Ito","IKr","IKs","IK1");


  time_t begin = time(NULL);
  int iout = 0;
  int imax = (p_cell->CONSTANTS[BCL] * p_param->num_pace1) / p_param->dt;
  cipa_result.init( p_cell->STATES[V], p_cell->STATES[cai] );
  qnet = 0.;
  pace = 0;
  is_eligible_AP = false;

  while( iout <= imax ){

    // solve ODE
    cvode_retval = CVode( cvode_mem, tnext, states_vec, &tcurr, CV_NORMAL  );
    // get the RATES array
    p_cell->computeRates(tnext, p_cell->CONSTANTS, p_cell->RATES, N_VGetArrayPointer_Serial(states_vec), p_cell->ALGEBRAIC);
    if( cvode_retval == CV_SUCCESS ){
      iout++;
      tnext += p_param->dt;
    }
    else{
      if( scenario_name == NULL || strlen(scenario_name) == 0 ) fprintf(stderr, "CVode error!!\n");
      else fprintf(stderr, "CVode error at scenario %s\n", scenario_name);
      break;
    }
    // calculate qnet and area under curve
    qnet += (p_cell->ALGEBRAIC[INaL]+p_cell->ALGEBRAIC[ICaL]+p_cell->ALGEBRAIC[Ito]+p_cell->ALGEBRAIC[IKr]+p_cell->ALGEBRAIC[IKs]+p_cell->ALGEBRAIC[IK1])*p_param->dt;

    if( pace >= p_param->num_pace1-250 && iout % print_freq == 0 ) {
      temp_result.cai_data.insert( std::pair<double, double> (tcurr, p_cell->STATES[cai]) );
      temp_result.vm_data.insert( std::pair<double, double> (tcurr, p_cell->STATES[V]) );
      snprintf( buffer, sizeof(buffer), "%lf %lf %lf %lf %lf %lf", p_cell->ALGEBRAIC[INaL], p_cell->ALGEBRAIC[ICaL], p_cell->ALGEBRAIC[Ito], p_cell->ALGEBRAIC[IKr], p_cell->ALGEBRAIC[IKs], p_cell->ALGEBRAIC[IK1] );
      temp_result.ires_data.insert( std::pair<double, string> (tcurr, string(buffer)) );
    }


    // executed when entering new pace
    if( iout > 0 && iout % (int)( p_cell->CONSTANTS[BCL] / p_param->dt ) == 0 ) {
      if( pace >= p_param->num_pace1-249 ){

        // get the cad50 and cad90
        for(std::multimap<double, double>::iterator itrmap = temp_result.cai_data.begin(); itrmap != temp_result.cai_data.end() ; itrmap++ ){
          // before the peak calcium
          if( itrmap->first < t_ca_peak ){
            if( itrmap->second < ca_amp50 ) cad50_prev = itrmap->first;
            if( itrmap->second < ca_amp90 ) cad90_prev = itrmap->first;
          }
          // after the peak calcium
          else{
            if( itrmap->second > ca_amp50 ) cad50_curr = itrmap->first;
            if( itrmap->second > ca_amp90 ) cad90_curr = itrmap->first;
          }
        }
        temp_result.cad50 = cad50_curr - cad50_prev;
        temp_result.cad90 = cad90_curr - cad90_prev;
        temp_result.qnet = qnet;
        temp_result.vm_dia = p_cell->STATES[V];
        temp_result.ca_dia = p_cell->STATES[cai];
        // replace result with steeper repolarization AP
        if( temp_result.dvmdt_repol > cipa_result.dvmdt_repol ) cipa_result = temp_result;
      }

      pace++;
      qnet = 0.;

      if( pace >= p_param->num_pace1-250 ) {
        temp_result.init( p_cell->STATES[V], p_cell->STATES[cai] );
        t_ca_peak = tcurr;
        t_depol = (p_cell->CONSTANTS[BCL]*pace)+p_cell->CONSTANTS[stim_start];
        is_eligible_AP = false;
      }

    }

    // entering the last 250 paces
    if( p_param->num_pace1 > 250 && pace >= p_param->num_pace1-250 ){
      // get maximum dvmdt
      if( temp_result.dvmdt_max < p_cell->RATES[V] ) temp_result.dvmdt_max = p_cell->RATES[V];
      // get the peak Vm 6 secs after depolarization (when Na channel just closed after bursting)
      if( iout % (int)(((p_cell->CONSTANTS[BCL]*pace)+(p_cell->CONSTANTS[stim_start]+6.)) / p_param->dt) == 0) {
        temp_result.vm_peak = p_cell->STATES[V];
        // only finding the vm_repol30, vm_repol_50 and vm_repol90 if the peak Vm more than 0 mV
        if( temp_result.vm_peak > 0. ){
          vm_repol30 = temp_result.vm_peak - (0.3 * (temp_result.vm_peak - temp_result.vm_valley));
          vm_repol50 = temp_result.vm_peak - (0.5 * (temp_result.vm_peak - temp_result.vm_valley));
          vm_repol90 = temp_result.vm_peak - (0.9 * (temp_result.vm_peak - temp_result.vm_valley));
          is_eligible_AP = true;
        }
      }
      // these operations will be executed if it's eligible AP and done after the beginning of repolarization
      else if( is_eligible_AP && tcurr > (p_cell->CONSTANTS[BCL]*pace)+(p_cell->CONSTANTS[stim_start]+6.) ){
        // get the steepest dvmdt_repol around vm_repol30 and vm_repol90 if AP shape is eligible
        if( vm_repol30 >= p_cell->STATES[V] && p_cell->STATES[V] >= vm_repol90 && temp_result.dvmdt_repol < p_cell->RATES[V] ){
            temp_result.dvmdt_repol = p_cell->RATES[V];
        }
        // get the APD50
        if( vm_repol50 > p_cell->STATES[V] && p_cell->STATES[V] > vm_repol50-2 )
          temp_result.apd50 = tcurr - t_depol;
        // get the APD90
        if( vm_repol90 > p_cell->STATES[V] && p_cell->STATES[V] > vm_repol90-2 )
          temp_result.apd90 = tcurr - t_depol;
        // get the peak calcium and amplitude 50% and 90%
        if( temp_result.ca_peak < p_cell->STATES[cai] ){
          temp_result.ca_peak = p_cell->STATES[cai];
          ca_amp50 = temp_result.ca_peak - (0.5 * (temp_result.ca_peak - temp_result.ca_valley));
          ca_amp90 = temp_result.ca_peak - (0.9 * (temp_result.ca_peak - temp_result.ca_valley));
          t_ca_peak = tcurr;
        }
      }

    }

  }

  // print last result
  if( cvode_retval == CV_SUCCESS ){
    if(p_param->is_print_graph == 1 ) {
      for(std::multimap<double, double>::iterator itrmap = cipa_result.cai_data.begin(); itrmap != cipa_result.cai_data.end() ; itrmap++ ){
        fprintf(fp_cai, "%lf %lf\n", itrmap->first, itrmap->second);
      }
      for(std::multimap<double, double>::iterator itrmap = cipa_result.vm_data.begin(); itrmap != cipa_result.vm_data.end() ; itrmap++ ){
        fprintf(fp_vm, "%lf %lf\n", itrmap->first, itrmap->second);
      }
      for(std::multimap<double, string>::iterator itrmap = cipa_result.ires_data.begin(); itrmap != cipa_result.ires_data.end() ; itrmap++ ){
        fprintf(fp_ires, "%lf %s\n", itrmap->first, (itrmap->second).c_str());
      }
    }
  }
  else{
    if( scenario_name == NULL || strlen(scenario_name) == 0 ) fprintf(stderr, "Simulation qnet_bench got error!!!\n");
    else fprintf(stderr, "Scenario %s got error!!!\n", scenario_name);
  }


  // Memory Cleanup
  fclose(fp_ires);
  fclose(fp_vm);
  fclose(fp_cai);
  N_VDestroy(states_vec);
  CVodeFree(&cvode_mem);
  if(is_local == true) delete p_cell;
  if(is_local == true) free(p_param);

  return cipa_result;

}
