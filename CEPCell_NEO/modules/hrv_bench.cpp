// Code located in this project directory
#include "hrv_bench.hpp"
#include "commons.hpp"
#include "globals.hpp"
#include "../cellmodels/Ohara_Rudy_2011.hpp"

// Standard libs
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <map>
#include <string>
#include <vector>
#include <mpi.h>

// SUNDIALs libs
#include <cvode/cvode.h>
#include <cvode/cvode_dense.h>
#include <nvector/nvector_serial.h>

// POSIX libs
#include <dirent.h>
#include <magic.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

using std::copy;
using std::map;
using std::multimap;
using std::pair;
using std::string;
using std::vector;

#define DEBUG

void hrv_bench( int argc, char **argv, param_t *p_param )
{
  // CellML object pointer
  patch_clamp *p_cell;

  // Supporting variables
  char buffer[150] = {'\0'};
  char *pch;
  int idx, jdx, group_id, sample_id, sim_counter;
  int print_freq;
  string foldername;
  vector<double> temp_data;

  // HRV-related variables
  FILE *fp_hrv;
  vector<double> hrv_data;

  // CiPA related variables
  cipa_t cipa_result, temp_result;
  double inal_auc_ctl, ical_auc_ctl, inal_auc_drug, ical_auc_drug, qnet, qinward, inal_auc, ical_auc;
  double hill_data_arr[14];
  FILE *fp_vm = NULL, *fp_ca = NULL;
  FILE *fp_hill, *fp_qni, *fp_debug;
  int pace, steepest_pace;
  vector< vector<double> > hill_data;  
  vector<double> concs;

  // AP related variables
  bool is_vm_peak_found;
  bool is_eligible_AP;
  double t_depol, vm_repol30, vm_repol90, vm_repol50;
  FILE *fp_ap_profile;

  // Ca related variables
  double ca_amp50, ca_amp90;
  double cad50_prev, cad90_prev, cad50_curr, cad90_curr;
  double t_ca_peak;
  FILE *fp_ca_profile, *fp_ca_debug;
  multimap<double, double> cai_data;

  // POSIX variable
  struct stat st = {0};
  
  // SUNDIALs variables
  bool cvode_firsttime;
  int cvode_retval, iout, imax;
  double tnext, tcurr = 0.0;
  void *cvode_mem;
  N_Vector states_vec;


  // Open the file and process the data inside it.
  fp_hrv = fopen(p_param->hrv_files, "r" );
  if( fp_hrv == NULL ) {fprintf(stderr, "Cannot open file %s\n", buffer ); return;}
  while( fgets(buffer, sizeof(buffer), fp_hrv) != NULL ) hrv_data.push_back(strtod(buffer, NULL));
  hrv_data.push_back(1000.);
  fclose(fp_hrv);

  // setup concentration vector
  concs = get_vector_from_string( std::string(p_param->concs) );
  if( concs.size() > 4+1 && mympi::rank == 0 ){  // +1 because we include concentration 0
    fprintf(stderr, "Only 4 concentrations is allowed in this program! Exiting...\n");
    return;
  }

  // make folders from the concentration
  // put MPI_BARRIER to avoid sync error
  if(mympi::rank == 0){
    for( idx = 0; idx < concs.size(); idx++ ){
      snprintf( buffer, sizeof(buffer), "result/%.2lf", concs[idx] );
      if(stat("result/", &st) == 0) mkdir(buffer, 0775);
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);

  // setup hill data vector
  fp_hill = fopen(p_param->hill_file, "r");
  if (fp_hill == NULL) {
    fprintf(stderr, "Cannot open Hill file %s\n", p_param->hill_file);
    return;
  }
  fgets(buffer, sizeof(buffer), fp_hill); // skip header
  temp_data.clear();
  while( fgets(buffer, sizeof(buffer), fp_hill) != NULL )
  { // begin data loop

    temp_data = get_vector_from_string( std::string(buffer) );
    hill_data.push_back( temp_data );

  } // end data loop
  fclose(fp_hill);
  if( hill_data.size() > 2000 && mympi::rank == 0 ){  // sample cannot exceed 2000
    fprintf(stderr, "Hill samples are more than 2000! Exiting...\n");
    return;
  }

  // variable initialization
  p_cell = new Ohara_Rudy_2011();
  cvode_mem = CVodeCreate( CV_BDF,CV_NEWTON );
  cvode_firsttime = true;
  print_freq = (1. / p_param->dt) * p_param->dt_write;


  // main loop
  double begin = MPI_Wtime();

  int sample_group = ( hill_data.size()%mympi::size != 0 ) ? (hill_data.size()/mympi::size)+1 : (hill_data.size()/mympi::size);
  for( group_id=0, sim_counter=mympi::rank;
       group_id < sample_group && sim_counter < hill_data.size();
       group_id++, sim_counter += mympi::size )
  { // begin hill_data loop

    sample_id = (group_id*mympi::size)+mympi::rank;
    copy( hill_data[sample_id].begin(), hill_data[sample_id].end(), hill_data_arr );
    printf("Sample_ID:%d  Count:%d Sample Group:%d Rank:%d Group_ID:%d\nData: ", sample_id, hill_data.size() ,sample_group, mympi::rank, group_id );
    for( std::vector<double>::iterator itrvec = hill_data[sample_id].begin(); itrvec != hill_data[sample_id].end(); itrvec++ ){
      printf("%lf ", *itrvec);
    }
    printf("\n");

    copy( hill_data[sample_id].begin(), hill_data[sample_id].end(), hill_data_arr );
  
    for( idx = 0; idx < concs.size(); idx++ ){ // begin concentration loop
      // init cell model
      p_cell->initConsts(p_param->celltype, concs[idx], hill_data_arr, true);
      p_cell->CONSTANTS[amp] *= p_param->stim_amp;

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

      // FILE pointer naming and initialization
      if( p_param->is_print_graph == 1 ){
        snprintf(buffer, sizeof(buffer), "result/%.2lf/%s_%.2lf_vmcheck_smp%d.plt", concs[idx], p_param->drug_name, concs[idx], sample_id );
        fp_vm = fopen( buffer, "w" );
        snprintf(buffer, sizeof(buffer), "result/%.2lf/%s_%.2lf_ca_i_smp%d.plt", concs[idx], p_param->drug_name, concs[idx], sample_id );
        fp_ca = fopen( buffer, "w" );
      }
      snprintf(buffer, sizeof(buffer), "result/%.2lf/%s_%.2lf_ap_profile_proc%d.plt", concs[idx], p_param->drug_name, concs[idx], mympi::rank );
      fp_ap_profile = fopen( buffer, "a" );
      snprintf(buffer, sizeof(buffer), "result/%.2lf/%s_%.2lf_ca_profile_proc%d.plt", concs[idx], p_param->drug_name, concs[idx], mympi::rank );
      fp_ca_profile = fopen( buffer, "a" );
      snprintf(buffer, sizeof(buffer), "result/%.2lf/%s_%.2lf_qni_proc%d.plt", concs[idx], p_param->drug_name, concs[idx], mympi::rank );
      fp_qni = fopen( buffer, "a" );

      // Write header for all of the result file
      if( group_id == 0 )
      {
        if( p_param->is_print_graph == 1 ){
          fprintf( fp_vm, "%s %s\n", "TIME", "Vm" );
          fprintf( fp_ca, "%s %s\n", "TIME", "cai" );
        }
        fprintf( fp_ap_profile, "%s %s %s %s %s %s %s %s %s %s\n",
                 "Sample_ID", "Dvm/Dt_Repol", "Max_Dvm/Dt", "Vm_Peak", "Vm_Resting","APD90", "APD50", "APDTri", "Steepest_Pace", "BCL");
        fprintf( fp_ca_profile, "%s %s %s %s %s %s %s\n",
                 "Sample_ID", "Ca_Peak", "Ca_Diastole", "CaD90", "CaD50","Catri", "Steepest_Pace");
        fprintf( fp_qni, "%s %s %s\n","Sample_ID", "Qnet", "Qinward");
      }

      // Initialization
      pace = 0;
      cipa_result.init( p_cell->STATES[V], p_cell->STATES[cai] );
      qnet = inal_auc = ical_auc = 0.;

      for( std::vector<double>::iterator itrvec = hrv_data.begin(); itrvec != hrv_data.end(); itrvec++ ){ // begin BCL loop

        // set the BCL value
        p_cell->CONSTANTS[BCL] = *itrvec;

        // mandatory initializations between paces (without this, expect some problems will happen)
        // THIS IS IMPORTANT!!!! PLEASE DON'T REMOVE THESE!!!
        qnet = inal_auc = ical_auc = 0.0;
        iout = 0;
        imax = p_cell->CONSTANTS[BCL] / p_param->dt;
        is_vm_peak_found = false;
        is_eligible_AP = false;
        t_ca_peak = tcurr;
        if(pace > 0){
          tcurr -= p_param->dt;
          tnext -= p_param->dt;
        }
        if(pace >= hrv_data.size()-250) temp_result.init( p_cell->STATES[V], p_cell->STATES[cai] );
        
        while(iout <= imax){ // begin one pace loop
          // solve ODE
          cvode_retval = CVode( cvode_mem, tnext, states_vec, &tcurr, CV_NORMAL  );
          // get the RATES array
          p_cell->computeRates(tnext, p_cell->CONSTANTS, p_cell->RATES, N_VGetArrayPointer_Serial(states_vec), p_cell->ALGEBRAIC);
          if( cvode_retval == CV_SUCCESS ){
            iout++;
            tnext += p_param->dt;
          }
          else{
            fprintf(stderr, "CVode calculation error at pace %d BCL %lf concentration %.2lf msec!!!\n", pace+1, p_cell->CONSTANTS[BCL], concs[idx]);
            break;
          }
          // calculation for CiPA parameters
          qnet += (p_cell->ALGEBRAIC[INaL]+p_cell->ALGEBRAIC[ICaL]+p_cell->ALGEBRAIC[Ito]+p_cell->ALGEBRAIC[IKr]+p_cell->ALGEBRAIC[IKs]+p_cell->ALGEBRAIC[IK1])*p_param->dt;
          inal_auc += (p_cell->ALGEBRAIC[INaL])*p_param->dt;
          ical_auc += (p_cell->ALGEBRAIC[ICaL])*p_param->dt;

          // save output for each pace
          if( pace >= hrv_data.size()-250 && iout % print_freq == 0 ) {
            temp_result.cai_data.insert( std::pair<double, double> (tcurr, p_cell->STATES[cai]) );
            temp_result.vm_data.insert( std::pair<double, double> (tcurr, p_cell->STATES[V]) );
          }

          if( pace >= hrv_data.size()-250 ){ // begin of last 250 paces operations
            // get maximum dvmdt
            if( temp_result.dvmdt_max < p_cell->RATES[V] ) temp_result.dvmdt_max = p_cell->RATES[V];

            // get the peak Vm 6 secs after depolarization (when Na channel just closed after bursting)
            if( iout % (int)((p_cell->CONSTANTS[stim_start]+6.) / p_param->dt) == 0 && is_vm_peak_found == false) {
              temp_result.vm_peak = p_cell->STATES[V];
              t_depol = tcurr;
              // only finding the vm_repol30, vm_repol_50 and vm_repol90 if the peak Vm more than 0 mV
              if( temp_result.vm_peak > 0. ){
                vm_repol30 = temp_result.vm_peak - (0.3 * (temp_result.vm_peak - temp_result.vm_valley));
                vm_repol50 = temp_result.vm_peak - (0.5 * (temp_result.vm_peak - temp_result.vm_valley));
                vm_repol90 = temp_result.vm_peak - (0.9 * (temp_result.vm_peak - temp_result.vm_valley));
                is_eligible_AP = true;
              }
              is_vm_peak_found = true;
            }
            // these operations will be executed if it's eligible AP and done after the beginning of repolarization
            else if( is_eligible_AP && tcurr > p_cell->CONSTANTS[cumm_BCL]+(p_cell->CONSTANTS[stim_start]+6.) ){

              //printf("DeTECTED MAX REPOL: %lf CURR REPOL: %lf VM30: %lf VM90: %lf\n", temp_result.dvmdt_repol, p_cell->RATES[V], vm_repol30, vm_repol90);
              // get the steepest dvmdt_repol around vm_repol30 and vm_repol90 if AP shape is eligible
              if( vm_repol30 >= p_cell->STATES[V] && p_cell->STATES[V] >= vm_repol90 && temp_result.dvmdt_repol < p_cell->RATES[V] ){
                temp_result.dvmdt_repol = p_cell->RATES[V];
              }
              // get the APD50
              if( vm_repol50 > p_cell->STATES[V] && p_cell->STATES[V] > vm_repol50-2 )temp_result.apd50 = tcurr - t_depol;
              // get the APD90
              if( vm_repol90 > p_cell->STATES[V] && p_cell->STATES[V] > vm_repol90-2 )temp_result.apd90 = tcurr - t_depol;
              // get the peak calcium and amplitude 50% and 90%
              if( temp_result.ca_peak < p_cell->STATES[cai] ){
                temp_result.ca_peak = p_cell->STATES[cai];
                ca_amp50 = temp_result.ca_peak - (0.5 * (temp_result.ca_peak - temp_result.ca_valley));
                ca_amp90 = temp_result.ca_peak - (0.9 * (temp_result.ca_peak - temp_result.ca_valley));
                t_ca_peak = tcurr;
              }
            }

          } // end of last 250 paces operations

        } // end one pace loop 

        if( pace >= hrv_data.size()-250 ){
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
          temp_result.inal_auc = inal_auc;
          temp_result.ical_auc = ical_auc;
          temp_result.vm_dia = p_cell->STATES[V];
          temp_result.ca_dia = p_cell->STATES[cai];
          temp_result.bcl = p_cell->CONSTANTS[BCL];

          if( temp_result.dvmdt_repol > cipa_result.dvmdt_repol ) {
            steepest_pace = pace;
            cipa_result = temp_result;
          }

        }

        pace++;
        p_cell->CONSTANTS[cumm_BCL] += p_cell->CONSTANTS[BCL];

        // initialize cipa_result with the first temp_result
        // (only executed once, to avoid error when STATES[V] < 0)
        // (reported by KITOX at Monday 16th March 2021)
        if( pace == p_param->num_pace1-249 ) {
          cipa_result = temp_result;
        }


      } // end BCL loop

      // print last result
      if( cvode_retval == CV_SUCCESS ){

        if(p_param->is_print_graph == 1 ) {
          for(std::multimap<double, double>::iterator itrmap = cipa_result.cai_data.begin(); itrmap != cipa_result.cai_data.end() ; itrmap++ ){
            fprintf(fp_ca, "%lf %lf\n", itrmap->first, itrmap->second);
          }
          for(std::multimap<double, double>::iterator itrmap = cipa_result.vm_data.begin(); itrmap != cipa_result.vm_data.end() ; itrmap++ ){
            fprintf(fp_vm, "%lf %lf\n", itrmap->first, itrmap->second);
          }
        }

        if( cipa_result.apd90 < cipa_result.apd50 ) cipa_result.apd90 = cipa_result.bcl;

        fprintf( fp_ap_profile, "%d %lf %lf %lf %lf %lf %lf %lf %d %d\n",
               sample_id,
               cipa_result.dvmdt_repol,
               cipa_result.dvmdt_max,
               cipa_result.vm_peak,
               cipa_result.vm_dia,
               cipa_result.apd90,
               cipa_result.apd50,
               cipa_result.apd90-cipa_result.apd50,
               steepest_pace,
               cipa_result.bcl);
        fprintf( fp_ca_profile, "%d %lf %lf %lf %lf %lf %d\n",
               sample_id,
               cipa_result.ca_peak,
               cipa_result.ca_dia,
               cipa_result.cad90,
               cipa_result.cad50,
               cipa_result.cad90-cipa_result.cad50,
               steepest_pace );
        if((int)ceil(concs[idx]) == 0){
          inal_auc_ctl = cipa_result.inal_auc;
          ical_auc_ctl = cipa_result.ical_auc;
          fprintf( fp_qni, "%d %lf %lf\n", sample_id, (cipa_result.qnet)/1000.0, 0. );
        }
        else{
          inal_auc_drug = cipa_result.inal_auc;
          ical_auc_drug = cipa_result.ical_auc;
          fprintf( fp_qni, "%d %lf %lf\n", sample_id, (cipa_result.qnet)/1000.0, ( (inal_auc_drug/inal_auc_ctl) + (ical_auc_drug/ical_auc_ctl) ) * 0.5 );
        }
      }
      else{
        fprintf( fp_ap_profile, "%d CVODE_ERR\n", sample_id);
        fprintf( fp_ca_profile, "%d CVODE_ERR\n", sample_id);
        fprintf( fp_qni, "%d CVODE_ERR\n", sample_id);
      }


      // memory cleanup
      if( p_param->is_print_graph == 1 ){
        fclose( fp_vm );
        fclose( fp_ca );
      }
      fclose( fp_ap_profile );
      fclose( fp_ca_profile );
      fclose( fp_qni );

      N_VDestroy_Serial(states_vec);
    } // end concentration loop

  } // end hill_data loop

  double end = MPI_Wtime();

  // wait for all nodes finished their job
  MPI_Barrier(MPI_COMM_WORLD);
/*
  if( mympi::rank == 0 ){
    // setup timestamp
    char buff_time[15] = {'\0'};
    FILE *fp_perf_log;
    struct tm* tm_info;
    time_t timer;
    timer = time(NULL);
    tm_info = localtime(&timer);
    strftime( buff_time, 15, "%Y%m%d%H%M%S", tm_info );

    // zip result folder
    snprintf( buffer, sizeof(buffer), "result_%s.zip", buff_time );
    create_zip(buffer, "result");

    // delete the content inside result folder
    remove_dir_content("result");

    // copy the zip into the result folder
    char buff_new[150] = {'\0'};
    snprintf( buffer, sizeof(buffer), "result_%s.zip", buff_time );
    snprintf( buff_new, sizeof(buff_new), "result/result_%s.zip", buff_time );
    if( rename( buffer, buff_new ) == 0 ){
      printf("Zip file %s moved succesfully to %s!!\n", buffer, buff_new);
    }
    else{
      printf("Error moving file %s to %s!!\n", buffer, buff_new);
      return;
    }

    // make a timestamp log
    snprintf( buffer, sizeof(buffer), "result/performance_%s.log", buff_time );
    fp_perf_log = fopen( buffer, "w" );
    fprintf( fp_perf_log, "BCL_Count CPU Concentration Time(Minutes)\n" );
    fprintf( fp_perf_log, "%d %d %d %.5lf \n", hrv_data.size(), mympi::size, concs.size(), (end-begin)/60. );
    fclose( fp_perf_log );
    printf("Computation Time: %lf minutes\n", (end-begin)/60.);
  }
*/
  CVodeFree(&cvode_mem);
  delete p_cell;
  free(p_param);

}
