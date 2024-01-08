#include "cipa_drug_last_bench.hpp"
#include "commons.hpp"
#include "globals.hpp"
#include "../cellmodels/Ohara_Rudy_2011.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <cvode/cvode.h>
#include <cvode/cvode_dense.h>
#include <map>
#include <mpi.h>
#include <nvector/nvector_serial.h>
#include <string>
#include <vector>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

using std::copy;
using std::multimap;
using std::vector;

void cipa_drug_last_bench( int argc, char **argv, param_t *p_param, bool is_dutta )
{
  // CellML object pointer
  patch_clamp *p_cell;

  // Supporting variables
  char buffer[150] = {'\0'};
  char *pch;
  double states_init[50];
  FILE *fp_hill;
  FILE *fp_states;
  int idx, jdx, group_id, sample_id, sim_counter, iout, imax;
  int print_freq;
  vector<double> temp_data;

  // CiPA related variables
  cipa_t cipa_result, temp_result;
  double INaL_auc_control, ICaL_auc_control, INaL_auc_drug, ICaL_auc_drug, qnet, inal_auc, ical_auc;
  double hill_data_arr[14];
  FILE* fp_vm = NULL, * fp_ca = NULL, * fp_ires = NULL;
  FILE *fp_qni, *fp_ap_profile, *fp_ca_profile;
  int pace, steepest_pace;
  vector< vector<double> > hill_data;
  vector<double> concs;

  // Calcium related variables
  double t_ca_peak, ca_amp50, cad50_curr, cad50_prev, ca_amp90, cad90_curr, cad90_prev;

  // AP related variables
  bool is_eligible_AP;
  double t_depol, vm_repol30, vm_repol90, vm_repol50;

  // SUNDIALs variables
  bool cvode_firsttime;
  int cvode_retval;
  double tnext, tcurr = 0.0;
  void *cvode_mem;
  N_Vector states_vec;

  // POSIX variables
  struct stat st = {0};

  // give limitations to the number of pacing
  assert( (p_param->num_pace1 > 2 && p_param->num_pace1 <= 200) && "Pace should be more than 2 and less than 1000");

  // setup concentration vector
  concs = get_vector_from_string( std::string(p_param->concs) );
  assert( concs.size() > 0 && concs.size() <= 4+1 && "Number of concentration should be less than or equal to 5");  // +1 because we include concentration 0

  // make folders from the concentration
  for( idx = 0; idx < concs.size(); idx++ ){
    snprintf( buffer, sizeof(buffer), "result/%.2lf", concs[idx] );
    if(mympi::rank == 0 && stat("result/", &st) == 0) mkdir(buffer, 0775);
  }

  // setup hill data vector
  fp_hill = fopen(p_param->hill_file, "r");
  assert (fp_hill != NULL && "Problem with hill input file");
  fgets(buffer, sizeof(buffer), fp_hill); // skip header
  while( fgets(buffer, sizeof(buffer), fp_hill) != NULL )
  { // begin data loop

    temp_data = get_vector_from_string( std::string(buffer) );
    hill_data.push_back( temp_data );

  } // end data loop
  fclose(fp_hill);
  assert( hill_data.size() >0 && hill_data.size() <= 2000 && "Hill samples should be less than or equal 2000" );

  // put initial steady state
  if( p_param->is_using_output > 0 ){
    fp_states = fopen("output_orudy.dat", "r");
    assert( fp_states != NULL && "Initial steady state file error");
    idx = 0;
    while(fgets( buffer, sizeof(buffer), fp_states) != NULL){
      states_init[idx] = strtod(buffer, NULL);
      idx++;
    }
    fclose(fp_states);
    printf("Using output steady states...\n");
  }

  p_cell = new Ohara_Rudy_2011();
  assert( p_cell != NULL && "Error when creating cell object" );

  cvode_mem = CVodeCreate( CV_BDF,CV_NEWTON );
  assert( cvode_mem != NULL && "Error when allocating CVode" );

  cvode_firsttime = true;
  print_freq = (1. / p_param->dt) * p_param->dt_write;

  // begin sample loop
  double begin = MPI_Wtime();

  int sample_group = ( hill_data.size()%mympi::size != 0 ) ? (hill_data.size()/mympi::size)+1 : (hill_data.size()/mympi::size);
  for( group_id=0, sim_counter=mympi::rank; 
       group_id < sample_group && sim_counter < hill_data.size();
       group_id++, sim_counter += mympi::size)
  { 
    sample_id = (group_id*mympi::size)+mympi::rank;
    printf("Sample_ID:%d  Count:%d Sample Group:%d Rank:%d Group_ID:%d\nData: ", sample_id, hill_data.size(), sample_group, mympi::rank, group_id );
    for( std::vector<double>::iterator itrvec = hill_data[sample_id].begin(); itrvec != hill_data[sample_id].end(); itrvec++ ){
      printf("%lf ", *itrvec);
    }
    printf("\n");

    copy( hill_data[sample_id].begin(), hill_data[sample_id].end(), hill_data_arr );

    for( idx = 0; idx < concs.size(); idx++ )
    { // begin concs loop

      // apply some cell initialization
      p_cell->initConsts(p_param->celltype, concs[idx], hill_data_arr, true);
      p_cell->CONSTANTS[BCL] = p_param->bcl_init;
      p_cell->CONSTANTS[amp] *= p_param->stim_amp;

      // replace states with the initial steady state, if possible
      if( p_param->is_using_output > 0 )  memcpy( p_cell->STATES, states_init, sizeof(p_cell->STATES) );

      // Time and CVode setup
      tcurr = 0.;
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

      if( p_param->is_print_graph == 1 ){
        snprintf(buffer, sizeof(buffer), "result/%.2lf/%s_%.2lf_vmcheck_smp%d.plt", concs[idx], p_param->drug_name, concs[idx], sample_id );
        fp_vm = fopen( buffer, "w" );
        assert(fp_vm != NULL);
        snprintf(buffer, sizeof(buffer), "result/%.2lf/%s_%.2lf_ca_i_smp%d.plt", concs[idx], p_param->drug_name, concs[idx], sample_id );
        fp_ca = fopen( buffer, "w" );
        snprintf(buffer, sizeof(buffer), "result/%.2lf/%s_%.2lf_ires_smp%d.plt", concs[idx], p_param->drug_name, concs[idx], sample_id );
        fp_ires = fopen( buffer, "w" );
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
          fprintf( fp_ires, "%s %s %s %s %s %s %s\n", "TIME", "INaL", "ICaL", "Ito", "IKr", "IKs", "IK1" );
        }
        fprintf( fp_ap_profile, "%s %s %s %s %s %s %s %s %s\n",
                 "Sample_ID", "Dvm/Dt_Repol","Max_Dvm/Dt", "Vm_Peak", "Vm_Resting","APD90", "APD50", "APDTri", "Steepest_Pace");
        fprintf( fp_ca_profile, "%s %s %s %s %s %s %s\n",
                 "Sample_ID", "Ca_Peak", "Ca_Diastole", "CaD90", "CaD50","Catri", "Steepest_Pace");
        fprintf( fp_qni, "%s %s %s\n","Sample_ID", "Qnet", "Qinward");
      }

      // initialize some variables
      cipa_result.init( p_cell->STATES[V], p_cell->STATES[cai] );
      qnet = inal_auc = ical_auc = 0.;
      iout = 0;
      imax = ( p_param->num_pace1 *  p_cell->CONSTANTS[BCL] ) / p_param->dt;
      pace = steepest_pace = 0;
      is_eligible_AP = false;
      t_depol = (p_cell->CONSTANTS[BCL]*pace)+p_cell->CONSTANTS[stim_start];

      while(iout < imax){ // begin paces loop

        // solve ODE
        cvode_retval = CVode( cvode_mem, tnext, states_vec, &tcurr, CV_NORMAL  );
        // get the RATES array
        p_cell->computeRates(tnext, p_cell->CONSTANTS, p_cell->RATES, N_VGetArrayPointer_Serial(states_vec), p_cell->ALGEBRAIC);
        if( cvode_retval == CV_SUCCESS ){
          iout++;
          tnext += p_param->dt;
        }
        else{
          fprintf(stderr, "CVode error at sample_ID %d and concentration %.2lf at rank %d\n", sample_id, concs[idx], mympi::rank);
          break;
        }
        // calculate qnet and area under curve
        qnet += (p_cell->ALGEBRAIC[INaL]+p_cell->ALGEBRAIC[ICaL]+p_cell->ALGEBRAIC[Ito]+p_cell->ALGEBRAIC[IKr]+p_cell->ALGEBRAIC[IKs]+p_cell->ALGEBRAIC[IK1])*p_param->dt;
        inal_auc += p_cell->ALGEBRAIC[INaL]*p_param->dt;
        ical_auc += p_cell->ALGEBRAIC[ICaL]*p_param->dt;

        // save output for each pace
        if( pace >= p_param->num_pace1-2 && iout % print_freq == 0 ) {
          temp_result.cai_data.insert( std::pair<double, double> (tcurr, p_cell->STATES[cai]) );
          temp_result.vm_data.insert( std::pair<double, double> (tcurr, p_cell->STATES[V]) );
          snprintf( buffer, sizeof(buffer), "%lf %lf %lf %lf %lf %lf", p_cell->ALGEBRAIC[INaL], p_cell->ALGEBRAIC[ICaL], p_cell->ALGEBRAIC[Ito], p_cell->ALGEBRAIC[IKr], p_cell->ALGEBRAIC[IKs], p_cell->ALGEBRAIC[IK1] );
          temp_result.ires_data.insert( std::pair<double, string> (tcurr, string(buffer)) );
        }

        // executed when entering new pace
        if( iout > 0 && iout % (int)( p_cell->CONSTANTS[BCL] / p_param->dt ) == 0 ) {
          // executed at the end of AP pacing
          if( is_eligible_AP && pace >= p_param->num_pace1-1) {
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
            // replace result with steeper repolarization AP
            if( temp_result.dvmdt_repol > cipa_result.dvmdt_repol ) {
              steepest_pace = pace;
              cipa_result = temp_result;
            }
          }

          pace++;
          qnet = inal_auc = ical_auc = 0.;

          // executed at the beginning of AP pacing
          if( pace >= p_param->num_pace1-2 ) {
            //if(mympi::rank == 0)printf("Time new cycle: %lf with pace: %d\n", tcurr, pace);
            temp_result.init( p_cell->STATES[V], p_cell->STATES[cai] );
            t_ca_peak = tcurr;
            t_depol = (p_cell->CONSTANTS[BCL]*pace)+p_cell->CONSTANTS[stim_start];
            is_eligible_AP = false;
          }

        }

        // entering the last 2 paces
        if( p_param->num_pace1 > 2 && pace >= p_param->num_pace1-2 ){
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

        }// end of last 2 paces operations

      } // end paces loop

      // print last result
      if( cvode_retval == CV_SUCCESS ){

        if(p_param->is_print_graph == 1 ) {
          for(std::multimap<double, double>::iterator itrmap = cipa_result.cai_data.begin(); itrmap != cipa_result.cai_data.end() ; itrmap++ ){
            fprintf(fp_ca, "%lf %lf\n", itrmap->first, itrmap->second); 
          }
          for(std::multimap<double, double>::iterator itrmap = cipa_result.vm_data.begin(); itrmap != cipa_result.vm_data.end() ; itrmap++ ){
            fprintf(fp_vm, "%lf %lf\n", itrmap->first, itrmap->second); 
          }
          for(std::multimap<double, string>::iterator itrmap = cipa_result.ires_data.begin(); itrmap != cipa_result.ires_data.end() ; itrmap++ ){
            fprintf(fp_ires, "%lf %s\n", itrmap->first, (itrmap->second).c_str()); 
          }

        }

        fprintf( fp_ap_profile, "%d %lf %lf %lf %lf %lf %lf %lf %d\n",
               sample_id,
               cipa_result.dvmdt_repol,
               cipa_result.dvmdt_max,
               cipa_result.vm_peak,
               cipa_result.vm_dia,
               cipa_result.apd90,
               cipa_result.apd50,
               cipa_result.apd90-cipa_result.apd50,
               steepest_pace);
        fprintf( fp_ca_profile, "%d %lf %lf %lf %lf %lf %d\n",
               sample_id,
               cipa_result.ca_peak,
               cipa_result.ca_dia,
               cipa_result.cad90,
               cipa_result.cad50,
               cipa_result.cad90-cipa_result.cad50,
               steepest_pace );
        if((int)ceil(concs[idx]) == 0){
          INaL_auc_control = cipa_result.inal_auc;
          ICaL_auc_control = cipa_result.ical_auc;
          fprintf( fp_qni, "%d %lf %lf\n", sample_id, (cipa_result.qnet)/1000.0, 0.0 );
        }
        else{
          INaL_auc_drug = cipa_result.inal_auc;
          ICaL_auc_drug = cipa_result.ical_auc;
          fprintf( fp_qni, "%d %lf %lf\n", sample_id, (cipa_result.qnet)/1000.0, ( (INaL_auc_drug/INaL_auc_control) + (ICaL_auc_drug/ICaL_auc_control) ) * 0.5 );
        }

      }
      else{
        fprintf( fp_ap_profile, "%d ERR ERR ERR ERR ERR ERR\n", sample_id);
        fprintf( fp_ca_profile, "%d ERR ERR ERR ERR ERR\n", sample_id);
        fprintf( fp_qni, "%d ERR ERR\n", sample_id);
      }

      if( p_param->is_print_graph == 1 ){
        fclose( fp_vm );
        fclose( fp_ca );
        fclose( fp_ires );
      }
      fclose( fp_ap_profile );
      fclose( fp_ca_profile );
      fclose( fp_qni );

      N_VDestroy(states_vec);


    }  //end concs loop
  } // end sample loop
  double end = MPI_Wtime();

  // wait for all nodes finished their job
  MPI_Barrier(MPI_COMM_WORLD);

  // master rank will make a performance log and clean all used memory
  if( mympi::rank == 0 ){
    char buff_time[15] = {'\0'};
    FILE *fp_perf_log;
    struct tm* tm_info;
    time_t timer;
    
    timer = time(NULL);
    tm_info = localtime(&timer);
    strftime( buff_time, 15, "%Y%m%d%H%M%S", tm_info );
    snprintf( buffer, sizeof(buffer), "result/performance_%s.log", buff_time );

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

    fp_perf_log = fopen( buffer, "w" );
    fprintf( fp_perf_log, "Iteration BCL Timestep CPU Sample Pacing Time(Minutes)\n" );
    fprintf( fp_perf_log, "%lld %lf %lf %d %d %d %lf\n", (int)(hill_data.size()*concs.size()*p_param->num_pace1*p_cell->CONSTANTS[BCL]/p_param->dt), p_cell->CONSTANTS[BCL], p_param->dt, mympi::size, hill_data.size(), p_param->num_pace1, (end-begin)/60. );

    printf("Computation Time: %lf minutes\n", (end-begin)/60.);

    fclose( fp_perf_log );
  }

  CVodeFree(&cvode_mem);
  delete p_cell;
  free(p_param);

}
