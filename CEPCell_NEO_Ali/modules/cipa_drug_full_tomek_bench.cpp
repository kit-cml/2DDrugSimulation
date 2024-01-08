#include "cipa_drug_full_tomek_bench.hpp"
#include "commons.hpp"
#include "globals.hpp"
#include "../cellmodels/Tomek_model.hpp"

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

void cipa_drug_full_tomek_bench( int argc, char **argv, param_t *p_param )
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
  int print_freq, print_freq_ca;
  vector<double> temp_data;
  double t_repol90, t_repol50;


  // CiPA related variables
  cipa_t cipa_result, temp_result;
  double INaL_auc_control, ICaL_auc_control, INaL_auc_drug, ICaL_auc_drug, qnet, inal_auc, ical_auc;
  double hill_data_arr[14];
  FILE *fp_vm, *fp_ca,*fp_feature, *fp_output;
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
  //assert( (p_param->num_pace1 > 250 && p_param->num_pace1 <= 1000) && "Pace should be more than 250 and less than 2000");

  // setup concentration vector
  concs = get_vector_from_string( std::string(p_param->concs) );
  assert( concs.size() > 0 && concs.size() <= 4+1 && "Number of concentration should be less than or equal to 5");  // +1 because we include concentration 0

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
    fp_states = fopen("output_tomek.dat", "r");
    assert( fp_states != NULL && "Initial steady state file error");
    idx = 0;
    while(fgets( buffer, sizeof(buffer), fp_states) != NULL){
      states_init[idx] = strtod(buffer, NULL);
      idx++;
    }
    fclose(fp_states);
    printf("Using output steady states...\n");
  }

  p_cell = new Tomek_model();
  assert( p_cell != NULL && "Error when creating cell object" );

  cvode_mem = CVodeCreate( CV_BDF,CV_NEWTON );
  assert( cvode_mem != NULL && "Error when allocating CVode" );

  cvode_firsttime = true;
  print_freq = (1. / p_param->dt) * p_param->dt_write;
  print_freq_ca = (1. / p_param->dt) * 5;

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
      p_cell->initConsts(p_param->celltype, concs[idx], hill_data_arr);
      p_cell->CONSTANTS[BCL] = p_param->bcl_init;
      p_cell->CONSTANTS[amp] *= p_param->stim_amp;
      p_cell->CONSTANTS[GNa] *=  p_param->gna_scale;

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
        snprintf(buffer, sizeof(buffer), "result/%.2lf/%s_%.2lf_ca_i_smp%d.plt", concs[idx], p_param->drug_name, concs[idx], sample_id );
        fp_ca = fopen( buffer, "w" );
      }
      snprintf(buffer, sizeof(buffer), "result/%.2lf/%s_%.2lf_feature_smp%d.plt", concs[idx], p_param->drug_name, concs[idx], sample_id );
      fp_feature = fopen( buffer, "w" );
      snprintf(buffer, sizeof(buffer), "result/%.2lf/%s_%.2lf_output_smp%d.plt", concs[idx], p_param->drug_name, concs[idx], sample_id );
      fp_output = fopen( buffer, "w" );


      // Write header for all of the result file
      if( p_param->is_print_graph == 1 ){
        fprintf( fp_vm, "%s %s\n", "TIME", "Vm" );
        fprintf( fp_ca, "%s %s\n", "TIME", "cai" );
      }
      fprintf( fp_feature, "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n",
               "Pace", "Dvm/Dt_Repol", "Max_Dvm/Dt", "Vm_Peak", "Vm_Resting","APD90", "APD50", "APDTri", "Ca_Peak", "Ca_Diastole", "CaD90", "CaD50","Catri", "Qnet", "Qinward");

      // initialize some variables
      cipa_result.init( p_cell->STATES[V], p_cell->STATES[cai] );
      temp_result.init( p_cell->STATES[V], p_cell->STATES[cai] );
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
          fprintf( fp_feature, "CVODE_ERROR at pace %d\n", pace);
          break;
        }

        // calculate qnet and area under curve
        qnet += (p_cell->ALGEBRAIC[INaL]+p_cell->ALGEBRAIC[ICaL]+p_cell->ALGEBRAIC[Ito]+p_cell->ALGEBRAIC[IKr]+p_cell->ALGEBRAIC[IKs]+p_cell->ALGEBRAIC[IK1])*p_param->dt;
        inal_auc += p_cell->ALGEBRAIC[INaL]*p_param->dt;
        ical_auc += p_cell->ALGEBRAIC[ICaL]*p_param->dt;

        // save output for each pace
        if( iout % print_freq == 0 )  temp_result.vm_data.insert( std::pair<double, double> (tcurr, p_cell->STATES[V]) );
        if( iout % print_freq_ca == 0 )  temp_result.cai_data.insert( std::pair<double, double> (tcurr, p_cell->STATES[cai]) );

        // executed when entering new pace
        if( iout > 0 && iout % (int)( p_cell->CONSTANTS[BCL] / p_param->dt ) == 0 ) {
          // executed at the end of AP pacing
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

          if( is_eligible_AP && pace >= p_param->num_pace1-249) {
            // replace result with steeper repolarization AP or first pace from the last 250 paces
            if( temp_result.dvmdt_repol > cipa_result.dvmdt_repol ) {
              steepest_pace = pace;
              cipa_result = temp_result;
            }
          }

          // print result for each pace
          fprintf( fp_feature, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf ",
                 pace,
                 temp_result.dvmdt_repol,
                 temp_result.dvmdt_max,
                 temp_result.vm_peak,
                 temp_result.vm_dia,
                 temp_result.apd90,
                 temp_result.apd50,
                 temp_result.apd90-temp_result.apd50,
                 temp_result.ca_peak,
                 temp_result.ca_dia,
                 temp_result.cad90,
                 temp_result.cad50,
                 temp_result.cad90-temp_result.cad50);
          if((int)ceil(concs[idx]) == 0){
            INaL_auc_control = temp_result.inal_auc;
            ICaL_auc_control = temp_result.ical_auc;
            fprintf( fp_feature, "%lf %lf\n", pace, qnet/1000.0, 0.0 );
          }
          else{
            INaL_auc_drug = temp_result.inal_auc;
            ICaL_auc_drug = temp_result.ical_auc;
            fprintf( fp_feature, "%lf %lf\n", pace, qnet/1000.0, ( (INaL_auc_drug/INaL_auc_control) + (ICaL_auc_drug/ICaL_auc_control) ) * 0.5 );
          }


          pace++;

          // initialize cipa_result with the first temp_result
          // (only executed once, to avoid error when STATES[V] < 0)
          // (reported by KITOX at Monday 16th March 2021)
          if( pace == p_param->num_pace1-249 ) {
            cipa_result = temp_result;
          }

          qnet = inal_auc = ical_auc = 0.;

          // executed at the beginning of AP pacing
          //if(mympi::rank == 0)printf("Time new cycle: %lf with pace: %d\n", tcurr, pace);
          temp_result.init( p_cell->STATES[V], p_cell->STATES[cai] );
          t_ca_peak = tcurr;
          t_depol = (p_cell->CONSTANTS[BCL]*pace)+p_cell->CONSTANTS[stim_start];
          is_eligible_AP = false;

        }

        // get maximum dvmdt
        if( temp_result.dvmdt_max < p_cell->RATES[V] ) temp_result.dvmdt_max = p_cell->RATES[V];

        // get the peak Vm 30 secs after depolarization (when Na channel just closed after bursting)
        if( iout % (int)(((p_cell->CONSTANTS[BCL]*pace)+(p_cell->CONSTANTS[stim_start]+30.)) / p_param->dt) == 0) {
          temp_result.vm_peak = p_cell->STATES[V];
          // only finding the vm_repol30, vm_repol_50 and vm_repol90 if the peak Vm more than 0 mV
          if( temp_result.vm_peak > 0. ){
            vm_repol30 = temp_result.vm_peak - (0.3 * (temp_result.vm_peak - temp_result.vm_valley));
            vm_repol50 = temp_result.vm_peak - (0.5 * (temp_result.vm_peak - temp_result.vm_valley));
            vm_repol90 = temp_result.vm_peak - (0.9 * (temp_result.vm_peak - temp_result.vm_valley));
            is_eligible_AP = true;
          }
        }
        // these operations will be executed if it's eligible AP and executed at the beginning of repolarization
        else if( is_eligible_AP && tcurr > (p_cell->CONSTANTS[BCL]*pace)+(p_cell->CONSTANTS[stim_start]+30.) ){
          // get the steepest dvmdt_repol around vm_repol30 and vm_repol90 if AP shape is eligible
          if( vm_repol30 >= p_cell->STATES[V] && p_cell->STATES[V] >= vm_repol90 && temp_result.dvmdt_repol < p_cell->RATES[V] ){
              temp_result.dvmdt_repol = p_cell->RATES[V];
          }
          // get the APD50
          if( vm_repol50 > p_cell->STATES[V] && p_cell->STATES[V] > vm_repol50-2 ) 
            t_repol50 = tcurr;
            temp_result.apd50 = t_repol50 - t_depol;
          // get the APD90
          if( vm_repol90 > p_cell->STATES[V] && p_cell->STATES[V] > vm_repol90-2 ) 
            t_repol90 = tcurr;
            temp_result.apd90 = t_repol90 - t_depol;
          // get the peak calcium and amplitude 50% and 90%
          if( temp_result.ca_peak < p_cell->STATES[cai] ){  
            temp_result.ca_peak = p_cell->STATES[cai];
            ca_amp50 = temp_result.ca_peak - (0.5 * (temp_result.ca_peak - temp_result.ca_valley));
            ca_amp90 = temp_result.ca_peak - (0.9 * (temp_result.ca_peak - temp_result.ca_valley));
            t_ca_peak = tcurr;
          } 
        }
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
        }

        for(jdx = 0; jdx < p_cell->states_size; jdx++ ) fprintf( fp_output, "%lf\n", p_cell->STATES[jdx] );
      }

      if( p_param->is_print_graph == 1 ){
        fclose( fp_vm );
        fclose( fp_ca );
      }
      fclose( fp_feature );
      fclose( fp_output );

      N_VDestroy(states_vec);

    }  //end concs loop
  } // end sample loop
  double end = MPI_Wtime();

  CVodeFree(&cvode_mem);

  // wait for all nodes finished their job
  MPI_Barrier(MPI_COMM_WORLD);

  // master rank will make a performance log and zip the result file
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
    snprintf( buffer, sizeof(buffer), "result_%s_%s.zip", p_param->drug_name, buff_time );
    create_zip(buffer, "result");

    // delete the content inside result folder
    remove_dir_content("result");

    // copy the zip into the result folder
    char buff_new[150] = {'\0'};
    snprintf( buffer, sizeof(buffer), "result_%s_%s.zip", p_param->drug_name, buff_time );
    snprintf( buff_new, sizeof(buff_new), "result/result_%s_%s.zip", p_param->drug_name, buff_time );
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
    fprintf( fp_perf_log, "Iteration BCL Timestep CPU Sample Pacing Time(Minutes)\n" );
    fprintf( fp_perf_log, "%lld %.2lf %.2lf %d %d %d %lf\n", (int)(hill_data.size()*concs.size()*p_param->num_pace1*p_cell->CONSTANTS[BCL]/p_param->dt), p_cell->CONSTANTS[BCL], p_param->dt, mympi::size, hill_data.size(), p_param->num_pace1, (end-begin)/60. );
    fclose( fp_perf_log );

    printf("Computation Time: %lf minutes\n", (end-begin)/60.);

  }

  // clean memory
  CVodeFree(&cvode_mem);
  delete p_cell;
  free(p_param);  

}
