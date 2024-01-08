#include "cipa_drug_tn2004_bench.hpp"
#include "commons.hpp"
#include "globals.hpp"
#include "../cellmodels/tentusscher_noble_noble_panfilov_2004_b.hpp"

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

void cipa_drug_tn2004_bench( int argc, char **argv, param_t *p_param, bool is_dutta )
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
  double hill_data_arr[14];
  FILE *fp_vm = NULL;
  int pace, steepest_pace;
  vector< vector<double> > hill_data;
  vector<double> concs;

  // SUNDIALs variables
  bool cvode_firsttime;
  int cvode_retval;
  double tnext, tcurr = 0.0;
  void *cvode_mem;
  N_Vector states_vec;

  // POSIX variables
  struct stat st = {0};

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

  p_cell = new tentusscher_noble_noble_panfilov_2004_b();
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
      p_cell->initConsts(p_param->celltype, concs[idx], hill_data_arr);
      p_cell->CONSTANTS[stim_period] = p_param->bcl_init;
      p_cell->CONSTANTS[stim_amplitude] *= p_param->stim_amp;

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
      }

      // Write header for all of the result file
      if( group_id == 0 )
      {
        if( p_param->is_print_graph == 1 ){
          fprintf( fp_vm, "%s %s\n", "TIME", "Vm" );
        }
      }

      // initialize some variables
      pace = 0;
      iout = 0;
      imax = ( p_param->num_pace1 *  p_cell->CONSTANTS[stim_period] ) / p_param->dt;

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

        if( iout > 0 && iout % (int)( p_cell->CONSTANTS[stim_period] / p_param->dt ) == 0 ) pace++;

        if(pace >= p_param->num_pace1-1 && iout % print_freq == 0) fprintf(fp_vm, "%lf %lf\n", tcurr, p_cell->STATES[V]);


      } // end paces loop

      fclose( fp_vm );
      N_VDestroy(states_vec);

    }  //end concs loop
  } // end sample loop
  double end = MPI_Wtime();

  CVodeFree(&cvode_mem);
  delete p_cell;

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
    fprintf( fp_perf_log, "Iteration CPU Sample Pacing Time(Minutes)\n" );
    fprintf( fp_perf_log, "%lld %d %d %d %lf\n", (int)(hill_data.size()*concs.size()*p_param->num_pace1*p_cell->CONSTANTS[stim_period]/p_param->dt), mympi::size, hill_data.size(), p_param->num_pace1, (end-begin)/60. );
    printf("Computation Time: %lf minutes\n", (end-begin)/60.);

    fclose( fp_perf_log );
  }

}
