#include "population_drug_bench.hpp"
#include "anm_drug_bench.hpp"
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

void population_drug_bench( int argc, char **argv, param_t *p_param, bool is_dutta )
{
  // CellML object pointer
  patch_clamp *p_cell;

  // Struct map to store anm result based on sample_id
  // first member is the concentration value
  // second member is the anm result data structure
  multimap< double, anm_result_t > sim_result_concs;
  // Struct map to store anm result based on sample_id
  // first member is the sample id
  // second member is the map of anm result data based on concentration
  multimap< int, multimap< double, anm_result_t > > sim_results_samples;


  
  // Supporting variables
  char buffer[150] = {'\0'};
  char *pch;
  double states_init[50];
  FILE *fp_hill;
  FILE *fp_states;
  int idx, jdx, group_id, sample_id, sim_counter, iout, imax;
  int print_freq;
  vector<double> temp_data;

  double hill_data_arr[14];
  FILE *fp_result;
  vector< vector<double> > hill_data;
  vector<double> concs;

  // POSIX variables
  struct stat st = {0};

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

  double begin = MPI_Wtime();

  int sample_group = ( hill_data.size()%mympi::size != 0 ) ? (hill_data.size()/mympi::size)+1 : (hill_data.size()/mympi::size);
  for( group_id=0, sim_counter=mympi::rank; 
       group_id < sample_group && sim_counter < hill_data.size();
       group_id++, sim_counter += mympi::size)
  { // begin sample loop
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
      p_cell->CONSTANTS[GNa] *=  p_param->gna_scale;


      
    }  //end concs loop
  } // end sample loop
  double end = MPI_Wtime();

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
    printf("Computation Time: %lf minutes\n", (end-begin)/60.);

    fclose( fp_perf_log );
  }

  // clean memory
  delete p_cell;
  free(p_param);
}
