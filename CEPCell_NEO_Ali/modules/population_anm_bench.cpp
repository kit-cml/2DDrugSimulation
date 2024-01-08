#include "population_anm_bench.hpp"
#include "anm_bench.hpp"
#include "commons.hpp"
#include "globals.hpp"
#include "../cellmodels/tentusscher_noble_noble_panfilov_2004_b.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>

#include <mpi.h>


void population_anm_bench(int argc, char **argv, param_t *p_param)
{
  // Struct to store result
  anm_result_t sim_result;

  // Cell object variables
  patch_clamp *p_cell;

  char scenario_name[255] = {'\0'};
  char buffer[255] = {'\0'};
  double scalars[10];

  FILE *fp_result;

  int idx,jdx;
  int scalar_base4;

  p_cell = new tentusscher_noble_noble_panfilov_2004_b();

  snprintf(buffer, sizeof(buffer), "result/mapping_proc%d.plt", mympi::rank);
  fp_result = fopen(buffer, "a");
  fprintf( fp_result, "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n",
           "g_Ks", "g_Kr", "g_K1", "g_Na", "g_bna",
           "g_CaL", "g_bca", "g_to", "g_pCa", "g_pK",
           "AOCL", "mean_APD", "AOCL/mean_APD", "mean_APD/AOCL", "ANM_AOCL", 
           "ANM_Alternant", "CL_Alternant", "Mean_APD_Alternant" );
  fclose(fp_result);


  double begin = MPI_Wtime();
  for( idx = mympi::rank; idx < 65536; idx += mympi::size ){  // 4^8 = 65536
  //for( idx = mympi::rank; idx < 6561; idx += mympi::size ){  // 3^8 = 6561
    // convert idx, which is base-10 number, into base-4 number.
    // Then, use the digits from scalar_base4 to populate the scalars.
    scalar_base4 = convert_base_n(idx, 4);
    //scalar_base4 = convert_base_n(idx, 3);
    for( jdx = 9; jdx >= 0; jdx-- ){
      scalars[jdx] = (double)(25 * (4-(scalar_base4 % 10)))/100.;
      //scalars[jdx] = (double)(50 * (3-(scalar_base4 % 10)))/100.;
      scalar_base4 /= 10;
    }
    scalars[1] = p_param->gkr_scale;
    scalars[0] = p_param->gks_scale;

    snprintf(scenario_name, sizeof(scenario_name),
             "g_Ks%.2fg_Kr%.2fg_K1%.2fg_Na%.2fg_bna%.2fg_CaL%.2fg_bca%.2fg_to%.2fg_pCa%.2fg_pK%.2f",
              scalars[0],scalars[1],scalars[2],scalars[3],scalars[4],
              scalars[5],scalars[6],scalars[7],scalars[8],scalars[9]);

    p_cell->initConsts();
    p_cell->CONSTANTS[g_Ks] *= scalars[0];
    p_cell->CONSTANTS[g_Kr] *= scalars[1];
    p_cell->CONSTANTS[g_K1] *= scalars[2];
    p_cell->CONSTANTS[g_Na] *= scalars[3];
    p_cell->CONSTANTS[g_bna] *= scalars[4];
    p_cell->CONSTANTS[g_CaL] *= scalars[5];
    p_cell->CONSTANTS[g_bca] *= scalars[6];
    p_cell->CONSTANTS[g_to] *= scalars[7];
    p_cell->CONSTANTS[g_pCa] *= scalars[8];
    p_cell->CONSTANTS[g_pK] *= scalars[9];

    sim_result = anm_bench(argc,argv,p_param,p_cell,scenario_name);

    fp_result = fopen(buffer, "a");
    fprintf( fp_result, "%.2lf %.2lf %.2lf %.2lf %.2lf %.2lf %.2lf %.2lf %.2lf %.2lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
             scalars[0],scalars[1],scalars[2],scalars[3],scalars[4],
             scalars[5],scalars[6],scalars[7],scalars[8],scalars[9],
             sim_result.aocl, sim_result.mean_apd, sim_result.aocl/sim_result.mean_apd, sim_result.mean_apd/sim_result.aocl, sim_result.anm_aocl,
             sim_result.anm_alternant, sim_result.cl_alternant, sim_result.mean_apd_alternant );
    fclose(fp_result);    

  }
  double end = MPI_Wtime();

  // wait for all nodes finished their job
  MPI_Barrier(MPI_COMM_WORLD);

  if(mympi::rank == 0) printf("Computation Time: %lf minutes\n", (end-begin)/60.);
#ifdef TESTESTES
  {
    char buff_time[15] = {'\0'};
    FILE *fp_perf_log;
    struct tm* tm_info;
    time_t timer;

    timer = time(NULL);
    tm_info = localtime(&timer);
    strftime( buff_time, 15, "%Y%m%d%H%M%S", tm_info );
    snprintf( buffer, sizeof(buffer), "result/performance_%s.log", buff_time );

    // zip result folder
    snprintf( buffer, sizeof(buffer), "result_%.0lf%.0lf_%s.zip", p_param->gks_scale*100, p_param->gkr_scale*100, buff_time );
    create_zip(buffer, "result");

    // delete the content inside result folder
    remove_dir_content("result");

    // copy the zip into the result folder
    char buff_new[150] = {'\0'};
    snprintf( buffer, sizeof(buffer), "result_%.0lf%.0lf_%s.zip", p_param->gks_scale*100, p_param->gkr_scale*100, buff_time );
    snprintf( buff_new, sizeof(buff_new), "result/result_%.0lf%.0lf_%s.zip", p_param->gks_scale*100, p_param->gkr_scale*100, buff_time );
    if( rename( buffer, buff_new ) == 0 ){
      printf("Zip file %s moved succesfully to %s!!\n", buffer, buff_new);
    }
    else{
      printf("Error moving file %s to %s!!\n", buffer, buff_new);
      return;
    }

  }
#endif
  delete p_cell;
}
