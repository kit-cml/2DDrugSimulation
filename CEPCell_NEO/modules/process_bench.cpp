#include "process_bench.hpp"

#include "globals.hpp"

#include "anm_bench.hpp"
#include "anm_bench_custom.hpp"
#include "ap_bench_orudy.hpp"
#include "ap_bench_tn2004.hpp"
#include "ap_bench_tomek.hpp"
#include "apdr_bench.hpp"
#include "bradycardia_bench.hpp"
#include "cipa_drug_bench.hpp"
#include "cipa_drug_tomek_bench.hpp"
#include "cipa_drug_full_bench.hpp"
#include "cipa_drug_full_tomek_bench.hpp"
#include "cipa_drug_last_bench.hpp"
#include "cipa_drug_tn2004_bench.hpp"
#include "patch_bench.hpp"
#include "population_anm_bench.hpp"
#include "population_drug_bench.hpp"
#include "population_qnet_bench.hpp"
#include "population_qnet_all_bench.hpp"
#include "qnet_bench.hpp"
#include "qnet_all_bench.hpp"
#include "hrv_bench.hpp"
#include "hrv_bench_combined.hpp"
#include "hrv_bench_optimal.hpp"

#include <cstdio>
#include <cstdlib>
#include <mpi.h>

#if defined( __linux__ )
  #include <sys/types.h>
  #include <sys/stat.h>
  #include <unistd.h>
#elif defined( _WIN32 )
  #include <direct.h>
#endif

void process(int argc, char* argv[])
{
  setvbuf( stdout, NULL, _IONBF, 0 );

  MPI_Init( &argc, &argv );
  MPI_Comm_size( MPI_COMM_WORLD, &mympi::size );
  MPI_Comm_rank( MPI_COMM_WORLD, &mympi::rank );

  param_t *p_param;

  p_param = (param_t*)malloc( sizeof( param_t ) );

  set_default_values(p_param);

  edison_assign_params_single(argc, argv, p_param);

  struct stat st = { 0 };
  // create result folder if not existed
  if( mympi::rank == 0 ){
    if (stat("result/", &st) == -1) mkdir("result", 0775);
  }

  /* check the modules folder to look the implementation */
  if( p_param->simulation_mode == 0 ){
    if(mympi::rank == 0) printf("Using AP Simulation with TN2004 model!\n");
    ap_bench_tn2004( argc, argv, p_param );
  }
  else if( p_param->simulation_mode == 10 ){
    if(mympi::rank == 0) printf("Using AP Simulation with ORudy model!\n");
    ap_bench_orudy( argc, argv, p_param );
  }
  else if( p_param->simulation_mode == 100 ){
    if(mympi::rank == 0) printf("Using AP Simulation with ToR-ORd model!\n");
    ap_bench_tomek( argc, argv, p_param );
  }
  else if( p_param->simulation_mode == 1 ){
    if(mympi::rank == 0) printf("Using APDR Simulation!\n");
    apdr_bench( argc, argv, p_param );
  }
  else if( p_param->simulation_mode == 2 ){
    if(mympi::rank == 0) printf("Using Drug_Hill Simulation!\n");
    cipa_drug_bench( argc, argv, p_param, false );
  }
  else if( p_param->simulation_mode == 3 ){
    if(mympi::rank == 0) {
      printf("Using Drug_Hill Simulation with Dutta Model!\n");
    }
    cipa_drug_bench( argc, argv, p_param, p_param->is_dutta );
  }
  else if( p_param->simulation_mode == 33 ){
    if(mympi::rank == 0) {
      printf("Using Drug_Hill Simulation with Dutta Model (Last)!\n");
    }
    cipa_drug_last_bench( argc, argv, p_param, p_param->is_dutta );
  }
  else if( p_param->simulation_mode == 333 ){
    if(mympi::rank == 0){
      printf("Using Drug_Hill Simulation with Tn2004 model!\n");
    }
    cipa_drug_tn2004_bench( argc, argv, p_param, p_param->is_dutta );
  }
  else if( p_param->simulation_mode == 3333 ){
    if(mympi::rank == 0){
      printf("Using Drug_Hill Simulation for ALL PACE!\n");
    }
    cipa_drug_full_bench( argc, argv, p_param, p_param->is_dutta );
  }
  else if( p_param->simulation_mode == 33333 ){
    if(mympi::rank == 0){
      printf("Using Drug_Hill Simulation using Tomek model!\n");
    }
    cipa_drug_tomek_bench( argc, argv, p_param );
  }
  else if( p_param->simulation_mode == 333333 ){
    if(mympi::rank == 0){
      printf("Using Drug_Hill Simulation using Tomek model for ALL PACE!\n");
    }
    cipa_drug_full_tomek_bench( argc, argv, p_param );
  }
  else if( p_param->simulation_mode == 4 ){
    if(mympi::rank == 0) printf("Using Drug_Herg Simulation!\n");
    //cipa_drug_herg_bench( argc, argv, p_param );
  }
  else if( p_param->simulation_mode == 5 ){
    if(mympi::rank == 0) printf("Using Patch Clamp Simulation!\n");
    patch_bench( argc, argv, p_param );
  }
  else if( p_param->simulation_mode == 6 ){
    if(mympi::rank == 0) printf("Using Bradycardia Simulation!\n");
    bradycardia_bench( argc, argv, p_param );
  }
  else if( p_param->simulation_mode == 7 ){
    glob_var::is_hrv = true;
    if(mympi::rank == 0) printf("Using HRV Simulation!\n");
    hrv_bench( argc, argv, p_param );
  }
  else if( p_param->simulation_mode == 8 ){
    glob_var::is_hrv = true;
    if(mympi::rank == 0) printf("Using HRV Optimal Simulation!\n");
    hrv_bench_optimal( argc, argv, p_param );
  }
  else if( p_param->simulation_mode == 9 ){
    anm_result_t sim_result;
    if(mympi::rank == 0) printf("Using ANM Simulation!\n");
    sim_result = anm_bench( argc, argv, p_param, NULL, NULL );
    printf("The AOCL is %lf with mean APD %lf\n", sim_result.aocl, sim_result.mean_apd);    
  }
  else if( p_param->simulation_mode == 99 ){
    anm_result_t sim_result;
    if(mympi::rank == 0) printf("Using ANM Simulation Custom Scale!\n");
    sim_result = anm_bench_custom( argc, argv, p_param, NULL, NULL );
    printf("The AOCL is %lf with mean APD %lf\n", sim_result.aocl, sim_result.mean_apd);    
  }
  else if( p_param->simulation_mode ==11 ){
    if(mympi::rank == 0) printf("Using ANM Population Simulation!\n");
    population_anm_bench( argc, argv, p_param  );
  }
  else if( p_param->simulation_mode == 12 ){
    cipa_t sim_result;
    if(mympi::rank == 0) printf("Using Qnet Simulation!\n");
    sim_result = qnet_bench( argc, argv, p_param, NULL, NULL );
    printf("The qnet is %lf with maximum dvm/dt repol %lf\n", sim_result.qnet, sim_result.dvmdt_repol);
  }
  else if( p_param->simulation_mode == 13 ){
    if(mympi::rank == 0) printf("Using Qnet Population Simulation!\n");
    population_qnet_bench( argc, argv, p_param );
  }
  else if( p_param->simulation_mode == 14 ){
    cipa_t sim_result;
    if(mympi::rank == 0) printf("Using Qnet (All result) Simulation!\n");
    sim_result = qnet_all_bench( argc, argv, p_param, NULL, NULL );
    printf("The qnet is %lf with maximum dvm/dt repol %lf\n", sim_result.qnet, sim_result.dvmdt_repol);
    printf("DVM/DT_MAX %lf VM_PEAK %lf VM_DIA %lf APD50 %lf APD90 %lf APD_TRI %lf\n", sim_result.dvmdt_max, sim_result.vm_peak, sim_result.vm_dia, sim_result.apd50, sim_result.apd90, sim_result.apd90 - sim_result.apd50);
    printf("CA_PEAK %lf CA_DIA %lf CAD50 %lf CAD90 %lf CA_TRI %lf\n", sim_result.ca_peak, sim_result.ca_dia, sim_result.cad50, sim_result.cad90, sim_result.cad90 - sim_result.cad50);
  }
  else if( p_param->simulation_mode == 15 ){
    if(mympi::rank == 0) printf("Using Qnet (All result) Population Simulation!\n");
    population_qnet_all_bench( argc, argv, p_param );
  }
  else if( p_param->simulation_mode == 777 ){
    glob_var::is_hrv = true;
    if(mympi::rank == 0) printf("Using HRV Multi-Sample Multi-HRV Simulation!\n");
    hrv_bench_combined( argc, argv, p_param );
  }



  else{
    if(mympi::rank == 0) printf("Selection unknown!!\n");
  }

  MPI_Finalize();


}
