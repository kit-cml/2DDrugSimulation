#include "process_bench.hpp"

#include "ap_bench.hpp"
#include "apdr_bench.hpp"
#include "cipa_drug_bench.hpp"
#include "cipa_drug_herg_bench.hpp"
#include "patch_bench.hpp"

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

  MPI_Init( &argc, &argv );
  MPI_Comm_size( MPI_COMM_WORLD, &mympi::size );
  MPI_Comm_rank( MPI_COMM_WORLD, &mympi::rank );

  param_t *p_param;

  p_param = (param_t*)malloc( sizeof( param_t ) );

  set_default_values(p_param);

  edison_assign_params_single(argc, argv, p_param);

#if defined( __linux__ )
  struct stat st = { 0 };
  // create result folder if not existed
  if( mympi::rank == 0 ){
    if (stat("result/", &st) == -1) mkdir("result", 0775);
  }
#elif defined( _WIN32 )
  _mkdir("result/");
#endif

  /* check the modules folder to look the implementation */
  if( p_param->simulation_mode == 0 ){
    if(mympi::rank == 0) printf("Using AP Simulation!\n");
//    ap_bench( argc, argv, p_param );
  }
  else if( p_param->simulation_mode == 1 ){
    if(mympi::rank == 0) printf("Using APDR Simulation!\n");
//    apdr_bench( argc, argv, p_param );
  }
  else if( p_param->simulation_mode == 2 ){
    if(mympi::rank == 0) printf("Using Drug_Hill Simulation!\n");
    cipa_drug_bench( argc, argv, p_param, false );
  }
  else if( p_param->simulation_mode == 3 ){
    if(mympi::rank == 0) printf("Using Drug_Hill Simulation with Dutta Model!\n");
    cipa_drug_bench( argc, argv, p_param, true );
  }
  else if( p_param->simulation_mode == 4 ){
    if(mympi::rank == 0) printf("Using Drug_Herg Simulation!\n");
    cipa_drug_herg_bench( argc, argv, p_param );
  }
  else if( p_param->simulation_mode == 5 ){
    if(mympi::rank == 0) printf("Using Patch Clamp Simulation!\n");
    patch_bench( argc, argv, p_param );
  }
  else{
    if(mympi::rank == 0) printf("Selection unknown!!\n");
  }

  MPI_Finalize();


}
