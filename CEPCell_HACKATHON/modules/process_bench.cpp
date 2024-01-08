#include "process_bench.hpp"

#include "globals.hpp"

#include "cipa_drug_bench.hpp"
#include "cipa_drug_shago_bench.hpp"
#include "cipa_drug_explicit_bench.hpp"

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

  if( p_param->simulation_mode == 3 ){
    if(mympi::rank == 0) {
      printf("Using Drug_Hill Simulation with Dutta Model CVODE!\n");
    }
    cipa_drug_bench( argc, argv, p_param, p_param->is_dutta );
  }
  else if( p_param->simulation_mode == 255 ){
    if(mympi::rank == 0) {
      printf("Using Drug_Hill Simulation with Dutta Model Shampine-Gordon!\n");
    }
    cipa_drug_shago_bench( argc, argv, p_param, p_param->is_dutta );
  }
  else if( p_param->simulation_mode == 510 ){
    if(mympi::rank == 0) {
      printf("Using Drug_Hill Simulation with Dutta Model Explicit!\n");
    }
    cipa_drug_explicit_bench( argc, argv, p_param, p_param->is_dutta );
  }


  else{
    if(mympi::rank == 0) printf("Selection unknown!!\n");
  }

  MPI_Finalize();


}
