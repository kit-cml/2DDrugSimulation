#include "cipa_drug_explicit_bench.hpp"
#include "commons.hpp"
#include "globals.hpp"
#include "../cellmodels/Ohara_Rudy_2011.hpp"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <map>
#include <mpi.h>
#include <string>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

void cipa_drug_explicit_bench( int argc, char **argv, param_t *p_param, bool is_dutta )
{
  char buffer[255];
  char *token;

  /* Extra variables */
  int idx, jdx;
  FILE *fp_vm;

  /* ODE-related variables */
  double abs_err;
  double rel_err;
  double steps;
  double tnext;
  double tcurr;
  int iout,imax,iflag;
  int pace;

  /* CiPA variables */
  FILE *fp_hill;
  double hill_data[2000][14];
  int hill_size;
  double concs[5];
  int concs_size;

  // POSIX variables
  struct stat st = {0};

  // CellML object pointer
  patch_clamp *p_cell;

  // setup concentration array
  token = strtok( p_param->concs, "," );
  idx = concs_size = 0;
  while( token != NULL ){
    concs[idx++] = strtod(token, NULL);
    token = strtok( NULL, "," );
    concs_size++;
  }
  if(mympi::rank == 0){
    for( idx = 0; idx < concs_size; idx++ ){
      snprintf( buffer, sizeof(buffer), "result/%.2lf", concs[idx] );
      if(stat("result/", &st) == 0) mkdir(buffer, 0775);
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  

  // setup hill data array
  fp_hill = fopen(p_param->hill_file, "r");
  fgets(buffer, sizeof(buffer), fp_hill); // skip header
  hill_size = 0;
  while( fgets(buffer, sizeof(buffer), fp_hill) != NULL )
  { // begin data loop
    token = strtok( buffer, "," );
    idx = 0;
    while( token != NULL ){
      hill_data[hill_size][idx++] = strtod(token, NULL);
      token = strtok( NULL, "," );
    }
    hill_size++;

  } // end data loop
  fclose(fp_hill);

/*
  if(mympi::rank == 0){
    for( idx = 0; idx < hill_size; idx++ ){
      for( jdx = 0; jdx < 14; jdx++ ){
        printf("%lf ", hill_data[idx][jdx]);
      }
      printf("\n");
    }
  }
*/ 
  
 

  p_cell = new Ohara_Rudy_2011(  );
  p_cell-> initConsts(p_param->celltype, concs[0], hill_data[0], p_param->is_dutta);
  p_cell->CONSTANTS[BCL] = p_param->bcl_init;
  p_cell->CONSTANTS[amp] *= p_param->stim_amp;
  snprintf(buffer, sizeof(buffer), "vmcheck_explicit_rank%d.plt", mympi::rank);
  fp_vm = fopen(buffer, "w");
  fprintf( fp_vm, "%s %s \n", "TIME", "Vm" );

  pace = p_param->num_pace1;
  steps = p_param->dt;
  iout = 0;
  imax = ( pace *  p_cell->CONSTANTS[BCL] ) / steps;
  iflag = 1;
  tcurr = 0.0;
  tnext = steps;
  abs_err = 10E-5;
  rel_err = 10E-7;
  double work[100+21*p_cell->states_size];
  int iwork[5];

  double begin = MPI_Wtime();
  while(iout < imax){
    p_cell->computeRates( tcurr, p_cell->CONSTANTS, p_cell->RATES, p_cell->STATES, p_cell->ALGEBRAIC );
    p_cell->solveAnalytical( steps );

    fprintf( fp_vm, "%lf %lf \n", tcurr, p_cell->STATES[0] );
    iout++;
    tcurr += steps;
    
  }
  double end = MPI_Wtime();

  // wait for all nodes finished their job
  MPI_Barrier(MPI_COMM_WORLD);

  if(mympi::rank == 0) printf("Time elapsed for simulation is: %lf min\n", (end - begin)/60. );

  delete p_cell;
  fclose(fp_vm);

}
