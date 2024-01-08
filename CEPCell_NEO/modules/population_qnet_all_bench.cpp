#include "population_qnet_all_bench.hpp"
#include "qnet_all_bench.hpp"
#include "commons.hpp"
#include "globals.hpp"
#include "../cellmodels/Ohara_Rudy_2011.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>

#include <mpi.h>


void population_qnet_all_bench(int argc, char **argv, param_t *p_param)
{
  // Struct to store result
  cipa_t sim_result;

  // Cell object variables
  patch_clamp *p_cell;

  char scenario_name[255] = {'\0'};
  char buffer[255] = {'\0'};
  double scalars[10];

  FILE *fp_result;

  int idx,jdx;
  int scalar_base4;

  p_cell = new Ohara_Rudy_2011();

  snprintf(buffer, sizeof(buffer), "result/mapping_proc%d.plt", mympi::rank);
  fp_result = fopen(buffer, "a");
  fprintf( fp_result, "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n",
           "GKs", "GKr", "GK1", "GNaL", "PCa", "Gto", "qnet", "dvm/dt_repol",
           "dvm/dt_max", "vm_peak", "vm_dia", "apd50", "apd90", "apd_tri",
           "ca_peak", "ca_dia", "cad50", "cad90", "cad_tri"  );
  fclose(fp_result);


  double begin = MPI_Wtime();

  for( idx = mympi::rank; idx < 4096; idx += mympi::size ){  // 4^6 = 4096
    // convert idx, which is base-10 number, into base-4 number.
    // Then, use the digits from scalar_base4 to populate the scalars.
    scalar_base4 = convert_base_n(idx, 4);
    for( jdx = 5; jdx >= 0; jdx-- ){
      scalars[jdx] = (double)( 25 * (4-(scalar_base4 % 10) ))/100.;
      scalar_base4 /= 10;
    }

    snprintf(scenario_name, sizeof(scenario_name),
             "GKs%.2fGKr%.2fGK1%.2fGNaL%.2fPCa%.2fGto%.2f",
              scalars[0],scalars[1],scalars[2],scalars[3],scalars[4],scalars[5]);

    p_cell->initConsts((double)p_param->celltype, !!p_param->is_dutta);
    p_cell->CONSTANTS[GKs] *= scalars[0];
    p_cell->CONSTANTS[GKr] *= scalars[1];
    p_cell->CONSTANTS[GK1] *= scalars[2];
    p_cell->CONSTANTS[GNaL] *= scalars[3];
    p_cell->CONSTANTS[PCa] *= scalars[4];
    p_cell->CONSTANTS[Gto] *= scalars[5];
    p_cell->CONSTANTS[BCL] = p_param->bcl_init;

    sim_result = qnet_all_bench(argc,argv,p_param,p_cell,scenario_name);

    fp_result = fopen(buffer, "a");
    fprintf( fp_result, "%.2lf %.2lf %.2lf %.2lf %.2lf %.2lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
             scalars[0],scalars[1],scalars[2],scalars[3],scalars[4],scalars[5],sim_result.qnet, sim_result.dvmdt_repol,
             sim_result.dvmdt_max, sim_result.vm_peak, sim_result.vm_dia, sim_result.apd50, sim_result.apd90, sim_result.apd90 - sim_result.apd50,
             sim_result.ca_peak, sim_result.ca_dia, sim_result.cad50, sim_result.cad90, sim_result.cad90 - sim_result.cad50 );
    fclose(fp_result);    

  }

  double end = MPI_Wtime();

  // wait for all nodes finished their job
  MPI_Barrier(MPI_COMM_WORLD);

  if(mympi::rank == 0) printf("Computation Time: %lf minutes\n", (end-begin)/60.);

  delete p_cell;

}
