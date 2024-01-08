#include "param.hpp"

#include <cstdio>
#include "glob_funct.hpp"

void param_t::init()
{
  simulation_mode = 0;
  is_dutta = true;
  is_print_graph = true;
  is_using_output = false;
  bcl = 2000.;
  pace_max = 1000;
  celltype = 0;
  last_drug_check_pace = 250;
  dt = 0.1;
  dt_write = 2.0;
  inet_vm_threshold = -88.0;
  snprintf(hill_file, sizeof(hill_file), "%s", "./drugs/bepridil/IC50_samples10.csv");
  snprintf(herg_file, sizeof(herg_file), "%s", "./herg_pars/bepridil/boot_pars10.csv");
  snprintf(drug_name, sizeof(drug_name), "%s", "bepridil");
  snprintf(concs, sizeof(concs), "%s", "99.0");
}

void param_t::show_val()
{
  mpi_printf( 0, "%s -- %s\n", "Simulation mode", simulation_mode ? "full-pace" : "sample-based" );
  mpi_printf( 0, "%s -- %s\n", "Hill File", hill_file );
  mpi_printf( 0, "%s -- %s\n", "hERG File", herg_file );
  mpi_printf( 0, "%s -- %hu\n", "Celltype", celltype);
  mpi_printf( 0, "%s -- %s\n", "Is_Dutta", is_dutta ? "true" : "false" );
  mpi_printf( 0, "%s -- %s\n", "Is_Print_Graph", is_print_graph ? "true" : "false" );
  mpi_printf( 0, "%s -- %s\n", "Is_Using_Output", is_using_output ? "true" : "false" );
  mpi_printf( 0, "%s -- %lf\n", "Basic_Cycle_Length", bcl);
  mpi_printf( 0, "%s -- %hu\n", "Number_of_Pacing", pace_max);
  mpi_printf( 0, "%s -- %hu\n", "Last_Drug_Check_Pace", last_drug_check_pace);
  mpi_printf( 0, "%s -- %lf\n", "Time_Step", dt);
  mpi_printf( 0, "%s -- %lf\n", "Inet_Vm_Threshold", inet_vm_threshold);
  mpi_printf( 0, "%s -- %lf\n", "Writing_Step", dt_write);
  mpi_printf( 0, "%s -- %s\n", "Drug_Name", drug_name);
  mpi_printf( 0, "%s -- %s\n", "Concentrations", concs);
}
