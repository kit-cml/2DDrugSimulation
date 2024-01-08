#include "commons.hpp"
#include "../cellmodels/patch_clamp.hpp"

#include <cstdio>
#include <cstdlib>
#include <cstring>

int mympi::rank, mympi::size;

int rhs_fn( realtype t, N_Vector y, N_Vector ydot, void *user_data )
{
  patch_clamp *data = (patch_clamp*)user_data;
  data->computeRates( t,
                      data->CONSTANTS,
                      N_VGetArrayPointer_Serial(ydot),
                      N_VGetArrayPointer_Serial(y),
                      data->ALGEBRAIC );
  return 0;
}

void set_default_values( param_t *p_param )
{
  p_param->celltype = 0;
  p_param->simulation_mode = 0;
  p_param->t_write_vtk = 0.0;
  p_param->dt = 0.01;
  p_param->dt_write = 1.0;
  p_param->stim_amp = 1.0;
  p_param->stim_dur = 0.5;
  p_param->num_pace1 = 3;
  p_param->bcl_decrement = 50.0;
  p_param->bcl_end = 2000000.0;
  p_param->bcl_init = 1000.0;
  p_param->mutation_scaling = 0.0;
  p_param->v_apd90 = -40.0;
  p_param->is_print_graph = 1;
  p_param->is_s1s2 = 0;
  p_param->is_using_output = 0;
  p_param->is_print_graph = 1;
  sprintf( p_param->variant, "%s", "ori" );
  sprintf( p_param->drug_name, "%s", "bepridil" );
  p_param->cmax = 33;
  sprintf( p_param->concs, "%s", "0.0" );

}


void edison_assign_params_single(int argc, char *args[], param_t *p_param)
{
  char buffer[100];
  char key[100];
  char value[100];
  char file_name[255];
  FILE *fp_inputdeck;

  for (int i = 1; i < argc; i += 2) {
    if (!strcmp(args[i], "-input_deck"))
      strcpy(file_name, args[i + 1]);
    else if (!strcmp(args[i], "-hill_file"))
      strcpy(p_param->hill_file, args[i + 1]);
    else if (!strcmp(args[i], "-herg_file"))
      strcpy(p_param->herg_file, args[i + 1]);
  }

  if(strlen(p_param->hill_file) != 0){
    if(mympi::rank == 0)  printf( "Hill File: %s\n", p_param->hill_file );
  }
  if(strlen(p_param->herg_file) != 0){
    if(mympi::rank == 0) printf( "hERG File: %s\n", p_param->herg_file );
  }

  fp_inputdeck = fopen( file_name, "r");
  if(fp_inputdeck == NULL){
    fprintf(stderr, "Cannot open file %s!!!\n", file_name);
    exit(0);
  }

  // read input_deck line by line
  // and store each line to the buffer
  while ( fgets( buffer, 100, fp_inputdeck ) != NULL ) {
    // parse the buffer and store it to key and value
    sscanf( buffer, "%s %*s %s", key, value );
    if (strcmp(key, "Celltype") == 0) {
      p_param->celltype = strtod( value, NULL );
      if(mympi::rank == 0) printf( "%s -- %lf\n", "Celltype", p_param->celltype);
    }
    else if (strcmp(key, "Variant") == 0) {
      strcpy( p_param->variant, value );
      if(mympi::rank == 0) printf( "%s -- %s\n", "Variant", p_param->variant);
    }
    else if (strcmp(key, "Simulation_Mode") == 0) {
      p_param->simulation_mode = strtol( value, NULL, 10 );
      if(mympi::rank == 0) printf( "%s -- %d\n", "Simulation_Mode", p_param->simulation_mode);
    }
    else if (strcmp(key, "Write_Time") == 0) {
      p_param->t_write_vtk = strtod( value, NULL );
      if(mympi::rank == 0) printf( "%s -- %lf\n", "Write_Time", p_param->t_write_vtk);
    }
    else if (strcmp(key, "Number_of_Pacing") == 0) {
      p_param->num_pace1 = strtol( value, NULL, 10 );
      if(mympi::rank == 0) printf( "%s -- %d\n", "Number_of_Pacing", p_param->num_pace1);
    }
    else if (strcmp(key, "Basic_Cycle_Length") == 0) {
      p_param->bcl_init = strtod( value, NULL );
      if(mympi::rank == 0) printf( "%s -- %lf\n", "Basic_Cycle_Length", p_param->bcl_init);
    }
    else if (strcmp(key, "Time_Step") == 0) {
      p_param->dt = strtod( value, NULL );
      if(mympi::rank == 0) printf( "%s -- %lf\n", "Time_Step", p_param->dt);
    }
    else if (strcmp(key, "Write_Step") == 0) {
      p_param->dt_write = strtod( value, NULL );
      if(mympi::rank == 0) printf( "%s -- %lf\n", "Writing_Step", p_param->dt_write);
    }
    else if (strcmp(key, "Stim_Amplitude") == 0) {
      p_param->stim_amp = strtod( value, NULL );
      if(mympi::rank == 0) printf( "%s -- %lf\n", "Stim_Amplitude", p_param->stim_amp);
    }
    else if (strcmp(key, "Stim_Duration") == 0) {
      p_param->stim_dur = strtod( value, NULL );
      if(mympi::rank == 0) printf( "%s -- %lf\n", "Stim_Duration", p_param->stim_dur);
    }
    else if (strcmp(key, "BCL_Step") == 0) {
      p_param->bcl_decrement = strtod( value, NULL );
      if(mympi::rank == 0) printf( "%s -- %lf\n", "BCL_Step", p_param->bcl_decrement);
    }
    else if (strcmp(key, "BCL_End") == 0) {
      p_param->bcl_end = strtod( value, NULL );
      if(mympi::rank == 0) printf( "%s -- %lf\n", "BCL_End", p_param->bcl_end);
    }
    else if (strcmp(key, "APD_90_Percent") == 0) {
      p_param->v_apd90 = strtod( value, NULL );
      if(mympi::rank == 0) printf( "%s -- %lf\n", "APD_90_Percent", p_param->v_apd90);
    }
    else if (strcmp(key, "APD_50_Percent") == 0) {
      p_param->v_apd50 = strtod( value, NULL );
      if(mympi::rank == 0) printf( "%s -- %lf\n", "APD_50_Percent", p_param->v_apd50);
    }
    else if (strcmp(key, "Scaled_Gate") == 0) {
      p_param->scaled_gate = strtol( value, NULL, 10 );
      if(mympi::rank == 0) printf( "%s -- %d\n", "Scaled_Gate", p_param->scaled_gate);
    }
    else if (strcmp(key, "Drug_Name") == 0) {
      strcpy( p_param->drug_name, value );
      if(mympi::rank == 0) printf( "%s -- %s\n", "Drug_Name", p_param->drug_name);
    }
    else if (strcmp(key, "Concentrations") == 0) {
      strcpy( p_param->concs, "0," );
      strcat( p_param->concs, value );
      if(mympi::rank == 0) printf( "%s -- %s\n", "Concentrations", p_param->concs);
    }
    else if (strcmp(key, "GKs_Mult") == 0) {
      p_param->gks_scale =  strtod(value, NULL) ;
      if(mympi::rank == 0) printf( "%s -- %lf\n", "GKs_Mult", p_param->gks_scale);
    }
    else if (strcmp(key, "GCaL_Mult") == 0) {
      p_param->gcal_scale =  strtod(value, NULL) ;
      if(mympi::rank == 0) printf( "%s -- %lf\n", "GCaL_Mult", p_param->gcal_scale);
    }
    else if (strcmp(key, "GNa_Mult") == 0) {
      p_param->gna_scale =  strtod(value, NULL) ;
      if(mympi::rank == 0) printf( "%s -- %lf\n", "GNa_Mult", p_param->gna_scale);
    }
    else if (strcmp(key, "Is_Print_Graph") == 0) {
      p_param->is_print_graph = strtol( value, NULL, 10 );
      if(mympi::rank == 0) printf( "%s -- %d\n", "Is_Print_Graph", p_param->is_print_graph );
    }
    else if (strcmp(key, "Is_Using_Output") == 0) {
      p_param->is_using_output = strtol( value, NULL, 10 );
      if(mympi::rank == 0) printf( "%s -- %d\n", "Is_Using_Output", p_param->is_using_output );
    }

  }
  fclose( fp_inputdeck );
}

