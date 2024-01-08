#include "patches/Ohara_Rudy_2011.hpp"

#include "param.hpp"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <cvode/cvode.h>
#include <cvode/cvode_dense.h>
#include <nvector/nvector_serial.h>

#define USE_CVODE

int rhs_fn( realtype t, N_Vector y, N_Vector ydot, void *user_data );
void edison_assign_params_single(int argc, char *args[], param_t *p_param);

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



int main( int argc, char **argv )
{
  // I/O variables
  char buff[255];
  char *pch;
  FILE *fp_curr;
  FILE *fp_curr10;
  FILE *fp_vm;
  FILE *fp_vm10;

  // Patch clamp variables
  patch_clamp *p_patch;
  double bcl,dt,pace;

  // SUNDIALs variables
  int retval;
  double tnext, tcurr;
  void *cvode_main;
  N_Vector states_vec;

  param_t *p_param;
  p_param = (param_t*)malloc( sizeof( param_t ) );
  edison_assign_params_single(argc,argv,p_param);

  p_patch = new Ohara_Rudy_2011();
  p_patch->initConsts();
  if((int)floor(p_param->bcl_init) == 0)bcl = 1000.0;
  else bcl = p_param->bcl_init;
  p_patch->CONSTANTS[stim_period] = bcl;
  if(p_param->dt > 0)dt = p_param->dt;
  else dt = 0.01;
  if(p_param->num_pace1 == 0) pace = 3;
  else pace = p_param->num_pace1;
  tcurr = 0.0;
  printf( "GK1 before: %lf\n", p_patch->CONSTANTS[GK1] );
  if(p_param->gk1_mult > 0) p_patch->CONSTANTS[GK1] *= p_param->gk1_mult;
  else p_patch->CONSTANTS[GK1] *= 1;
  printf( "GK1 after: %lf\n", p_patch->CONSTANTS[GK1] );

#if defined USE_CVODE
  // Create CVODE solver
  cvode_main = CVodeCreate( CV_BDF,CV_NEWTON );
  // Give p_patch as CVode User Data
  CVodeSetUserData( cvode_main, p_patch );
  // Create the states vector based on the STATES array
  states_vec = N_VMake_Serial( p_patch->states_size, p_patch->STATES );
  // Initalize CVODE solver
  CVodeInit( cvode_main, rhs_fn, 0.0, states_vec );
  // Set up the future time (must be more than 0.0)
  tnext = dt;
  // Set up the maximum step size (dt)
  CVodeSetMaxStep( cvode_main, dt );
  // Set up the linear solver type
  CVDense( cvode_main, p_patch->states_size );
  // Set up the numerical error tolerances
  CVodeSStolerances( cvode_main, 1.0e-7, 1.0e-7 );
#endif
  fp_vm = fopen( "vmcheck.plt", "w");

  fprintf(fp_vm, "%s %s\n", "TIME", "Vm");

  int iout = 0;
  int imax = ( pace * p_patch->CONSTANTS[stim_period] ) / dt;
  time_t begin = time(NULL);

#if defined USE_CVODE
  while( iout <= imax ){
    retval = CVode( cvode_main, tnext, states_vec, &tcurr, CV_NORMAL  );
    if( retval == CV_SUCCESS ){
      iout++;
      tnext += dt;
    }
    else{
      printf("CVode calculation error at %lf msec!!!\n", tcurr);
      break;
    }
    fprintf(fp_vm, "%lf %lf\n", tcurr, p_patch->STATES[V]);

  }
  // Memory Cleanup
  N_VDestroy(states_vec);
  CVodeFree(&cvode_main);
#else
  while( iout <= imax ){
    p_patch->computeRates( tcurr,
      p_patch->CONSTANTS,
      p_patch->RATES,
      p_patch->STATES,
      p_patch->ALGEBRAIC );
    p_patch->solveAnalytical(dt);

    tcurr +=dt;
    iout++;

    fprintf(fp_vm, "%lf %lf\n", tcurr, p_patch->STATES[V]); 

  }

#endif

  time_t end = time(NULL);
  // end simulation loop
  printf("Time elapsed for simulation is: %ld sec\n", end - begin);


  fclose(fp_vm);

  delete p_patch;
  free(p_param);
  printf("FINISHED!\n");

  return 0;
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
    printf( "Hill File: %s\n", p_param->hill_file );
  }
  if(strlen(p_param->herg_file) != 0){
    printf( "hERG File: %s\n", p_param->herg_file );
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
    if (strcasecmp(key, "Celltype") == 0) {
      p_param->celltype = strtod( value, NULL );
      printf( "%s -- %lf\n", "Celltype", p_param->celltype);
    }
    else if (strcasecmp(key, "Variant") == 0) {
      strcpy( p_param->variant, value );
      printf( "%s -- %s\n", "Variant", p_param->variant);
    }
    else if (strcasecmp(key, "Simulation_Mode") == 0) {
      p_param->simulation_mode = strtol( value, NULL, 10 );
      printf( "%s -- %d\n", "Simulation_Mode", p_param->simulation_mode);
    }
    else if (strcasecmp(key, "Write_Time") == 0) {
      p_param->t_write_vtk = strtod( value, NULL );
      printf( "%s -- %lf\n", "Write_Time", p_param->t_write_vtk);
    }
    else if (strcasecmp(key, "Number_Of_Pacing") == 0) {
      p_param->num_pace1 = strtol( value, NULL, 10 );
      printf( "%s -- %d\n", "Number_of_Pacing", p_param->num_pace1);
    }
    else if (strcasecmp(key, "Basic_Cycle_Length") == 0) {
      p_param->bcl_init = strtod( value, NULL );
      printf( "%s -- %lf\n", "Basic_Cycle_Length", p_param->bcl_init);
    }
    else if (strcasecmp(key, "Time_Step") == 0) {
      p_param->dt = strtod( value, NULL );
      printf( "%s -- %lf\n", "Time_Step", p_param->dt);
    }
    else if (strcasecmp(key, "Write_Step") == 0) {
      p_param->dt_write = strtod( value, NULL );
      printf( "%s -- %lf\n", "Writing_Step", p_param->dt_write);
    }
    else if (strcasecmp(key, "Stim_Amplitude") == 0) {
      p_param->stim_amp = strtod( value, NULL );
      printf( "%s -- %lf\n", "Stim_Amplitude", p_param->stim_amp);
    }
    else if (strcasecmp(key, "Stim_Duration") == 0) {
      p_param->stim_dur = strtod( value, NULL );
      printf( "%s -- %lf\n", "Stim_Duration", p_param->stim_dur);
    }
    else if (strcasecmp(key, "BCL_Step") == 0) {
      p_param->bcl_decrement = strtod( value, NULL );
      printf( "%s -- %lf\n", "BCL_Step", p_param->bcl_decrement);
    }
    else if (strcasecmp(key, "BCL_End") == 0) {
      p_param->bcl_end = strtod( value, NULL );
      printf( "%s -- %lf\n", "BCL_End", p_param->bcl_end);
    }
    else if (strcasecmp(key, "APD_90_Percent") == 0) {
      p_param->v_apd90 = strtod( value, NULL );
      printf( "%s -- %lf\n", "APD_90_Percent", p_param->v_apd90);
    }
    else if (strcasecmp(key, "APD_50_Percent") == 0) {
      p_param->v_apd50 = strtod( value, NULL );
      printf( "%s -- %lf\n", "APD_50_Percent", p_param->v_apd50);
    }
    else if (strcasecmp(key, "Scaled_Gate") == 0) {
      p_param->scaled_gate = strtol( value, NULL, 10 );
      printf( "%s -- %d\n", "Scaled_Gate", p_param->scaled_gate);
    }
    else if (strcasecmp(key, "Drug_Name") == 0) {
      strcpy( p_param->drug_name, value );
      printf( "%s -- %s\n", "Drug_Name", p_param->drug_name);
    }
    else if (strcasecmp(key, "Concentrations") == 0) {
      strcpy( p_param->concs, "0," );
      strcat( p_param->concs, value );
      printf( "%s -- %s\n", "Concentrations", p_param->concs);
    }
    else if (strcasecmp(key, "GK1_Mult") == 0) {
      p_param->gk1_mult =  strtod(value, NULL) ;
      printf( "%s -- %lf\n", "GK1_Mult", p_param->gk1_mult);
    }
    else if (strcasecmp(key, "GCaL_Mult") == 0) {
      p_param->gcal_mult =  strtod(value, NULL) ;
      printf( "%s -- %lf\n", "GCaL_Mult", p_param->gcal_mult);
    }
    else if (strcasecmp(key, "GNa_Mult") == 0) {
      p_param->gna_mult =  strtod(value, NULL) ;
      printf( "%s -- %lf\n", "GNa_Mult", p_param->gna_mult);
    }
    else if (strcasecmp(key, "Is_Print_Vm") == 0) {
      p_param->is_print_vm = strtol( value, NULL, 10 );
      printf( "%s -- %d\n", "Is_Print_Vmcheck", p_param->is_print_vm );
    }

  }
  fclose( fp_inputdeck );
}

