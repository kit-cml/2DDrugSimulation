#include "patches/tentusscher_noble_noble_panfilov_2004_b.hpp"

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
void create_apdr_output(int num_pace1, int number_of_lines);
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
  FILE *fp_curr_all;
  FILE *fp_vm;
  FILE *fp_vm_all;
  FILE *fp_conc;
  FILE *fp_conc_all;
  FILE *fp_apd90;
  FILE *fp_apdr;
  FILE *fp_apdur;

  // Patch clamp variables
  patch_clamp *p_patch;
  double bcl,dt,pace;

  // SUNDIALs variables
  int retval;
  double tnext, t_curr;
  void *cvode_main;
  N_Vector states_vec;

  // Support variables
  char file_name[255];
  double t_prev;
  double t_stim_start;
  double t_stim_end;
  double t_di;
  double t_di_prev;
  double t_depolarize;
  double t_depolarize_prev;
  double t_repolarize;
  double bcl_decrement;
  bool is_ap_increasing;
  bool is_action_potential;
  int depolarization_counter;
  int apdr_file_lines;
  double apd;
  double v_prev;
  double v_top;
  double v_valley;
  double v_apd90;



  param_t *p_param;
  p_param = (param_t*)malloc( sizeof( param_t ) );
  edison_assign_params_single(argc,argv,p_param);

  p_patch = new tentusscher_noble_noble_panfilov_2004_b();
  p_patch->initConsts();
  if((int)floor(p_param->bcl_init) == 0)bcl = 1000.0;
  else bcl = p_param->bcl_init;
  p_patch->CONSTANTS[stim_period] = bcl;
  if(p_param->dt > 0)dt = p_param->dt;
  else dt = 0.01;
  if(p_param->num_pace1 == 0) pace = 3;
  else pace = p_param->num_pace1;
  t_curr = 0.0;

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

  bcl = p_param->bcl_init;
  bcl_decrement = p_param->bcl_decrement;
  t_curr = 0.0;
  t_prev = 0.0;
  t_stim_start = 0.0;
  t_stim_end = t_stim_start + p_patch->CONSTANTS[stim_duration];
  t_di = 0.0;
  t_di_prev = 0.0;
  t_depolarize_prev = 0.0;
  t_depolarize = 0.0;
  apd = 0.0;
  depolarization_counter = 0;
  apdr_file_lines = 0;
  v_prev = 0.0;
  v_top = 0.0;
  v_valley = 0.0;
  v_apd90 = 0.0;


  fp_apdr = fopen( "apdr.plt", "w" );
  fp_apd90 = fopen( "debug.plt", "w" );
  sprintf(file_name,"ires_%.0lf.plt", bcl);
  fp_curr = fopen(file_name, "w");
  fp_curr_all = fopen("ires.plt", "w");
  sprintf(file_name,"vmcheck_%.0lf.plt", bcl);
  fp_vm = fopen( file_name, "w");
  fp_vm_all = fopen( "vmcheck.plt", "w");
  sprintf(file_name,"conc_%.0lf.plt", bcl);
  fp_conc = fopen( file_name, "w");
  fp_conc_all = fopen( "conc.plt", "w");
  fp_apdur = fopen("apdur_log.plt", "w");

  fprintf(fp_vm, "%s %s\n", "TIME", "Vm");
  fprintf(fp_vm_all, "%s %s\n", "TIME", "Vm");
  fprintf(fp_curr, "%s %s %s %s %s %s %s %s\n",
           "TIME","INa","ICaL", "INaCa", "Ito","IKr","IKs","IK1");
  fprintf(fp_curr_all, "%s %s %s %s %s %s %s %s\n",
           "TIME","INa","ICaL", "INaCa","Ito","IKr","IKs","IK1");
  fprintf(fp_conc, "%s %s %s\n", "TIME", "Na_i", "Ca_i");
  fprintf(fp_conc_all, "%s %s %s\n", "TIME", "Na_i", "Ca_i");

  fprintf(fp_apdur, "%s %s %s %s %s %s\n",
          "Depol_Count", "T_Depol", "APD", "V_APD90", "DELTA","BCL");

  int iout = 0;
  int imax = ( pace * p_patch->CONSTANTS[stim_period] ) / dt;
  time_t begin = time(NULL);

  printf("BCL: %lf\n", bcl);
  while( bcl >= p_param->bcl_end ){
    v_prev = p_patch->STATES[0];
    if (bcl <= 300.0) bcl_decrement = 10;
    if (t_curr <= t_prev + bcl * p_param->num_pace1) {
      if (t_curr > t_stim_end) {
/*        printf( "%s = %d %s = %lf %s = %lf %s = %lf\n",
                "depolarization_counter", depolarization_counter,
                "t_stim_start", t_stim_start,
                "t_stim_end", t_stim_end,
                "bcl", bcl );
*/
        t_stim_start += bcl;
        t_stim_end = t_stim_start + p_patch->CONSTANTS[stim_duration];
      }
    }
    else {
      //printf("BCL: %lf LAST APD: %lf\n", bcl, apd);
      bcl -= p_param->bcl_decrement;
      p_patch->CONSTANTS[stim_period] = bcl;

      fclose(fp_vm);
      fclose(fp_curr);
      fclose(fp_conc);
      sprintf(file_name,"ires_%.0lf.plt", bcl);
      fp_curr = fopen(file_name, "w");
      sprintf(file_name,"vmcheck_%.0lf.plt", bcl);
      fp_vm = fopen( file_name, "w");
      sprintf(file_name,"conc_%.0lf.plt", bcl);
      fp_conc = fopen( file_name, "w");

      fprintf(fp_vm, "%s %s\n", "TIME", "Vm");
      fprintf(fp_vm_all, "%s %s\n", "TIME", "Vm");
      fprintf(fp_curr, "%s %s %s %s %s %s %s %s\n",
           "TIME","INa","ICaL", "INaCa","Ito","IKr","IKs","IK1");
      fprintf(fp_curr_all, "%s %s %s %s %s %s %s %s\n",
           "TIME","INa","ICaL", "INaCa","Ito","IKr","IKs","IK1");
      fprintf(fp_conc, "%s %s %s\n", "TIME", "Na_i", "Ca_i");
      fprintf(fp_conc_all, "%s %s %s\n", "TIME", "Na_i", "Ca_i");

      t_prev = t_curr;
      depolarization_counter = 0;
    }

    /* solve ODE */
    retval = CVode( cvode_main, tnext, states_vec, &t_curr, CV_NORMAL  );
    if( retval == CV_SUCCESS ){
      iout++;
      tnext += dt;
    }
    else{
      printf("CVode calculation error at %lf msec!!!\n", t_curr);
      break;
    }

    /* find out whether the simulation can produce ACTION POTENTIAL or not */
    if (p_patch->STATES[0] > v_prev && !is_ap_increasing) {
      is_ap_increasing = true;
      v_valley = v_prev;
      is_action_potential = v_top - v_valley > v_apd90;
      if (is_action_potential && v_prev < 0) {
        v_apd90 = v_top - (0.9 * (v_top - v_valley));
        fprintf( fp_apd90, "%lf %lf %lf %lf %lf\n",
                 round(t_curr), v_top, v_valley, v_apd90, bcl );
      }
    }
    else if (p_patch->STATES[0] < v_prev && is_ap_increasing) {
      is_ap_increasing = false;
      v_top = v_prev;
    }


    // repolarization, when the potential moves from HIGH to LOW
    if ( v_prev > v_apd90 && p_patch->STATES[0] <= v_apd90 ) {
      t_repolarize = t_curr;
      apd = t_repolarize - t_depolarize;
      fprintf(fp_apdur, "%d %lf %lf %lf %lf\n", depolarization_counter, t_depolarize, apd, v_apd90, bcl);
    }
    // depolarization, when the potential moves from LOW to HIGH
    if ( v_prev <= v_apd90 && p_patch->STATES[0] > v_apd90 ) {
      depolarization_counter++;
      apdr_file_lines++;
      t_depolarize_prev = t_depolarize;
      t_depolarize = t_curr;
      if (depolarization_counter >= 3) {
        t_di_prev = t_di;
        t_di = t_depolarize - t_repolarize;
        if(depolarization_counter > 3 && t_di-t_di_prev > 100.){
          apd = t_curr - t_depolarize_prev;
        }
        fprintf( fp_apdr, "%d %lf %lf %lf %lf\n",
                 depolarization_counter, bcl, apd, t_di, t_curr);
      }
    }


    fprintf(fp_vm, "%lf %lf\n", t_curr, p_patch->STATES[V]);
    fprintf(fp_vm_all, "%lf %lf\n", t_curr, p_patch->STATES[V]);
    fprintf(fp_curr, "%lf %lf %lf %lf %lf %lf %lf %lf\n",
           t_curr, p_patch->ALGEBRAIC[i_Na], p_patch->ALGEBRAIC[i_CaL], p_patch->ALGEBRAIC[i_NaCa], p_patch->ALGEBRAIC[i_to], p_patch->ALGEBRAIC[i_Kr], p_patch->ALGEBRAIC[i_Ks], p_patch->ALGEBRAIC[i_K1]);
    fprintf(fp_curr_all, "%lf %lf %lf %lf %lf %lf %lf %lf\n",
           t_curr, p_patch->ALGEBRAIC[i_Na], p_patch->ALGEBRAIC[i_CaL], p_patch->ALGEBRAIC[i_NaCa], p_patch->ALGEBRAIC[i_to], p_patch->ALGEBRAIC[i_Kr], p_patch->ALGEBRAIC[i_Ks], p_patch->ALGEBRAIC[i_K1]);
    fprintf(fp_conc, "%lf %lf %f\n", t_curr, p_patch->STATES[Na_i], p_patch->STATES[Ca_i]);
    fprintf(fp_conc_all, "%lf %lf %f\n", t_curr, p_patch->STATES[Na_i], p_patch->STATES[Ca_i]);


  }
  // Memory Cleanup
  N_VDestroy(states_vec);
  CVodeFree(&cvode_main);


  time_t end = time(NULL);
  // end simulation loop
  printf("Time elapsed for simulation is: %ld sec\n", end - begin);


  fclose(fp_vm_all);
  fclose(fp_curr_all);
  fclose(fp_conc_all);
  fclose( fp_apd90 );
  fclose( fp_apdr );
  fclose( fp_apdur );


  delete p_patch;
  free(p_param);
  /* generate APD vs DiastolicInterval and APD vs BCL result files */
  create_apdr_output(p_param->num_pace1, apdr_file_lines + 1);
 

 printf("FINISHED!\n");

  return 0;
}

void create_apdr_output(int num_pace1, int number_of_lines)
{
  long int CP = 0;

  int* icnt = (int*)malloc(number_of_lines * sizeof(int));

  printf("Entering APDR\n%s: %d\n", "Number of Lines", number_of_lines);

  double time;
  double* ap_duration_arr = (double*)malloc(number_of_lines * sizeof(double));
  double* di_arr = (double*)malloc(number_of_lines * sizeof(double));
  double* bcl_arr = (double*)malloc(number_of_lines * sizeof(double));

  FILE *fp_apdr;
  FILE *fp_apdr_last_2_cycles;
  FILE *fp_apdr_di;
  FILE *fp_apdr_bcl;

  fp_apdr = fopen( "apdr.plt", "r" );
  fp_apdr_last_2_cycles = fopen( "apdr2.plt", "w" );
  fp_apdr_di = fopen( "apdr_di.plt", "w" );
  fp_apdr_bcl = fopen( "apdr_bcl.plt", "w" );

  fscanf(fp_apdr, "%d\t%lf\t%lf\t%lf\t%lf\n",
         &icnt[1], &bcl_arr[1], &ap_duration_arr[1], &di_arr[1], &time);
  fprintf(fp_apdr_last_2_cycles, "%s %s %s %s %s\n",
          "bcl", "apd90-1", "apd90-2", "di-1", "di-2");
  fprintf(fp_apdr_bcl, "%s %s %s\n", "bcl", "apd90-1", "apd90-2");
  fprintf(fp_apdr_di, "%s %s %s\n", "di", "apd90-1", "apd90-2");

  for (int i = 2; i < number_of_lines; i++) {
    fscanf(fp_apdr, "%d\t%lf\t%lf\t%lf\t%lf\n",
           &icnt[i], &bcl_arr[i], &ap_duration_arr[i], &di_arr[i], &time);
    if (di_arr[i] == -1) break;
    if (bcl_arr[i] != bcl_arr[i - 1] && (int)bcl_arr[i] != 0 ) {
      CP += bcl_arr[i - 2] * num_pace1;
      fprintf(fp_apdr_last_2_cycles, "%lf %lf %lf %lf %lf %ld\n",
              bcl_arr[i - 2], ap_duration_arr[i - 3],
              ap_duration_arr[i - 2], di_arr[i - 3], di_arr[i - 2], CP);
      fprintf(fp_apdr_di, "%lf %lf %lf\n",
              di_arr[i - 3], ap_duration_arr[i - 3], ap_duration_arr[i - 2]);
      fprintf(fp_apdr_bcl, "%lf %lf %lf\n",
              bcl_arr[i - 3], ap_duration_arr[i - 3], ap_duration_arr[i - 2]);
    }
  }

  fclose(fp_apdr);
  fclose(fp_apdr_last_2_cycles);
  fclose(fp_apdr_di);
  fclose(fp_apdr_bcl);

  free(icnt);
  free(ap_duration_arr);
  free(di_arr);
  free(bcl_arr);
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
    else if (strcasecmp(key, "APD_90_Percent_Max") == 0) {
      p_param->v_apd90_max = strtod( value, NULL );
      printf( "%s -- %lf\n", "APD_90_Percent_Max", p_param->v_apd90_max);
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
    else if (strcasecmp(key, "GKs_Mult") == 0) {
      p_param->gks_mult =  strtod(value, NULL) ;
      printf( "%s -- %lf\n", "GKs_Mult", p_param->gks_mult);
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

