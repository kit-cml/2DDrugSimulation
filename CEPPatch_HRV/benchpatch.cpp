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


void edison_assign_params_single(int argc, char *args[], param_t *p_param);
int rhs_fn( realtype t, N_Vector y, N_Vector ydot, void *user_data );

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
  char concs_str[255];
  char filename[255];
  char *pch;
  FILE *fp_curr;
  FILE *fp_vm;
  FILE *fp_hrv;
  FILE *fp_hill;
  FILE *fp_qni;

  // Patch clamp variables
  patch_clamp *p_patch;
  double bcl,dt;

  // Drug related variables
  bool is_drug;
  double hill[14];
  double conc;
  double qnet_prev, inal_auc_prev, ical_auc_prev, inal_auc_ctl, ical_auc_ctl,qinward;

  // SUNDIALs variables
  int retval;
  double tnext, tcurr;
  void *cvode_main;
  N_Vector states_vec;

  // Supporting variables
  int iout;
  int imax;
  int idx;

  param_t *p_param;
  p_param = (param_t*)malloc( sizeof( param_t ) );
  edison_assign_params_single(argc,argv,p_param);

  fp_hill = fopen(p_param->hill_file,"r");
  if(fp_hill == 0){
    printf("Hill data cannot be opened!!!\n");
    exit(0);
  }

  /* Skip the header in the sample file */
  fgets(buff, sizeof(buff), fp_hill);
  /* Read the actual data and put it inside buff string */
  fgets(buff, sizeof(buff), fp_hill);
  pch = strtok(buff, ",");
  if(pch == NULL){
    printf("Incorrect drug sample file. Program aborted!!\n");
    exit(0);
  }
  else{
    idx = 0;
    /* begin splitting sample loop */
    do
    {
      hill[idx] = strtod(pch, NULL);
      printf("hill[%d] = %lf\n", idx, hill[idx]);
      pch = strtok(NULL, ",");
      idx++;
    }
    while( pch != NULL );
    /* end splitting sample loop */
  }

  // assigning concentration value
  strcat( concs_str, "0," );
  strcat( concs_str, p_param->concs );

  int print_freq = (1. / p_param->dt) * 2.0;
  pch = strtok(concs_str, ",");

  time_t begin = time(NULL);
  do{

    conc = strtod(pch, NULL);
    printf("Concentration: %lf\n", conc);
    p_patch = new Ohara_Rudy_2011();
    if( is_drug == true ) p_patch->initConsts(0, conc, hill, true);
    else p_patch->initConsts();
    dt = 0.01;
    tcurr = 0.0;

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

    qnet_prev = 0.0;
    inal_auc_prev = 0.0;
    ical_auc_prev = 0.0;

    sprintf(filename, "vmcheck_%.0lfnM.plt", conc);
    filename[sizeof(filename) - 1] = '\0';
    fp_vm = fopen( filename, "w");
    sprintf(filename, "qnet_qinward_%.0lfnM.plt", conc);
    filename[sizeof(filename) - 1] = '\0';
    fp_qni = fopen( filename, "a");
    fp_hrv = fopen( p_param->bcl_file, "r");
    if(fp_hrv == NULL){
      printf("Error opening BCL list file");
      exit(0);
    };

    fprintf(fp_vm, "%s %s\n", "TIME", "Vm");
    fprintf(fp_qni, "%s %s %s\n", "BCL", "Qnet", "Qinward");

    while( fgets( buff, sizeof(buff), fp_hrv ) != NULL){

      p_patch->CONSTANTS[stim_period] = strtod(buff, NULL);
      if(p_patch->CONSTANTS[stim_period] < 10E-14){
        printf("WRONG BCL VALUE!!!\n");
        continue;
      }
      printf("BCL: %s\n", buff);
      iout = 0;
      imax = p_patch->CONSTANTS[stim_period] / dt;
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
        if(iout % print_freq == 0 ){
          fprintf(fp_vm, "%lf %lf\n", tcurr, p_patch->STATES[V]);
        }
      }


      if((int)ceil(conc) == 0){
        inal_auc_ctl = p_patch->STATES[INaL_AUC]-inal_auc_prev;
        ical_auc_ctl = p_patch->STATES[ICaL_AUC]-ical_auc_prev;
      }
      else{
        qinward = ( ((p_patch->STATES[INaL_AUC]-inal_auc_prev)/inal_auc_ctl) + ((p_patch->STATES[ICaL_AUC]-ical_auc_prev)/ical_auc_ctl) ) * 0.5;
      }
      if((int)ceil(conc) != 0)  fprintf(fp_qni, "%.2lf %lf %lf\n", p_patch->CONSTANTS[stim_period], (p_patch->STATES[qnet] - qnet_prev)/1000.0, qinward);
      else   fprintf(fp_qni, "%.2lf %lf %.2lf\n", p_patch->CONSTANTS[stim_period], (p_patch->STATES[qnet] - qnet_prev)/1000.0, 0.);
      p_patch->CONSTANTS[last_stim_period] += p_patch->CONSTANTS[stim_period];
      qnet_prev = p_patch->STATES[qnet];
      inal_auc_prev = p_patch->STATES[INaL_AUC];
      ical_auc_prev = p_patch->STATES[ICaL_AUC];
    }
    // CVode Memory Cleanup
    N_VDestroy(states_vec);
    CVodeFree(&cvode_main);
    // FILE pointer cleanup
    fclose(fp_vm);
    fclose(fp_hrv);
    fclose(fp_qni);
    if( fp_hill != NULL )fclose(fp_hill);

    delete p_patch;

    pch = strtok(NULL, ",");
  }while(pch != NULL);

  printf("FINISHED!\n");

  time_t end = time(NULL);
  printf("Time elapsed for simulation is: %ld sec\n", end - begin);

  free(p_param);

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
    else if (!strcmp(args[i], "-bcl_file"))
      strcpy(p_param->bcl_file, args[i + 1]);
  }

  if(strlen(p_param->hill_file) != 0){
    printf( "Hill File: %s\n", p_param->hill_file );
  }
  if(strlen(p_param->bcl_file) != 0){
    printf( "BCL File: %s\n", p_param->bcl_file );
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
    else if (strcasecmp(key, "GKs_Scale") == 0) {
      p_param->gks_scale =  strtod(value, NULL) ;
      printf( "%s -- %lf\n", "GKs_Scale", p_param->gks_scale);
    }
    else if (strcasecmp(key, "GCaL_Scale") == 0) {
      p_param->gcal_scale =  strtod(value, NULL) ;
      printf( "%s -- %lf\n", "GCaL_Scale", p_param->gcal_scale);
    }
    else if (strcasecmp(key, "GNa_Scale") == 0) {
      p_param->gna_scale =  strtod(value, NULL) ;
      printf( "%s -- %lf\n", "GNa_Scale", p_param->gna_scale);
    }
    else if (strcasecmp(key, "GNaL_Scale") == 0) {
      p_param->gnal_scale =  strtod(value, NULL) ;
      printf( "%s -- %lf\n", "GNaL_Scale", p_param->gnal_scale);
    }
    else if (strcasecmp(key, "Is_Print_Vm") == 0) {
      p_param->is_print_vm = strtol( value, NULL, 10 );
      printf( "%s -- %d\n", "Is_Print_Vmcheck", p_param->is_print_vm );
    }

  }
  fclose( fp_inputdeck );
}
