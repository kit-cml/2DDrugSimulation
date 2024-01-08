#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <cvode/cvode.h>
#include <cvode/cvode_dense.h>
#include <nvector/nvector_serial.h>

#define IS_EDISON 1
#include "commons/helper.hpp"

static int rhs_fn( realtype t, N_Vector y, N_Vector ydot, void *user_data );
static double average( double *apd_arr, int start_idx, int end_idx );
static double get_ANM( double *arr, int length );
static void create_apdr_output(int num_pace1, int number_of_lines);
static int AP(int argc, char **argv, param_t *p_param);
static int APDR(int argc, char **argv, param_t *p_param);
static int ANM(int argc, char **argv, param_t *p_param);
static int gate_variation_AP(int argc, char **argv, param_t *p_param);
#if (defined TN2006ENDO) || (defined TN2006M) || (defined TN2006M)
static int multi_population_AP(int argc, char **argv, param_t *p_param);
#endif
#if (defined OHARA_RUDY2011)
static int population_hill( int argc, char **argv, param_t *p_param, bool is_dutta );
static int population_hill_ANM( int argc, char **argv, param_t *p_param, bool is_dutta );
#endif
#if (defined ORUDY_CIPA2017)
static int population_herg( int argc, char **argv, param_t *p_param );
static int population_herg_ANM( int argc, char **argv, param_t *p_param );
#endif

int main( int argc, char **argv )
{
  int err_code;
  param_t *p_param;


  p_param = (param_t*)malloc( sizeof( param_t ) );
  err_code = 0;

  set_default_values(p_param);

 
#if IS_EDISON == 1
  edison_assign_params_single(argc, argv, p_param);
#else
//  assign_params(argc, argv, p_param);
#endif

  if( p_param->simulation_mode == 0 ){
    err_code = AP( argc, argv, p_param );
  }
#if defined TN2006ENDO || defined TN2006EPI || defined TN2006M
  else if( p_param->simulation_mode == 1 ){
    err_code = multi_population_AP( argc, argv, p_param );
  }
#endif
  else if( p_param->simulation_mode == 2 ){
    err_code = gate_variation_AP( argc, argv, p_param );
  }
  else if( p_param->simulation_mode == 3 ){
    err_code = APDR( argc, argv, p_param );
  }
  else if( p_param->simulation_mode == 13 ){
    err_code = ANM( argc, argv, p_param );
  }
#if defined OHARA_RUDY2011
  else if( p_param->simulation_mode == 4 ){
    err_code = population_hill( argc, argv, p_param, false);
  }
  else if( p_param->simulation_mode == 8 ){
    err_code = population_hill( argc, argv, p_param, true );
  }
  else if( p_param->simulation_mode == 9 ){
    err_code = population_hill_ANM( argc, argv, p_param, false );
  }
#endif
#if defined ORUDY_CIPA2017
  else if( p_param->simulation_mode == 5 ){
    err_code = population_herg( argc, argv, p_param );
  }
  else if( p_param->simulation_mode == 10 ){
    err_code = population_herg_ANM( argc, argv, p_param );
  }
#endif
  else{
    printf("Selection unknown!!\n");
  }


  return err_code;
}

#if (defined OHARA_RUDY2011)
static int population_hill( int argc, char **argv, param_t *p_param, bool is_dutta )
{
  /* CellML object pointer */
  CellML *p_cell;
  /* time related variables */
  int print_freq;
  int err_code;
  int i;
  char file_name[100];
  FILE *fp_apd;
  FILE *fp_vm;
  FILE *fp_ca;
  FILE *fp_ires;

  /* CiPA related variables */
  FILE *fp_hill;
  FILE *fp_qnet;
  FILE *fp_qinward;
  double hill[14];
  double dosage;
  double conc;
  double cmax;
  char *pch;
  char dosages_str[100];
  char hill_sample[128];
  char result_file[128];
  int sim_id;
  double qnet_prev;
  double qnet_curr;
  double INaL_auc_prev;
  double INaL_auc_curr;
  double ICaL_auc_prev;
  double ICaL_auc_curr;
  double INaL_auc_control;
  double ICaL_auc_control;
  double INaL_auc_drug;
  double ICaL_auc_drug;
  double qinward;

  /* APD related variables */
  bool is_ap_increasing;
  bool is_action_potential;
  bool is_apd_written;
  double apd;
  double v_prev;
  double v_top;
  double v_valley;
  double v_apd90;
  double t_depolarize;
  double t_repolarize;


  /* SUNDIALs variables */
  int retval;
  double t_out, t;
  void *cvode_mem;
  N_Vector states_vec;

  /* Initialize p_cell based on the ionmodel value */
  p_cell = init_cell( );
  if( p_cell == NULL ){
    fprintf(stderr, "Problem when initializing cellmodel\n");
    return EXIT_FAILURE;
  }

  /* Initialize mutation object of p_cell if possible */
  p_cell->isMutated = 0;
  if( strcmp(p_param->variant, "ori") != 0 ){
    p_cell = init_mutation( p_param->variant, 0.5, p_cell);
    if( p_cell->mutation != NULL ){
      p_cell->isMutated = 1;
    }
  }


  /* Open the sample file */
  fp_hill = fopen(p_param->hill_file, "r");
  if(fp_hill == NULL){
    fprintf(stderr, "Cannot open file %s\n", p_param->hill_file);
    return EXIT_FAILURE;
  }

  /* Define the value of cmax based on the drug name */
  if( (int)round(p_param->cmax) == 0 ){
    cmax = get_cmax( p_param->drug_name );
  }
  else{
    cmax = p_param->cmax;
  }

  /* Skip the header in the sample file */
  fgets(hill_sample, sizeof(hill_sample), fp_hill);

  /* begin population loop */
  sim_id = 1; 
  while(fgets(hill_sample, sizeof(hill_sample), fp_hill) != NULL)
  {
    i = 0;
    pch = strtok(hill_sample, ",");
    /* begin splitting sample loop */
    do
    {
      hill[i++] = strtod(pch, NULL);
      pch = strtok(NULL, ",");
    }
    while( pch != NULL );
    /* end splitting sample loop */

    strncpy(dosages_str, p_param->dosages, sizeof(dosages_str));
    pch = strtok(dosages_str, ",");
    /* begin splitting the dosages_str string loop */
    do{
      dosage = strtod(pch, NULL);
      conc = cmax * dosage;

      p_cell->initConsts(p_param->celltype, conc, hill, is_dutta, p_param->drug_name);
      if( p_param->is_using_output == 1 ){
        load_last_state(p_cell);
      }
      p_cell->isS1 = 1;
      p_cell->CONSTANTS[stim_period] = p_param->bcl_init;
      p_cell->CONSTANTS[stim_amplitude] *= p_param->stim_amp;

      cvode_mem = CVodeCreate( CV_BDF,CV_NEWTON );
      CVodeSetUserData( cvode_mem, p_cell );
      states_vec = N_VMake_Serial( p_cell->states_size, p_cell->STATES );
      CVodeInit( cvode_mem, rhs_fn, T0, states_vec );
      t_out = T1;
      CVodeSetMaxStep( cvode_mem, p_param->dt );
      CVDense( cvode_mem, p_cell->states_size );
      CVodeSStolerances( cvode_mem, 1.0e-7, 1.0e-7 );

      int iout = 0;
      int imax = ( p_param->num_pace1 *  p_cell->CONSTANTS[stim_period] ) / p_param->dt;
      print_freq = (1. / p_param->dt) * p_param->dt_write;
      if( p_param->is_print_vm == 1 ){
        sprintf(result_file, "%s_%.1lf_vmcheck_#%d.plt", p_param->drug_name, dosage, sim_id );
        result_file[sizeof(result_file) - 1] = '\0';
        fp_vm = fopen( result_file, "w" );
        sprintf(result_file, "%s_%.1lf_ca_i_#%d.plt", p_param->drug_name, dosage, sim_id );
        result_file[sizeof(result_file) - 1] = '\0';
        fp_ca = fopen( result_file, "w" );
        sprintf(result_file, "%s_%.1lf_ires_#%d.plt", p_param->drug_name, dosage, sim_id );
        result_file[sizeof(result_file) - 1] = '\0';
        fp_ires = fopen( result_file, "w" );
      }
      sprintf(result_file, "%s_%.1lf_apd.plt", p_param->drug_name, dosage );
      result_file[sizeof(result_file) - 1] = '\0';
      fp_apd = fopen( result_file, "a" );
      sprintf(result_file, "%s_%.1lf_qnet.plt", p_param->drug_name, dosage );
      result_file[sizeof(result_file) - 1] = '\0';
      fp_qnet = fopen( result_file, "a" );
      sprintf(result_file, "%s_%.1lf_qinward.plt", p_param->drug_name, dosage );
      result_file[sizeof(result_file) - 1] = '\0';
      fp_qinward = fopen( result_file, "a" );

     
      is_apd_written = false;
      is_ap_increasing = false;
      is_action_potential = false;
      v_prev = 0.0;
      v_top = 0.0;
      v_valley = 0.0;
      v_apd90 = p_param->v_apd90;

#ifdef DEBUG_COMPLETE
      if( p_param->is_print_vm == 1 ){
        fprintf( fp_vm, "%s %s\n", 
                 "TIME", 
                 "Vm" );
      fprintf( fp_ca, "%s %s\n", 
               "TIME", 
               "cai" );

      fprintf( fp_ires, "%s %s %s %s %s %s %s %s %s\n",
               "TIME",
               "INaL",
               "ICaL",
               "Ito",
               "IKr",
               "IKs",
               "IK1",
               "INaCa",
               "INaK" );
      }
#endif
      qnet_prev = 0.0;
      qnet_curr = 0.0;
      INaL_auc_prev = 0.0;
      ICaL_auc_prev = 0.0;
      INaL_auc_curr = 0.0;
      ICaL_auc_curr = 0.0;

      /* begin simulation loop */
      while(1){
        v_prev = p_cell->STATES[V];

        retval = CVode( cvode_mem, t_out, states_vec, &t, CV_NORMAL  );
        if( retval == CV_SUCCESS ){
          iout++;
          t_out += p_param->dt;
        }
        if (iout >= imax) break;


        // find out whether the simulation can produce ACTION POTENTIAL or not
        if (p_cell->STATES[V] > v_prev && !is_ap_increasing) {
          is_ap_increasing = true;
          v_valley = v_prev;
          is_action_potential = v_top - v_valley > -40;
          if (is_action_potential) {
            //v_apd90 = v_top - (0.9 * (v_top - v_valley));
            //v_apd90 = -75;
          }
        }
        else if (p_cell->STATES[V] < v_prev && is_ap_increasing) {
          is_ap_increasing = false;
          v_top = v_prev;
        }


        // repolarization, when the potential moves from HIGH to LOW
        if (v_prev > v_apd90 && p_cell->STATES[V] <= v_apd90) {
          t_repolarize = t;
          apd = t_repolarize - t_depolarize;
          if( t >= (p_param->num_pace1 * p_param->bcl_init) - ( 1 * p_param->bcl_init ) && is_apd_written == false){
            fprintf( fp_apd, "%lf\n", apd  );
           is_apd_written = true;
          }
        }
        // depolarization, when the potential moves from LOW to HIGH
        if (v_prev <= v_apd90 && p_cell->STATES[V] > v_apd90) {
          t_depolarize = t;
       }

#ifdef DEBUG_COMPLETE
       if(iout % print_freq == 0){
         if( p_param->is_print_vm == 1 ){
           fprintf( fp_vm, "%lf %lf\n",
                    t,
                    p_cell->STATES[V] );
           fprintf( fp_ca, "%lf %lf\n",
                    t,
                    p_cell->STATES[Ca_i] );
           fprintf( fp_ires, "%lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
                    t,
                    p_cell->ALGEBRAIC[i_NaL],
                    p_cell->ALGEBRAIC[i_CaL],
                    p_cell->ALGEBRAIC[i_to],
                    p_cell->ALGEBRAIC[i_Kr],
                    p_cell->ALGEBRAIC[i_Ks],
                    p_cell->ALGEBRAIC[i_K1],
                    p_cell->ALGEBRAIC[i_NaCa_i],
                    p_cell->ALGEBRAIC[i_NaK] );
         }
         save_last_state(p_cell, t);
       }
#endif
           // save qnet and auc(area under curve) for the last and before last cycles
       if( iout % (int)(p_cell->CONSTANTS[stim_period]/p_param->dt) == 0 ){
           if( iout / (int)(p_cell->CONSTANTS[stim_period]/p_param->dt) 
              == p_param->num_pace1-2 ){
              qnet_prev = p_cell->STATES[qnet];
              INaL_auc_prev = p_cell->STATES[INaL_AUC];
              ICaL_auc_prev = p_cell->STATES[ICaL_AUC];
              printf("TEST1: %d %lf %lf %lf\n", iout, qnet_prev, INaL_auc_prev, ICaL_auc_prev);
           }
           else if( iout / (int)(p_cell->CONSTANTS[stim_period]/p_param->dt) 
              == p_param->num_pace1-1 ){
              qnet_curr = p_cell->STATES[qnet];
              INaL_auc_curr = p_cell->STATES[INaL_AUC];
              ICaL_auc_curr = p_cell->STATES[ICaL_AUC];
              printf("TEST2: %d %lf %lf %lf\n", iout, qnet_curr, INaL_auc_curr, ICaL_auc_curr);
           }
       }


    }
    /* end simulation loop */
  
    /* saves qnet result to the file */
    fprintf( fp_qnet, "%lf\n", (qnet_curr - qnet_prev)/1000.0 );

    /* calculate qinward for dosages > 0. */
    if((int)round(dosage) == 0){
      INaL_auc_control = INaL_auc_curr - INaL_auc_prev;
      ICaL_auc_control = ICaL_auc_curr - ICaL_auc_prev;
    }
    else{
      INaL_auc_drug = INaL_auc_curr - INaL_auc_prev;
      ICaL_auc_drug = ICaL_auc_curr - ICaL_auc_prev;
      qinward = ( (INaL_auc_drug/INaL_auc_control) + (ICaL_auc_drug/ICaL_auc_control) ) * 0.5;
      fprintf( fp_qinward, "%lf\n", qinward );
    }

    
    fclose( fp_qinward );
    fclose( fp_qnet );
    fclose( fp_apd );
#ifdef DEBUG_COMPLETE
    if(fp_ires != NULL) fclose( fp_ires );
    if(fp_vm != NULL) fclose( fp_vm );
    if(fp_ca != NULL)fclose( fp_ca );
#endif

    N_VDestroy(states_vec);
    CVodeFree(&cvode_mem);

    pch = strtok(NULL, ","); 
    }
    while( pch != NULL );
    /* end splitting the dosages_str string loop */
    sim_id++;

 
  }
  /* end population loop */

  // generate EDISON output
  system("rm -rf result");
  system("mkdir result");
  system("mv *.plt result");
  system("mv output.dat result");

  delete p_cell;
  return 0;
}
#endif





#if (defined ORUDY_CIPA2017)
static int population_herg( int argc, char **argv, param_t *p_param )
{
  /* CellML object pointer */
  CellML *p_cell;
  /* time related variables */
  int print_freq;
  int err_code;
  int i;
  char file_name[100];
  FILE *fp_apd;
  FILE *fp_vm;
  FILE *fp_ca;
  FILE *fp_ires;

  /* CiPA related variables */
  FILE *fp_herg;
  FILE *fp_hill;
  FILE *fp_qnet;
  FILE *fp_qinward;
  double herg[20];
  double dosage;
  double conc;
  double cmax;
  char *pch;
  char dosages_str[100];
  char hill_sample[128];
  char herg_sample[128];
  char result_file[128];
  int sim_id;
  double qnet_prev;
  double qnet_curr;
  double INaL_auc_prev;
  double INaL_auc_curr;
  double ICaL_auc_prev;
  double ICaL_auc_curr;
  double INaL_auc_control;
  double ICaL_auc_control;
  double INaL_auc_drug;
  double ICaL_auc_drug;
  double qinward;

  /* APD related variables */
  bool is_ap_increasing;
  bool is_action_potential;
  bool is_apd_written;
  double apd;
  double v_prev;
  double v_top;
  double v_valley;
  double v_apd90;
  double t_depolarize;
  double t_repolarize;


  /* SUNDIALs variables */
  int retval;
  double t_out, t;
  void *cvode_mem;
  N_Vector states_vec;

  /* Initialize p_cell based on the ionmodel value */
  p_cell = init_cell( );
  if( p_cell == NULL ){
    fprintf(stderr, "Problem when initializing cellmodel\n");
    return EXIT_FAILURE;
  }

  /* Initialize mutation object of p_cell if possible */
  p_cell->isMutated = 0;
  if( strcmp(p_param->variant, "ori") != 0 ){
    p_cell = init_mutation( p_param->variant, 0.5, p_cell);
    if( p_cell->mutation != NULL ){
      p_cell->isMutated = 1;
    }
  }


  /* Open the sample file */
  fp_herg = fopen(p_param->herg_file, "r");
  if(fp_herg == NULL){
    fprintf(stderr, "Cannot open file %s\n", p_param->herg_file);
    return EXIT_FAILURE;
  }
  /* Open the sample file */
  fp_hill = fopen(p_param->hill_file, "r");
  if(fp_hill == NULL){
    fprintf(stderr, "Cannot open file %s\n", p_param->hill_file);
    return EXIT_FAILURE;
  }


  /* Define the value of cmax based on the drug name */
  if( (int)round(p_param->cmax) == 0 ){
    cmax = get_cmax( p_param->drug_name );
  }
  else{
    cmax = p_param->cmax;
  }
  printf("Drug name: %s\t\tCmax: %lf", p_param->drug_name, cmax);

  /* Skip the header in the sample file */
  fgets(herg_sample, sizeof(herg_sample), fp_herg);
  /* Skip the header in the sample file */
  fgets(hill_sample, sizeof(hill_sample), fp_hill);


  /* begin population loop */
  sim_id = 1; 
  while(fgets(herg_sample, sizeof(herg_sample), fp_herg) != NULL)
  {
    i = 0;
    
    fgets(hill_sample, sizeof(hill_sample), fp_hill);
    pch = strtok(hill_sample, ",");
    /* begin splitting sample loop */
    do
    {
      herg[i++] = strtod(pch, NULL);
      pch = strtok(NULL, ",");
    }
    while( pch != NULL );
    /* end splitting sample loop */


    pch = strtok(herg_sample, ",");
    /* begin splitting sample loop */
    do
    {
      herg[i++] = strtod(pch, NULL);
      pch = strtok(NULL, ",");
    }
    while( pch != NULL );
    /* end splitting sample loop */

    strncpy(dosages_str, p_param->dosages, sizeof(dosages_str));
    pch = strtok(dosages_str, ",");
    /* begin splitting the dosages_str string loop */
    do{
      dosage = strtod(pch, NULL);
      conc = cmax * dosage;

      p_cell->initConsts(0., conc, herg, p_param->drug_name);
      if( p_param->is_using_output == 1 ){
        load_last_state(p_cell);
      }
      p_cell->isS1 = 1;
      p_cell->CONSTANTS[stim_period] = p_param->bcl_init;
      p_cell->CONSTANTS[stim_amplitude] *= p_param->stim_amp;

      cvode_mem = CVodeCreate( CV_BDF,CV_NEWTON );
      CVodeSetUserData( cvode_mem, p_cell );
      states_vec = N_VMake_Serial( p_cell->states_size, p_cell->STATES );
      CVodeInit( cvode_mem, rhs_fn, T0, states_vec );
      t_out = T1;
      CVodeSetMaxStep( cvode_mem, p_param->dt );
      CVDense( cvode_mem, p_cell->states_size );
      CVodeSStolerances( cvode_mem, 1.0e-7, 1.0e-7 );

      int iout = 0;
      int imax = ( p_param->num_pace1 *  p_cell->CONSTANTS[stim_period] ) / p_param->dt;
      print_freq = (1. / p_param->dt) * p_param->dt_write;
      sprintf(result_file, "%s_%.1lf_vmcheck_#%d.plt", p_param->drug_name, dosage, sim_id );
      result_file[sizeof(result_file) - 1] = '\0';
      fp_vm = fopen( result_file, "w" );
      sprintf(result_file, "%s_%.1lf_ca_i_#%d.plt", p_param->drug_name, dosage, sim_id );
      result_file[sizeof(result_file) - 1] = '\0';
      fp_ca = fopen( result_file, "w" );
      sprintf(result_file, "%s_%.1lf_ires_#%d.plt", p_param->drug_name, dosage, sim_id );
      result_file[sizeof(result_file) - 1] = '\0';
      fp_ires = fopen( result_file, "w" );
      sprintf(result_file, "%s_%.1lf_apd.plt", p_param->drug_name, dosage );
      result_file[sizeof(result_file) - 1] = '\0';
      fp_apd = fopen( result_file, "a" );
      sprintf(result_file, "%s_%.1lf_qnet.plt", p_param->drug_name, dosage );
      result_file[sizeof(result_file) - 1] = '\0';
      fp_qnet = fopen( result_file, "a" );
      sprintf(result_file, "%s_%.1lf_qinward.plt", p_param->drug_name, dosage );
      result_file[sizeof(result_file) - 1] = '\0';
      fp_qinward = fopen( result_file, "a" );

     
      is_apd_written = false;
      is_ap_increasing = false;
      is_action_potential = false;
      v_prev = 0.0;
      v_top = 0.0;
      v_valley = 0.0;
      v_apd90 = 0.0;

      fprintf( fp_vm, "%s %s\n", 
               "TIME", 
               "Vm" );
      fprintf( fp_ca, "%s %s\n", 
               "TIME", 
               "cai" );

      fprintf( fp_ires, "%s %s %s %s %s %s %s %s %s\n",
               "TIME",
               "INaL",
               "ICaL",
               "Ito",
               "IKr",
               "IKs",
               "IK1",
               "INaCa",
               "INaK" );

     
      /* begin simulation loop */
      while(1){
        v_prev = p_cell->STATES[V];

        retval = CVode( cvode_mem, t_out, states_vec, &t, CV_NORMAL  );
        if( retval == CV_SUCCESS ){
          iout++;
          t_out += p_param->dt;
        }
        if (iout >= imax) break;

        // find out whether the simulation can produce ACTION POTENTIAL or not
        if (p_cell->STATES[V] > v_prev && !is_ap_increasing) {
          is_ap_increasing = true;
          v_valley = v_prev;
          is_action_potential = v_top - v_valley > -40;
          if (is_action_potential) {
            //v_apd90 = v_top - (0.9 * (v_top - v_valley));
            v_apd90 = -70;
          }
        }
        else if (p_cell->STATES[V] < v_prev && is_ap_increasing) {
          is_ap_increasing = false;
          v_top = v_prev;
        }


        // repolarization, when the potential moves from HIGH to LOW
        if (v_prev > v_apd90 && p_cell->STATES[V] <= v_apd90) {
          t_repolarize = t;
          apd = t_repolarize - t_depolarize;
          if( t >= (p_param->num_pace1 * p_param->bcl_init) - p_param->bcl_init && is_apd_written == false){
            fprintf( fp_apd, "%d %lf\n", 
                     sim_id,
                     apd  );
           is_apd_written = true;
          }
        }
        // depolarization, when the potential moves from LOW to HIGH
        if (v_prev <= v_apd90 && p_cell->STATES[V] > v_apd90) {
          t_depolarize = t;
       }

       if(iout % print_freq == 0){
         if( p_param->is_print_vm == 1 ){
           fprintf( fp_vm, "%lf %lf\n",
                    t,
                    p_cell->STATES[V] );
         }
         fprintf( fp_ca, "%lf %lf\n",
                  t,
                  p_cell->STATES[Ca_i] );
         fprintf( fp_ires, "%lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
                  t,
                  p_cell->ALGEBRAIC[i_NaL],
                  p_cell->ALGEBRAIC[i_CaL],
                  p_cell->ALGEBRAIC[i_to],
                  p_cell->ALGEBRAIC[i_Kr],
                  p_cell->ALGEBRAIC[i_Ks],
                  p_cell->ALGEBRAIC[i_K1],
                  p_cell->ALGEBRAIC[i_NaCa_i],
                  p_cell->ALGEBRAIC[i_NaK] );

         save_last_state(p_cell, t);
       }

       // save qnet and auc(area under curve) for the last and before last cycles
       if( iout % (int)(p_cell->CONSTANTS[stim_period]/p_param->dt) == 0 ){
           if( iout / (int)(p_cell->CONSTANTS[stim_period]/p_param->dt) 
              == p_param->num_pace1-2 ){
              qnet_prev = p_cell->STATES[qnet];
              INaL_auc_prev = p_cell->STATES[INaL_AUC];
              ICaL_auc_prev = p_cell->STATES[ICaL_AUC];
              printf("TEST1: %d %lf %lf %lf\n", iout, qnet_prev, INaL_auc_prev, ICaL_auc_prev);
           }
           else if( iout / (int)(p_cell->CONSTANTS[stim_period]/p_param->dt) 
              == p_param->num_pace1-1 ){
              qnet_curr = p_cell->STATES[qnet];
              INaL_auc_curr = p_cell->STATES[INaL_AUC];
              ICaL_auc_curr = p_cell->STATES[ICaL_AUC];
              printf("TEST2: %d %lf %lf %lf\n", iout, qnet_curr, INaL_auc_curr, ICaL_auc_curr);
           }
       }

    }
    /* end simulation loop */

    /* saves qnet result to the file */
    fprintf( fp_qnet, "%lf\n", (qnet_curr - qnet_prev)/1000.0 );

    /* calculate qinward for dosages > 0. */
    if((int)round(dosage) == 0){
      INaL_auc_control = INaL_auc_curr - INaL_auc_prev;
      ICaL_auc_control = ICaL_auc_curr - ICaL_auc_prev;
    }
    else{
      
      INaL_auc_drug = INaL_auc_curr - INaL_auc_prev;
      ICaL_auc_drug = ICaL_auc_curr - ICaL_auc_prev;
      qinward = ( (INaL_auc_drug/INaL_auc_control) + (ICaL_auc_drug/ICaL_auc_control) ) * 0.5;
      fprintf( fp_qinward, "%lf\n", qinward );
    }


    fclose( fp_qinward );
    fclose( fp_qnet );
    fclose( fp_ires );
    if( p_param->is_print_vm == 1 ) fclose( fp_vm );
    fclose( fp_ca );
    fclose( fp_apd );

    N_VDestroy(states_vec);
    CVodeFree(&cvode_mem);

    pch = strtok(NULL, ","); 
    }
    while( pch != NULL );
    /* end splitting the dosages_str string loop */
    sim_id++;
  }
  /* end population loop */

  // generate EDISON output
#ifdef IS_EDISON == 1
  system("rm -rf result");
  system("mkdir result");
  system("mv *.plt result");
  system("mv output.dat result");
#endif

  delete p_cell;

  return 0;
}
#endif


#if (defined OHARA_RUDY2011)
static int population_hill_ANM( int argc, char **argv, param_t *p_param, bool is_dutta )
{
  /* CellML object pointer */
  CellML *p_cell;
  /* time related variables */
  int print_freq;
  int err_code;
  int i;
  char file_name[100];
  FILE *fp_apd;
  FILE *fp_vm;
  FILE *fp_ca;
  FILE *fp_ires;

  /* CiPA related variables */
  FILE *fp_hill;
  FILE *fp_qnet;
  FILE *fp_qinward;
  double hill[14];
  double dosage;
  double conc;
  double cmax;
  char *pch;
  char dosages_str[100];
  char hill_sample[128];
  char result_file[128];
  int sim_id;

  /* APD related variables */
  bool is_ap_increasing;
  bool is_action_potential;
  bool is_apd_written;
  double apd;
  double v_prev;
  double v_top;
  double v_valley;
  double v_apd90;
  double t_depolarize;
  double t_repolarize;
  double bcl_curr;
  double bcl_dec;
  double t_stim_start;
  double t_stim_end;
  double t_di;
  double t_prev;
  double anm;
  double *apd_vec;
  int depolarization_counter;
  int apd_vec_size;
  bool is_alternant;
  int alternant_counter;
  FILE *fp_anm;
  FILE *fp_apd90;
  FILE *fp_apdr;

  /* SUNDIALs variables */
  int retval;
  double t_out, t_curr;
  void *cvode_mem;
  N_Vector states_vec;



  /* Initialize p_cell based on the ionmodel value */
  p_cell = init_cell( );
  if( p_cell == NULL ){
    fprintf(stderr, "Problem when initializing cellmodel\n");
    return EXIT_FAILURE;
  }

  /* Initialize mutation object of p_cell if possible */
  p_cell->isMutated = 0;
  if( strcmp(p_param->variant, "ori") != 0 ){
    p_cell = init_mutation( p_param->variant, 0.5, p_cell);
    if( p_cell->mutation != NULL ){
      p_cell->isMutated = 1;
    }
  }


  /* Open the sample file */
  fp_hill = fopen(p_param->hill_file, "r");
  if(fp_hill == NULL){
    fprintf(stderr, "Cannot open file %s\n", p_param->hill_file);
    return EXIT_FAILURE;
  }

  /* Define the value of cmax based on the drug name */
  cmax = get_cmax( p_param->drug_name );

  /* Skip the header in the sample file */
  fgets(hill_sample, sizeof(hill_sample), fp_hill);

  /* begin population loop */
  sim_id = 1; 
  while(fgets(hill_sample, sizeof(hill_sample), fp_hill) != NULL)
  {
    i = 0;
    printf("Sample:\n%s\n", hill_sample);
    pch = strtok(hill_sample, ",");
    /* begin splitting sample loop */
    do
    {
      hill[i++] = strtod(pch, NULL);
      pch = strtok(NULL, ",");
    }
    while( pch != NULL );
    /* end splitting sample loop */

    strncpy(dosages_str, p_param->dosages, sizeof(dosages_str));
    pch = strtok(dosages_str, ",");
    /* begin splitting the dosages_str string loop */
    do{
      dosage = strtod(pch, NULL);
      conc = cmax * dosage;

      p_cell->initConsts(p_param->celltype, conc, hill, is_dutta, p_param->drug_name);
      p_cell->isS1 = 1;
      p_cell->CONSTANTS[stim_period] = p_param->bcl_init;
      p_cell->CONSTANTS[stim_amplitude] *= p_param->stim_amp;

      cvode_mem = CVodeCreate( CV_BDF,CV_NEWTON );
      CVodeSetUserData( cvode_mem, p_cell );
      states_vec = N_VMake_Serial( p_cell->states_size, p_cell->STATES );
      CVodeInit( cvode_mem, rhs_fn, T0, states_vec );
      t_out = T1;
      t_curr = T0;
      CVodeSetMaxStep( cvode_mem, p_param->dt );
      CVDense( cvode_mem, p_cell->states_size );
      CVodeSStolerances( cvode_mem, 1.0e-7, 1.0e-7 );

      int iout = 0;
      int imax = ( p_param->num_pace1 *  p_cell->CONSTANTS[stim_period] ) / p_param->dt;
      print_freq = (1. / p_param->dt) * p_param->dt_write;
      if( p_param->is_print_vm == 1 ){
        sprintf(result_file, "%s_%.1lf_vmcheck_#%d.plt", p_param->drug_name, dosage, sim_id );
        result_file[sizeof(result_file) - 1] = '\0';
        fp_vm = fopen( result_file, "w" );
        sprintf(result_file, "%s_%.1lf_ca_i_#%d.plt", p_param->drug_name, dosage, sim_id );
        result_file[sizeof(result_file) - 1] = '\0';
        fp_ca = fopen( result_file, "w" );
      }
      sprintf(result_file, "%s_%.1lf_ires_#%d.plt", p_param->drug_name, dosage, sim_id );
      result_file[sizeof(result_file) - 1] = '\0';
      fp_ires = fopen( result_file, "w" );
      sprintf(result_file, "%s_%.1lf_apd.plt", p_param->drug_name, dosage );
      result_file[sizeof(result_file) - 1] = '\0';
      fp_apd = fopen( result_file, "a" );
      sprintf(result_file, "%s_%.1lf_qnet.plt", p_param->drug_name, dosage );
      result_file[sizeof(result_file) - 1] = '\0';
      fp_qnet = fopen( result_file, "a" );
      sprintf(result_file, "%s_%.1lf_qinward.plt", p_param->drug_name, dosage );
      result_file[sizeof(result_file) - 1] = '\0';
      fp_qinward = fopen( result_file, "a" );
      sprintf(result_file, "%s_%.1lf_anm#%d.plt", p_param->drug_name, dosage, sim_id );
      result_file[sizeof(result_file) - 1] = '\0';
      fp_anm = fopen( result_file, "a" );
      fp_apdr = fopen( "apdr.plt", "w" );
      fp_apd90 = fopen( "debug.plt", "w" ); 

     
      is_apd_written = false;
      is_ap_increasing = false;
      is_action_potential = false;
      is_alternant = false;
      v_prev = 0.0;
      v_top = 0.0;
      v_valley = 0.0;
      v_apd90 = p_param->v_apd90;
      bcl_curr = p_cell->CONSTANTS[stim_period];
      bcl_dec = p_param->bcl_decrement;
      apd_vec = new double[p_param->num_pace1];
      apd_vec_size = 0;
      depolarization_counter = 0;
      alternant_counter = 0;

      if( p_param->is_print_vm == 1 ){
        fprintf( fp_vm, "%s %s\n", 
                 "TIME", 
                 "Vm" );
      }
      fprintf( fp_ca, "%s %s\n", 
               "TIME", 
               "cai" );

      fprintf( fp_ires, "%s %s %s %s %s %s %s %s %s\n",
               "TIME",
               "INaL",
               "ICaL",
               "Ito",
               "IKr",
               "IKs",
               "IK1",
               "INaCa",
               "INaK" );

     
 
      /* begin simulation loop */
      while(1){
        v_prev = p_cell->STATES[V];
        if (t_curr <= t_prev + bcl_curr * p_param->num_pace1) {
          if (t_curr > t_stim_end) {
            printf( "%s = %d %s = %lf %s = %lf %s = %lf\n",
                    "depolarization_counter", depolarization_counter,
                    "t_stim_start", t_stim_start,
                    "t_stim_end", t_stim_end,
                    "bcl", bcl_curr );
            t_stim_start += bcl_curr;
            t_stim_end = t_stim_start + p_cell->CONSTANTS[stim_duration];
          }
        }
        else {
          printf("%s: %lf %s: %d\n", "BCL", bcl_curr, "Size APD", apd_vec_size);
          if( apd_vec_size < 11 ) {
            anm = -1000.0;
            is_alternant = true;
          }
          else{
            anm = get_ANM(apd_vec, apd_vec_size);
            if(anm > 0.05) is_alternant = true;
            if(alternant_counter > 4) is_alternant = true;
          }
          fprintf(fp_anm, "BCL: %lf ANM: %lf\n", bcl_curr, anm);

          if( is_alternant == true ) break; 
          bcl_curr -= bcl_dec;
          if (bcl_curr >= 750.0) bcl_curr = 500.0;
          if (bcl_curr <= 300.0) bcl_dec = 10;
          p_cell->CONSTANTS[stim_period] = bcl_curr;
          t_prev = t_curr;
          apd_vec_size = 0;
          alternant_counter = 0;
          if (bcl_curr < p_param->bcl_end) break;
        }


        retval = CVode( cvode_mem, t_out, states_vec, &t_curr, CV_NORMAL  );
        if( retval == CV_SUCCESS ){
          iout++;
          t_out += p_param->dt;
        }

        if(iout % (int)(p_cell->CONSTANTS[stim_period]/p_param->dt) == 0 ){
          p_cell->STATES[qnet] = 0.;
        }


        // find out whether the simulation can produce ACTION POTENTIAL or not
        if (p_cell->STATES[V] > v_prev && !is_ap_increasing) {
          is_ap_increasing = true;
          v_valley = v_prev;
          is_action_potential = v_top - v_valley > -40;
          if (is_action_potential) {
            //v_apd90 = v_top - (0.9 * (v_top - v_valley));
            //v_apd90 = -70;
           fprintf( fp_apd90, "%lf %lf %lf %lf\n",
                 round(t_curr), v_top, v_prev, v_apd90 );
         }
        }
        else if (p_cell->STATES[V] < v_prev && is_ap_increasing) {
          is_ap_increasing = false;
          v_top = v_prev;
        }


        // repolarization, when the potential moves from HIGH to LOW
        if (v_prev > v_apd90 && p_cell->STATES[V] <= v_apd90) {
          t_repolarize = t_curr;
          apd = t_repolarize - t_depolarize;
          if(apd > bcl_curr) {
             printf("WARNING!! APD is larger than BCL!! Alternant detected!!\n");
             alternant_counter++;
          }
          if( t_curr >= (p_param->num_pace1 * p_param->bcl_init) - p_param->bcl_init && is_apd_written == false){
            fprintf( fp_apd, "%d %lf\n", 
                     sim_id,
                     apd  );
           is_apd_written = true;
          }
          apd_vec[apd_vec_size] = apd;
          printf("APD VEC %d: %lf\n", apd_vec_size,  apd_vec[apd_vec_size]);  
          apd_vec_size++;
        }
        // depolarization, when the potential moves from LOW to HIGH
        if (v_prev <= v_apd90 && p_cell->STATES[V] > v_apd90) {
          t_depolarize = t_curr;
          depolarization_counter++;
          if (depolarization_counter >= 2) {
            t_di = t_depolarize - t_repolarize;
            fprintf( fp_apdr, "%d %lf %lf %lf %lf\n",
                   depolarization_counter, bcl_curr, apd,
                   t_di, t_curr );
         }
       }

       if(iout % print_freq == 0){
         if( p_param->is_print_vm == 1 ){
           fprintf( fp_vm, "%lf %lf\n",
                    t_curr,
                    p_cell->STATES[V] );
         }
         fprintf( fp_ca, "%lf %lf\n",
                  t_curr,
                  p_cell->STATES[Ca_i] );
         fprintf( fp_ires, "%lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
                  t_curr,
                  p_cell->ALGEBRAIC[i_NaL],
                  p_cell->ALGEBRAIC[i_CaL],
                  p_cell->ALGEBRAIC[i_to],
                  p_cell->ALGEBRAIC[i_Kr],
                  p_cell->ALGEBRAIC[i_Ks],
                  p_cell->ALGEBRAIC[i_K1],
                  p_cell->ALGEBRAIC[i_NaCa_i],
                  p_cell->ALGEBRAIC[i_NaK] );
       }

    }
    /* end simulation loop */

    fprintf( fp_qnet, "%lf\n",
             p_cell->STATES[qnet] );

    fclose( fp_qinward );
    fclose( fp_qnet );
    fclose( fp_ires );
    if( p_param->is_print_vm == 1 ) fclose( fp_vm );
    fclose( fp_ca );
    fclose( fp_apd );
    fclose( fp_anm );
    fclose( fp_apd90 );
    fclose( fp_apdr );

    /* generate APD vs DiastolicInterval and APD vs BCL result files */
    create_apdr_output(p_param->num_pace1, depolarization_counter + 1);


    N_VDestroy(states_vec);
    CVodeFree(&cvode_mem);
    delete []apd_vec;

    printf("END IT ALL\n");
    pch = strtok(NULL, ","); 
    }
    while( pch != NULL );
    /* end splitting the dosages_str string loop */
    sim_id++;
  }
  // generate EDISON output
#ifdef IS_EDISON == 1
  system("rm -rf result");
  system("mkdir result");
  system("mv *.plt result");
  exit(0);
#endif

  delete p_cell;

  return 0;
}
#endif



#if (defined ORUDY_CIPA2017)
static int population_herg_ANM( int argc, char **argv, param_t *p_param )
{
  /* CellML object pointer */
  CellML *p_cell;
  /* time related variables */
  int print_freq;
  int err_code;
  int i;
  char file_name[100];
  FILE *fp_apd;
  FILE *fp_vm;
  FILE *fp_ca;
  FILE *fp_ires;

  /* CiPA related variables */
  FILE *fp_herg;
  FILE *fp_hill;
  FILE *fp_qnet;
  FILE *fp_qinward;
  double herg[20];
  double dosage;
  double conc;
  double cmax;
  char *pch;
  char dosages_str[100];
  char hill_sample[128];
  char herg_sample[128];
  char result_file[128];
  int sim_id;

  /* APD related variables */
  bool is_ap_increasing;
  bool is_action_potential;
  bool is_apd_written;
  double apd;
  double v_prev;
  double v_top;
  double v_valley;
  double v_apd90;
  double t_depolarize;
  double t_repolarize;
  double bcl_curr;
  double bcl_dec;
  double t_stim_start;
  double t_stim_end;
  double t_di;
  double t_prev;
  double anm;
  double *apd_vec;
  int depolarization_counter;
  int apd_vec_size;
  bool is_alternant;
  int alternant_counter;
  FILE *fp_anm;
  FILE *fp_apd90;
  FILE *fp_apdr;


  /* SUNDIALs variables */
  int retval;
  double t_out, t_curr;
  void *cvode_mem;
  N_Vector states_vec;

  /* Initialize p_cell based on the ionmodel value */
  p_cell = init_cell( );
  if( p_cell == NULL ){
    fprintf(stderr, "Problem when initializing cellmodel\n");
    return EXIT_FAILURE;
  }

  /* Initialize mutation object of p_cell if possible */
  p_cell->isMutated = 0;
  if( strcmp(p_param->variant, "ori") != 0 ){
    p_cell = init_mutation( p_param->variant, 0.5, p_cell);
    if( p_cell->mutation != NULL ){
      p_cell->isMutated = 1;
    }
  }


  /* Open the sample file */
  fp_herg = fopen(p_param->herg_file, "r");
  if(fp_herg == NULL){
    fprintf(stderr, "Cannot open file %s\n", p_param->herg_file);
    return EXIT_FAILURE;
  }
  /* Open the sample file */
  fp_hill = fopen(p_param->hill_file, "r");
  if(fp_hill == NULL){
    fprintf(stderr, "Cannot open file %s\n", p_param->hill_file);
    return EXIT_FAILURE;
  }


  /* Define the value of cmax based on the drug name */
  cmax = get_cmax( p_param->drug_name );
  printf("Drug name: %s\t\tCmax: %lf", p_param->drug_name, cmax);

  /* Skip the header in the sample file */
  fgets(herg_sample, sizeof(herg_sample), fp_herg);
  /* Skip the header in the sample file */
  fgets(hill_sample, sizeof(hill_sample), fp_hill);


  /* begin population loop */
  sim_id = 1; 
  while(fgets(herg_sample, sizeof(herg_sample), fp_herg) != NULL)
  {
    i = 0;
    
    fgets(hill_sample, sizeof(hill_sample), fp_hill);
    pch = strtok(hill_sample, ",");
    /* begin splitting sample loop */
    do
    {
      herg[i++] = strtod(pch, NULL);
      pch = strtok(NULL, ",");
    }
    while( pch != NULL );
    /* end splitting sample loop */


    pch = strtok(herg_sample, ",");
    /* begin splitting sample loop */
    do
    {
      herg[i++] = strtod(pch, NULL);
      pch = strtok(NULL, ",");
    }
    while( pch != NULL );
    /* end splitting sample loop */

    strncpy(dosages_str, p_param->dosages, sizeof(dosages_str));
    pch = strtok(dosages_str, ",");
    /* begin splitting the dosages_str string loop */
    do{
      dosage = strtod(pch, NULL);
      conc = cmax * dosage;

      p_cell->initConsts(0., conc, herg, p_param->drug_name);
      p_cell->isS1 = 1;
      p_cell->CONSTANTS[stim_period] = p_param->bcl_init;
      p_cell->CONSTANTS[stim_amplitude] *= p_param->stim_amp;

      cvode_mem = CVodeCreate( CV_BDF,CV_NEWTON );
      CVodeSetUserData( cvode_mem, p_cell );
      states_vec = N_VMake_Serial( p_cell->states_size, p_cell->STATES );
      CVodeInit( cvode_mem, rhs_fn, T0, states_vec );
      t_out = T1;
      CVodeSetMaxStep( cvode_mem, p_param->dt );
      CVDense( cvode_mem, p_cell->states_size );
      CVodeSStolerances( cvode_mem, 1.0e-7, 1.0e-7 );

      int iout = 0;
      int imax = ( p_param->num_pace1 *  p_cell->CONSTANTS[stim_period] ) / p_param->dt;
      print_freq = (1. / p_param->dt) * p_param->dt_write;
      if( p_param->is_print_vm == 1 ){
        sprintf(result_file, "%s_%.1lf_vmcheck_#%d.plt", p_param->drug_name, dosage, sim_id );
        result_file[sizeof(result_file) - 1] = '\0';
        fp_vm = fopen( result_file, "w" );
        sprintf(result_file, "%s_%.1lf_ca_i_#%d.plt", p_param->drug_name, dosage, sim_id );
        result_file[sizeof(result_file) - 1] = '\0';
        fp_ca = fopen( result_file, "w" );
      }
      sprintf(result_file, "%s_%.1lf_ires_#%d.plt", p_param->drug_name, dosage, sim_id );
      result_file[sizeof(result_file) - 1] = '\0';
      fp_ires = fopen( result_file, "w" );
      sprintf(result_file, "%s_%.1lf_apd.plt", p_param->drug_name, dosage );
      result_file[sizeof(result_file) - 1] = '\0';
      fp_apd = fopen( result_file, "a" );
      sprintf(result_file, "%s_%.1lf_qnet.plt", p_param->drug_name, dosage );
      result_file[sizeof(result_file) - 1] = '\0';
      fp_qnet = fopen( result_file, "a" );
      sprintf(result_file, "%s_%.1lf_qinward.plt", p_param->drug_name, dosage );
      result_file[sizeof(result_file) - 1] = '\0';
      fp_qinward = fopen( result_file, "a" );
      sprintf(result_file, "%s_%.1lf_anm#%d.plt", p_param->drug_name, dosage, sim_id );
      result_file[sizeof(result_file) - 1] = '\0';
      fp_anm = fopen( result_file, "a" );
      fp_apdr = fopen( "apdr.plt", "w" );
      fp_apd90 = fopen( "debug.plt", "w" );
     
      is_apd_written = false;
      is_ap_increasing = false;
      is_action_potential = false;
      v_prev = 0.0;
      v_top = 0.0;
      v_valley = 0.0;
      v_apd90 = p_param->v_apd90;
      bcl_curr = p_cell->CONSTANTS[stim_period];
      bcl_dec = p_param->bcl_decrement;
      apd_vec = new double[p_param->num_pace1];
      apd_vec_size = 0;
      alternant_counter = 0;
      depolarization_counter = 0;

      if( p_param->is_print_vm == 1 ){
        fprintf( fp_vm, "%s %s\n", 
                 "TIME", 
                 "Vm" );
      }
      fprintf( fp_ca, "%s %s\n", 
               "TIME", 
               "cai" );

      fprintf( fp_ires, "%s %s %s %s %s %s %s %s %s\n",
               "TIME",
               "INaL",
               "ICaL",
               "Ito",
               "IKr",
               "IKs",
               "IK1",
               "INaCa",
               "INaK" );

     
 
      /* begin simulation loop */
      while(1){
        v_prev = p_cell->STATES[V];
        if (t_curr <= t_prev + bcl_curr * p_param->num_pace1) {
          if (t_curr > t_stim_end) {
            printf( "%s = %d %s = %lf %s = %lf %s = %lf\n",
                    "depolarization_counter", depolarization_counter,
                    "t_stim_start", t_stim_start,
                    "t_stim_end", t_stim_end,
                    "bcl", bcl_curr );
            t_stim_start += bcl_curr;
            t_stim_end = t_stim_start + p_cell->CONSTANTS[stim_duration];
          }
        }
        else {
          printf("%s: %lf %s: %d\n", "BCL", bcl_curr, "Size APD", apd_vec_size);
          if( apd_vec_size < 11 ) {
            anm = -1000.0;
            is_alternant = true;
          }
          else{
            anm = get_ANM(apd_vec, apd_vec_size);
            if(anm > 0.05) is_alternant = true;
            if(alternant_counter > 4) is_alternant = true;
          }
          fprintf(fp_anm, "BCL: %lf ANM: %lf\n", bcl_curr, anm);

          if( is_alternant == true ) break; 
          bcl_curr -= bcl_dec;
          if (bcl_curr >= 750.0) bcl_curr = 500.0;
          if (bcl_curr <= 300.0) bcl_dec = 10;
          p_cell->CONSTANTS[stim_period] = bcl_curr;
          t_prev = t_curr;
          apd_vec_size = 0;
          alternant_counter = 0;
          if (bcl_curr < p_param->bcl_end) break;
        }


        retval = CVode( cvode_mem, t_out, states_vec, &t_curr, CV_NORMAL  );
        if( retval == CV_SUCCESS ){
          iout++;
          t_out += p_param->dt;
        }

        // find out whether the simulation can produce ACTION POTENTIAL or not
        if (p_cell->STATES[V] > v_prev && !is_ap_increasing) {
          is_ap_increasing = true;
          v_valley = v_prev;
          is_action_potential = v_top - v_valley > -40;
          if (is_action_potential) {
            //v_apd90 = v_top - (0.9 * (v_top - v_valley));
            //v_apd90 = -70;
            fprintf( fp_apd90, "%lf %lf %lf %lf\n",
                   round(t_curr), v_top, v_prev, v_apd90 );
          }
        }
        else if (p_cell->STATES[V] < v_prev && is_ap_increasing) {
          is_ap_increasing = false;
          v_top = v_prev;
        }


        // repolarization, when the potential moves from HIGH to LOW
        if (v_prev > v_apd90 && p_cell->STATES[V] <= v_apd90) {
          t_repolarize = t_curr;
          apd = t_repolarize - t_depolarize;
           if(apd > bcl_curr) {
             printf("WARNING!! APD is larger than BCL!! Alternant detected!!\n");
             alternant_counter++;
          }
          if( t_curr >= (p_param->num_pace1 * p_param->bcl_init) - p_param->bcl_init && is_apd_written == false){
            fprintf( fp_apd, "%d %lf\n", 
                     sim_id,
                     apd  );
           is_apd_written = true;
          }
          apd_vec[apd_vec_size] = apd;
          printf("APD VEC %d: %lf\n", apd_vec_size,  apd_vec[apd_vec_size]);  
          apd_vec_size++;
        }
        // depolarization, when the potential moves from LOW to HIGH
        if (v_prev <= v_apd90 && p_cell->STATES[V] > v_apd90) {
          t_depolarize = t_curr;
          depolarization_counter++;
          if (depolarization_counter >= 2) {
            t_di = t_depolarize - t_repolarize;
            fprintf( fp_apdr, "%d %lf %lf %lf %lf\n",
                   depolarization_counter, bcl_curr, apd,
                   t_di, t_curr );
         }
       }

       if(iout % print_freq == 0){
         if( p_param->is_print_vm == 1 ){
           fprintf( fp_vm, "%lf %lf\n",
                    t_curr,
                    p_cell->STATES[V] );
         }
         fprintf( fp_ca, "%lf %lf\n",
                  t_curr,
                  p_cell->STATES[Ca_i] );
         fprintf( fp_ires, "%lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
                  t_curr,
                  p_cell->ALGEBRAIC[i_NaL],
                  p_cell->ALGEBRAIC[i_CaL],
                  p_cell->ALGEBRAIC[i_to],
                  p_cell->ALGEBRAIC[i_Kr],
                  p_cell->ALGEBRAIC[i_Ks],
                  p_cell->ALGEBRAIC[i_K1],
                  p_cell->ALGEBRAIC[i_NaCa_i],
                  p_cell->ALGEBRAIC[i_NaK] );
       }
    }
    /* end simulation loop */

    fprintf( fp_qnet, "%lf\n",
             p_cell->STATES[qnet] );

    fclose( fp_qinward );
    fclose( fp_qnet );
    fclose( fp_ires );
    if( p_param->is_print_vm == 1 ) fclose( fp_vm );
    fclose( fp_ca );
    fclose( fp_apd );
    fclose( fp_anm );
    fclose( fp_apd90 );
    fclose( fp_apdr );

    /* generate APD vs DiastolicInterval and APD vs BCL result files */
    create_apdr_output(p_param->num_pace1, depolarization_counter + 1);

    N_VDestroy(states_vec);
    CVodeFree(&cvode_mem);
    delete []apd_vec;

    pch = strtok(NULL, ","); 
    }
    while( pch != NULL );
    /* end splitting the dosages_str string loop */
    sim_id++;
  }
  /* end population loop */

  // generate EDISON output
#ifdef IS_EDISON == 1
  system("rm -rf result");
  system("mkdir result");
  system("mv *.plt result");
#endif

  delete p_cell;

  return 0;
}
#endif



static int AP( int argc, char **argv, param_t *p_param )
{
  /* CellML object pointer */
  CellML *p_cell;
  /* time related variables */
  int print_freq;
  int err_code;
  FILE *fp_ca;
  FILE *fp_vm;
  FILE *fp_ires;

  /* SUNDIALs variables */
  int retval;
  double tnext, tcurr;
  void *cvode_main;
  N_Vector states_vec;

  /* Initialize p_cell based on the ionmodel value */
  p_cell = init_cell( );
  if( p_cell == NULL ){
    return EXIT_FAILURE;
  }

  /* Initialize mutation object of p_cell if possible */
  p_cell->isMutated = 0;
  if( strcmp(p_param->variant, "ori") != 0 ){
    p_cell = init_mutation( p_param->variant, 0.5, p_cell);
    if( p_cell->mutation != NULL ){
      p_cell->isMutated = 1;
    }
  }


  /*
     Invoke the states' initial condition of the cell model
     References:
     https://www.cellml.org/getting-started/tutorials/usinggeneratedc
  */
  p_cell->initConsts();
  p_cell->isS1 = 1;
  p_cell->CONSTANTS[stim_period] = p_param->bcl_init;
  p_cell->CONSTANTS[stim_amplitude] *= p_param->stim_amp;

  // Create CVODE solver
  cvode_main = CVodeCreate( CV_BDF,CV_NEWTON );
  // Give p_cell as CVode User Data
  CVodeSetUserData( cvode_main, p_cell );
  // Create the states vector based on the STATES array
  states_vec = N_VMake_Serial( p_cell->states_size, p_cell->STATES );
  // Initalize CVODE solver
  CVodeInit( cvode_main, rhs_fn, T0, states_vec );
  // Set up the future time (must be more than 0.0)
  tnext = T1;
  // Set up the maximum step size (dt)
  CVodeSetMaxStep( cvode_main, p_param->dt );
  // Set up the linear solver
  CVDense( cvode_main, p_cell->states_size );
  // Set up the tolerances
  CVodeSStolerances( cvode_main, 1.0e-7, 1.0e-7 );

  int iout = 0;
  int imax = ( p_param->num_pace1 * p_cell->CONSTANTS[stim_period] ) / p_param->dt;
  print_freq = (1. / p_param->dt) * p_param->dt_write;
  fp_ca = fopen( "ca_i_result.plt", "w" );
  fp_vm = fopen( "vmcheck.plt", "w" );
  fp_ires = fopen( "ires.plt", "w" );
  fprintf( fp_ca, "%s %s\n", "TIME", "Ca_i" );
  fprintf( fp_vm, "%s %s \n", 
           "TIME", 
           "Vm" );
#if defined OHARA_RUDY2011 || defined ORUDY_CIPA2017
  fprintf( fp_ires, "%s %s %s %s %s %s %s %s %s\n",
           "TIME",
           "INaL",
           "ICaL",
           "Ito",
           "IKr",
           "IKs",
           "IK1",
           "INaCa",
           "INaK" );
#endif
  // core simulation loop
  time_t begin = time(NULL);
  while(1){

    retval = CVode( cvode_main, tnext, states_vec, &tcurr, CV_NORMAL  );
    if( retval == CV_SUCCESS ){
      iout++;
      tnext += p_param->dt;
    }
   
    if(iout % print_freq == 0){
      fprintf( fp_vm, "%lf %lf\n", tcurr, p_cell->STATES[V]);
      fprintf( fp_ca, "%lf %lf\n",
               tcurr,
               p_cell->STATES[Ca_i] );
#if defined OHARA_RUDY2011 || defined ORUDY_CIPA2017
      fprintf( fp_ires, "%lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
               tcurr,
               p_cell->ALGEBRAIC[i_NaL],
               p_cell->ALGEBRAIC[i_CaL],
               p_cell->ALGEBRAIC[i_to],
               p_cell->ALGEBRAIC[i_Kr],
               p_cell->ALGEBRAIC[i_Ks],
               p_cell->ALGEBRAIC[i_K1],
               p_cell->ALGEBRAIC[i_NaCa_i],
               p_cell->ALGEBRAIC[i_NaK] );
#endif

    }

   if (iout >= imax) break;
  }

  time_t end = time(NULL);
  // end simulation loop
  printf("Time elapsed for simulation is: %ld sec", end - begin);

  fclose( fp_vm );
  fclose( fp_ca );
  fclose( fp_ires );

  // Memory Cleanup
  N_VDestroy(states_vec);
  CVodeFree(&cvode_main);
  delete p_cell;
  return 0;
}

static int gate_variation_AP( int argc, char **argv, param_t *p_param )
{
  /* CellML object pointer */
  CellML *p_cell;
  /* time related variables */
  double t_curr;
  int print_freq;
  int err_code;
  char file_name[20];
  char command[255];
  FILE *fp_ca;
  FILE *fp_vm;

  /* Initialize p_cell based on the ionmodel value */
  p_cell = init_cell( );
  if( p_cell == NULL ){
    return 20;
  }

  /* Initialize mutation object of p_cell if possible */
  p_cell->isMutated = 0;
  if( strcmp(p_param->variant, "ori") != 0 ){
    p_cell = init_mutation( p_param->variant, p_param->mutation_scaling, p_cell);
    if( p_cell->mutation != NULL ){
      p_cell->isMutated = 1;
    }
  }

  p_cell->isS1 = 1;
  p_cell->initConsts();
  p_cell->CONSTANTS[stim_period] = p_param->bcl_init;
  p_cell->CONSTANTS[stim_amplitude] *= p_param->stim_amp;
  p_cell->CONSTANTS[stim_duration] = p_param->stim_dur;
  t_curr = 0.0;
  print_freq = (1. / p_param->dt) * 2;


  /* core simulation loop */
  int iout;
  int imax = ( p_param->num_pace1 *  p_cell->CONSTANTS[stim_period] ) / p_param->dt;
  //generate result folder for EDISON
#ifdef IS_EDISON == 1
  system("rm -rf result");
  system("mkdir result");
#endif
  for(double mult = 0.1; mult <= 10; mult += 0.05  ){
    iout = 0;
    p_cell->initConsts();
    p_cell->CONSTANTS[stim_period] = p_param->bcl_init;
    p_cell->CONSTANTS[stim_amplitude] *= p_param->stim_amp;
    p_cell->CONSTANTS[stim_duration] = p_param->stim_dur;
    p_cell->CONSTANTS[p_param->scaled_gate] *= mult;
    
    sprintf(file_name, "vmcheck_%.2lf.plt", mult );
    fp_vm = fopen( file_name, "w" );
    fprintf( fp_vm, "%s %s\n", "TIME", "Vm" );
    while(1){
      t_curr = iout * p_param->dt;
      p_cell->computeRates(t_curr, 
                           p_cell->CONSTANTS, 
                           p_cell->RATES, 
                           p_cell->STATES,
                           p_cell->ALGEBRAIC);
      p_cell->solveAnalytical(p_param->dt);

     /* print the result */
    if (iout % print_freq == 0 && t_curr >= p_param->t_write_vtk){
        fprintf( fp_vm, "%lf %lf\n", t_curr, p_cell->STATES[V] );
    }
   iout++;
      if( iout >= imax ) break; 
    }

    p_cell->CONSTANTS[p_param->scaled_gate] /= mult;

    fclose( fp_vm );
    // generate EDISON output
#ifdef IS_EDISON == 1
    sprintf(command, "./convert-to-edison-format.sh vmcheck_%.2lf.plt %s > vmcheck_%.2lf_edison.plt", mult, "ms",mult);
    system(command);
    system("mv *.plt result");
#endif
 }
 /* end simulation loop */


  delete p_cell;
  free( p_param );
  return 0;
}



#ifdef TN2006ENDO || TN2006EPI || TN2006M
static int multi_population_AP(int argc, char **argv, param_t *p_param)
{
  /* CellML object pointer */
  CellML *p_cell;
  /* time related variables */
  double t_curr;
  double t_prev;
  double t_stim_start;
  double t_stim_end;
  double t_di;
  double t_depolarize;
  double t_repolarize;
  double bcl;
  double bcl_decrement; 
  double ap_duration;
  /* volatge related variables */
  double v_prev;
  double v_top;
  double v_valley;
  double v_apd90;
  /* supporting variables */
  bool is_ap_increasing;
  bool is_action_potential;
  int depolarization_counter;
  int print_freq;
  int err_code;
  FILE *fp_sim;
  FILE *fp_res;

  double apd3[3 + 1][3 + 1];
  double apd9[9 + 1][9 + 1];
  double apd27[27 + 1][27 + 1];
  double apd81[81 + 1][81 + 1];
  double apd243[243 + 1][243 + 1];

  double mean_ca_i3[3 + 1][3 + 1];
  double mean_ca_i9[9 + 1][9 + 1];
  double mean_ca_i27[27 + 1][27 + 1];
  double mean_ca_i81[81 + 1][81 + 1];
  double mean_ca_i243[243 + 1][243 + 1];

  double mean_i_nak3[3 + 1][3 + 1];
  double mean_i_nak9[9 + 1][9 + 1];
  double mean_i_nak27[27 + 1][27 + 1];
  double mean_i_nak81[81 + 1][81 + 1];
  double mean_i_nak243[243 + 1][243 + 1];


  double last_ca_i[10000];
  double last_i_nak[10000];
  double sum_ca_i;
  double sum_i_nak;
  double mean_ca_i;
  double mean_i_nak;
  
  /* Initialize p_cell based on the ionmodel value */
  p_cell = init_cell( );
  if( p_cell == NULL ){
    return 20;
  }

  char scenario_name[255] = "";
  char result_file[100] = "";
  char command[255] = "";
  int a;
  int b;
  int c;
  int d;
  int e;
  int f;
  int g;
  int h;
  int i;
  int j;
  int m;

#ifdef IS_EDISON == 1
  system("rm -rf result");
  system("mkdir result");
#endif
  for( a = 5; a <= 15; a += 5) {
    for( b = 5; b <= 15; b += 5) {
      for( c = 5; c <= 15; c += 5) {
        for( d = 5; d <= 15; d += 5) {
          for( e = 5; e <= 15; e += 5) {
            for( f = 5; f <= 15; f += 5) {
              for( g = 5; g <= 15; g += 5) {
                for( h = 5; h <= 15; h += 5) {
                  for( i = 5; i <= 15; i += 5) {
                    for( j = 5; j <= 15; j += 5) { 

                      p_cell->initConsts();

                      p_cell->CONSTANTS[g_Ks] *= (double)a/(double)10;
                      p_cell->CONSTANTS[g_Kr] *= (double)b/(double)10;      
                      p_cell->CONSTANTS[g_K1] *= (double)c/(double)10;
                      p_cell->CONSTANTS[g_Na] *= (double)d/(double)10;
                      p_cell->CONSTANTS[g_bna] *= (double)e/(double)10;
                      p_cell->CONSTANTS[g_CaL] *= (double)f/(double)10;
                      p_cell->CONSTANTS[g_bca] *= (double)g/(double)10;
                      p_cell->CONSTANTS[g_to] *= (double)h/(double)10;
                      p_cell->CONSTANTS[g_pCa] *= (double)i/(double)10;
                      p_cell->CONSTANTS[g_pK] *= (double)j/(double)10;
	
                      bcl = p_param->bcl_init;
                      bcl_decrement = p_param->bcl_decrement;
                      p_cell->CONSTANTS[stim_period] = p_param->bcl_init;
                      p_cell->CONSTANTS[stim_amplitude] *= p_param->stim_amp;
                      p_cell->CONSTANTS[stim_duration] = p_param->stim_dur;
                      is_ap_increasing = false;
                      is_action_potential = false;
                      t_curr = 0.0;
                      t_prev = 0.0;
                      t_stim_start = 0.0;
                      t_stim_end = t_stim_start + p_cell->CONSTANTS[stim_duration];
                      ap_duration = 0.0;
                      depolarization_counter = 0;
                      v_prev = 0.0;
                      v_top = 0.0;
                      v_valley = 0.0;
                      v_apd90 = 0.0;
                      t_di = 0.0;
                      print_freq = (1. / p_param->dt) * 2;
		
                      sum_ca_i = 0.0;
                      sum_i_nak = 0.0;	
                      m = 0;
					
                      sprintf(scenario_name, 
                      "g_Ks%.1fg_Kr%.1fg_K1%.1fg_Na%.1fg_bna%.1fg_CaL%.1fg_bca%.1fg_to%.1fg_pCa%.1fg_pK%.1f", 
                      (double)a/(double)10, 
                      (double)b/(double)10, 
                      (double)c/(double)10, 
                      (double)d/(double)10, 
                      (double)e/(double)10, 
                      (double)f/(double)10, 
                      (double)g/(double)10, 
                      (double)h/(double)10, 
                      (double)i/(double)10, 
                      (double)j/(double)10);
                      sprintf( result_file, "%s-%s.plt", "Result", scenario_name );	
                      fp_sim = fopen( result_file, "w" );
                      fprintf( fp_sim, "%s %s %s %s\n",
                           "TIME","Vm","Ca_i","i_NaK" );

                      /* START SIMULATION CORE */
                      for (int k = 0; ; k++) {
                        t_curr = k * p_param->dt;
                        v_prev = p_cell->STATES[0];
                        /* 
                          Control the flow of simulation.
                          In case of APDR simulation, this part will also 
                          decrease the bcl based on the bcl_decrement value.
                        */
                       if (t_curr <= t_prev + bcl * p_param->num_pace1) {
                         if (t_curr > t_stim_end) {
                           /*printf( "%s = %d %s = %lf %s = %lf %s = %lf\n",
                             "depolarization_counter", depolarization_counter,
                             "t_stim_start", t_stim_start,
                             "t_stim_end", t_stim_end,
                             "bcl", bcl );*/
                          t_stim_start += bcl;
                          t_stim_end = t_stim_start + p_cell->CONSTANTS[stim_duration];
                         }
                       }
                       else {
                         // since no BCL decrement, just break it when enter here
                         fclose( fp_sim );
                         // generate EDISON output
#ifdef IS_EDISON == 1
                         sprintf(command, 
                              "./convert-to-edison-format.sh Result-%s.plt %s > Result_%s_edison.plt", 
                              scenario_name, "ms", scenario_name);
                         system(command);
                         system("mv *.plt result");
#endif
                         break;
                       }
                       /* calculate ALGEBRAIC and RATES */
                       p_cell->computeRates( t_curr, 
                                            p_cell->CONSTANTS,
                                            p_cell->RATES,
                                            p_cell->STATES,
                                            p_cell->ALGEBRAIC );
                       /* solve the RATES using Analytical Method */
                       p_cell->solveAnalytical(p_param->dt);
 
                       /* find out whether the simulation can produce ACTION POTENTIAL or not */
                       if (p_cell->STATES[0] > v_prev && !is_ap_increasing) {
                         is_ap_increasing = true;
                         v_valley = v_prev;
                         is_action_potential = v_top - v_valley > p_param->v_apd90;
                         if (is_action_potential) {
                           //v_apd90 = v_top - (0.9 * (v_top - v_valley));
                           v_apd90 = -70;
                         }
                       }
                       else if (p_cell->STATES[0] < v_prev && is_ap_increasing) {
                         is_ap_increasing = false;
                         v_top = v_prev;
                       }
                       // repolarization, when the potential moves from HIGH to LOW
                       if (v_prev > v_apd90 && p_cell->STATES[0] <= v_apd90) {
                         t_repolarize = t_curr;
                         ap_duration = t_repolarize - t_depolarize;
                       }
                       // depolarization, when the potential moves from LOW to HIGH
                       if (v_prev <= v_apd90 && p_cell->STATES[0] > v_apd90) {
                         depolarization_counter++;
                         t_depolarize = t_curr;
                         if (depolarization_counter >= 2) {
                           t_di = t_depolarize - t_repolarize;
                         }
                       }

                       // retrieve some informations from the last action potential
                       if( t_curr >= (p_param->num_pace1 * bcl) - bcl ){
                         fprintf( fp_sim, "%lf %lf %lf %lf\n",
                           t_curr,
                           p_cell->STATES[0],
                           p_cell->STATES[Ca_i],
                           p_cell->ALGEBRAIC[i_NaK] );
                        
                       }
                       sum_ca_i += p_cell->STATES[Ca_i];
                       sum_i_nak += p_cell->ALGEBRAIC[i_NaK];
                       m++;
                     }
                      /* END SIMULATION CORE */

                      // calculate both of the average Ca_i and i_NaK
                      m--;
                      mean_ca_i = sum_ca_i / m;
                      mean_i_nak = sum_i_nak / m;

                      apd3[i/5][j/5] = ap_duration;
                      mean_ca_i3[i/5][j/5] = mean_ca_i;
                      mean_i_nak3[i/5][j/5] = mean_i_nak;

                      fp_res = fopen( "Result-AP-Mean.plt", "a" );
                      if( ftell(fp_res) == 0 ){
                        fprintf( fp_res, "%s %s %s %s\n",
                                 "scenario_name", "ap_duration", "mean_ca_i", "mean_i_nak" );
                      }
                      fprintf( fp_res, "%s %lf %lf %lf\n",
                               scenario_name, 
                               ap_duration, 
                               mean_ca_i, 
                               mean_i_nak );
                      fclose(fp_res);

                      p_cell->CONSTANTS[g_Ks] /= (double)a/(double)10;
                      p_cell->CONSTANTS[g_Kr] /= (double)b/(double)10;
                      p_cell->CONSTANTS[g_K1] /= (double)c/(double)10;
                      p_cell->CONSTANTS[g_Na] /= (double)d/(double)10;
                      p_cell->CONSTANTS[g_bna] /= (double)e/(double)10;
                      p_cell->CONSTANTS[g_CaL] /= (double)f/(double)10;
                      p_cell->CONSTANTS[g_bca] /= (double)g/(double)10;
                      p_cell->CONSTANTS[g_to] /= (double)h/(double)10;
                      p_cell->CONSTANTS[g_pCa] /= (double)i/(double)10;
                      p_cell->CONSTANTS[g_pK] /= (double)j/(double)10;

                    }
                  }
                  for(int x = 1; x < 3 + 1; x++){
                    for (int y = 1; y < 4; y++) {
                      apd9[(g/5 - 1) * 3 + x][(h/5 - 1) * 3 + y] = apd3[x][y];
                      mean_ca_i9[(g/5 - 1) * 3 + x][(h/5 - 1) * 3 + y] = mean_ca_i3[x][y];
                      mean_i_nak9[(g/5 - 1) * 3 + x][(h/5 - 1) * 3 + y] = mean_i_nak3[x][y];
                    }
                  }
                  printf("kimoo: mapd9\n");
                }
              }
              for(int x=1; x<9+1; x++){ 
                for(int y = 1; y < 9 + 1; y++) {
                  apd27[(e/5 - 1) * 9 + x][(f/5 - 1) * 9 + y] = apd9[x][y];
                  mean_ca_i27[(e/5 - 1) * 9 + x][(f/5 - 1) * 9 + y] = mean_ca_i9[x][y];
                  mean_i_nak27[(e/5 - 1) * 9 + x][(f/5 - 1) * 9 + y] = mean_i_nak9[x][y];
                }
              }
              printf("kimoo: mapd27\n");
            }
          }
          for(int x=1; x<27+1; x++){ 
            for (int y = 1; y < 27 + 1; y++) {
              apd81[(c/5 - 1) * 27 + x][(d/5 - 1) * 27 + y] = apd27[x][y];
              mean_ca_i81[(c/5 - 1) * 27 + x][(d/5 - 1) * 27 + y] = mean_ca_i27[x][y];
              mean_i_nak81[(c/5 - 1) * 27 + x][(d/5 - 1) * 27 + y] = mean_i_nak27[x][y];
            }  
          }
          printf("kimoo: mapd81\n");
        }
      }
      for(int x=1; x<81+1; x++){ 
        for (int y = 1; y < 81 + 1; y++) {
          apd243[(a/5 - 1) * 81 + x][(b/5 - 1) * 81 + y] = apd81[x][y];
          mean_ca_i243[(a/5 - 1) * 81 + x][(b/5 - 1) * 81 + y] = mean_ca_i81[x][y];
          mean_i_nak243[(a/5 - 1) * 81 + x][(b/5 - 1) * 81 + y] = mean_i_nak81[x][y];
        }
      }
      printf("kimoo: mapd243\n");
    }
  }  printf("kimoo: the end of the loop\n");


  printf("Create matrices\n");
  FILE *fp_matrix;
  fp_matrix = fopen( "apd_matrix.txt", "w" );
  for (int y = 1; y < 243 + 1; y++) {
    for (int x = 1; x < 243 + 1; x++) {
      fprintf( fp_matrix, "%lf\n", apd243[y][x] );
    }
  }
  fclose( fp_matrix );
  fp_matrix = fopen( "mean_ca_i_matrix.txt", "w" );
  for (int y = 1; y < 243 + 1; y++) {
    for (int x = 1; x < 243 + 1; x++) {
      fprintf( fp_matrix, "%lf\n", mean_ca_i243[y][x] );
    }
  }
  fclose( fp_matrix );
  fp_matrix = fopen( "mean_i_nak_matrix.txt", "w" );
  for (int y = 1; y < 243 + 1; y++) {
    for (int x = 1; x < 243 + 1; x++) {
      fprintf( fp_matrix, "%lf\n", mean_i_nak243[y][x] );
    }
  }
  fclose( fp_matrix );
  printf("Finished creating matrices\n");


  delete p_cell;
  free( p_param );

  return 0;
}
#endif

static int APDR( int argc, char **argv, param_t *p_param )
{
  /* CellML object pointer */
  CellML *p_cell;
  /* time related variables */
  double t_prev;
  double t_stim_start;
  double t_stim_end;
  double t_di;
  double t_depolarize;
  double t_repolarize;
  double bcl;
  double bcl_decrement; 
  double apd;
  /* volatge related variables */
  double v_prev;
  double v_top;
  double v_valley;
  double v_apd90;
  /* supporting variables */
  bool is_ap_increasing;
  bool is_action_potential;
  int depolarization_counter;
  int print_freq;
  int err_code;
  FILE *fp_apd90;
  FILE *fp_apdr;
  FILE *fp_vm;
  /* SUNDIALs variables */
  int retval;
  double t_out, t_curr;
  void *cvode_mem;
  N_Vector states_vec;
  
  /* Initialize p_cell based on the ionmodel value */
  p_cell = init_cell();
  if( p_cell == NULL ){
    return 20;
  }

  /* Initialize mutation object of p_cell if possible */
  p_cell->isMutated = 0;
  if( strcmp(p_param->variant, "ori") != 0 ){
    p_cell = init_mutation( p_param->variant, p_param->mutation_scaling, p_cell);
    if( p_cell->mutation != NULL ){
      p_cell->isMutated = 1;
    }
  }

  /* 
     Invoke the states' initial condition of the cell model 
     References:
     https://www.cellml.org/getting-started/tutorials/usinggeneratedc
  */
  p_cell->initConsts();

  cvode_mem = CVodeCreate( CV_BDF,CV_NEWTON );
  CVodeSetUserData( cvode_mem, p_cell );
  states_vec = N_VMake_Serial( p_cell->states_size, p_cell->STATES );
  CVodeInit( cvode_mem, rhs_fn, T0, states_vec );
  t_out = T1;
  CVodeSetMaxStep( cvode_mem, p_param->dt );
  CVDense( cvode_mem, p_cell->states_size );
  CVodeSStolerances( cvode_mem, 1.0e-7, 1.0e-7 );

  bcl = p_param->bcl_init;
  bcl_decrement = p_param->bcl_decrement;
  p_cell->CONSTANTS[stim_period] = p_param->bcl_init;
  p_cell->CONSTANTS[stim_amplitude] *= p_param->stim_amp;
  p_cell->CONSTANTS[stim_duration] = p_param->stim_dur;
  is_ap_increasing = false;
  is_action_potential = false;
  t_curr = 0.0;
  t_prev = 0.0;
  t_stim_start = 0.0;
  t_stim_end = t_stim_start + p_cell->CONSTANTS[stim_duration];
  apd = 0.0;
  depolarization_counter = 0;
  v_prev = 0.0;
  v_top = 0.0;
  v_valley = 0.0;
  v_apd90 = 0.0;
  t_di = 0.0;
  print_freq = (1. / p_param->dt) * 2;
 
  fp_apdr = fopen( "apdr.plt", "w" );
  fp_apd90 = fopen( "debug.plt", "w" );
  fp_vm = fopen( "vmcheck.plt", "w" );

  int iout = 0;

  /* Print the header of the output files */
  fprintf( fp_vm, "%s %s\n", "TIME", "Vm" );
  /* Extra flag for S1 stimulus */
  p_cell->isS1 = true;

  /* core simulation loop */
  for (int i = 0; ; i++) {
    v_prev = p_cell->STATES[0];
    /* 
      Control the flow of simulation.
      In case of APDR simulation, this part will also 
      decrease the bcl based on the bcl_decrement value.
    */
    if (bcl <= 300.0) bcl_decrement = 10;
    if (t_curr <= t_prev + bcl * p_param->num_pace1) {
      if (t_curr > t_stim_end) {
        printf( "%s = %d %s = %lf %s = %lf %s = %lf\n",
                "depolarization_counter", depolarization_counter,
                "t_stim_start", t_stim_start,
                "t_stim_end", t_stim_end,
                "bcl", bcl );
        t_stim_start += bcl;
        t_stim_end = t_stim_start + p_cell->CONSTANTS[stim_duration];
      }
    }
    else {
      bcl -= p_param->bcl_decrement;
      p_cell->CONSTANTS[stim_period] = bcl;
      t_prev = t_curr;
      if (bcl < p_param->bcl_end) break;
    }

    retval = CVode( cvode_mem, t_out, states_vec, &t_curr, CV_NORMAL  );
    if( retval == CV_SUCCESS ){
      iout++;
      t_out += p_param->dt;
    }
 
    /* find out whether the simulation can produce ACTION POTENTIAL or not */
    if (p_cell->STATES[0] > v_prev && !is_ap_increasing) {
      is_ap_increasing = true;
      v_valley = v_prev;
      is_action_potential = v_top - v_valley > p_param->v_apd90;
      if (is_action_potential) {
        //v_apd90 = v_top - (0.9 * (v_top - v_valley));
        fprintf( fp_apd90, "%lf %lf %lf %lf\n",
                 round(t_curr), v_top, v_prev, v_apd90 );
      }
    }
    else if (p_cell->STATES[0] < v_prev && is_ap_increasing) {
      is_ap_increasing = false;
      v_top = v_prev;
    }
    // repolarization, when the potential moves from HIGH to LOW
    if (v_prev > v_apd90 && p_cell->STATES[0] <= v_apd90) {
      t_repolarize = t_curr;
      apd = t_repolarize - t_depolarize;
    }
    // depolarization, when the potential moves from LOW to HIGH
    if (v_prev <= v_apd90 && p_cell->STATES[0] > v_apd90) {
      depolarization_counter++;
      t_depolarize = t_curr;
      if (depolarization_counter >= 2) {
        t_di = t_depolarize - t_repolarize;
        fprintf( fp_apdr, "%d %lf %lf %lf %lf\n",
                 depolarization_counter, bcl, apd,
                 t_di, t_curr );
     }
    }
    /* print the result */
    if (i % print_freq == 0 ){ 
      fprintf( fp_vm, "%lf %lf\n", t_curr, p_cell->STATES[V] );
    }
  }


  fclose( fp_vm );
  fclose( fp_apd90 );
  fclose( fp_apdr );

  N_VDestroy(states_vec);
  CVodeFree(&cvode_mem);
 
  /* generate APD vs DiastolicInterval and APD vs BCL result files */
  create_apdr_output(p_param->num_pace1, depolarization_counter + 1);

  // generate output
#ifdef IS_EDISON == 1
    char command[255];
    sprintf(command, "./convert-to-edison-format.sh vmcehck.plt %s > vmcheck_edison.plt", "ms");
    system(command);
    system("mv *.plt result");
#endif


  delete p_cell;
  free( p_param );
  return 0;
}

static int ANM( int argc, char **argv, param_t *p_param )
{
  /* CellML object pointer */
  CellML *p_cell;
  /* time related variables */
  double t_prev;
  double t_stim_start;
  double t_stim_end;
  double t_di;
  double t_depolarize;
  double t_repolarize;
  double bcl_curr;
  double bcl_dec; 
  double apd;
  double anm;
  bool is_alternant;
  int alternant_counter;
  double *apd_vec;
  int apd_vec_size;
  double apd_mean;
  /* volatge related variables */
  double v_prev;
  double v_top;
  double v_valley;
  double v_apd90;
  /* supporting variables */
  bool is_ap_increasing;
  bool is_action_potential;
  int depolarization_counter;
  int print_freq;
  int err_code;
  FILE *fp_apd90;
  FILE *fp_apdr;
  FILE *fp_anm;
  FILE *fp_vm;
  /* SUNDIALs variables */
  int retval;
  double t_out, t_curr;
  void *cvode_mem;
  N_Vector states_vec;
  
  /* Initialize p_cell based on the ionmodel value */
  p_cell = init_cell();
  if( p_cell == NULL ){
    return 20;
  }

  /* Initialize mutation object of p_cell if possible */
  p_cell->isMutated = 0;
  if( strcmp(p_param->variant, "ori") != 0 ){
    p_cell = init_mutation( p_param->variant, p_param->mutation_scaling, p_cell);
    if( p_cell->mutation != NULL ){
      p_cell->isMutated = 1;
    }
  }

  /* 
     Invoke the states' initial condition of the cell model 
     References:
     https://www.cellml.org/getting-started/tutorials/usinggeneratedc
  */
  p_cell->initConsts();

  cvode_mem = CVodeCreate( CV_BDF,CV_NEWTON );
  CVodeSetUserData( cvode_mem, p_cell );
  states_vec = N_VMake_Serial( p_cell->states_size, p_cell->STATES );
  CVodeInit( cvode_mem, rhs_fn, T0, states_vec );
  t_out = T1;
  CVodeSetMaxStep( cvode_mem, p_param->dt );
  CVDense( cvode_mem, p_cell->states_size );
  CVodeSStolerances( cvode_mem, 1.0e-7, 1.0e-7 );

  bcl_curr = p_param->bcl_init;
  bcl_dec = p_param->bcl_decrement;
  p_cell->CONSTANTS[stim_period] = p_param->bcl_init;
  p_cell->CONSTANTS[stim_amplitude] *= p_param->stim_amp;
  p_cell->CONSTANTS[stim_duration] = p_param->stim_dur;
  is_ap_increasing = false;
  is_action_potential = false;
  t_curr = 0.0;
  t_prev = 0.0;
  t_stim_start = 0.0;
  t_stim_end = t_stim_start + p_cell->CONSTANTS[stim_duration];
  apd = 0.0;
  depolarization_counter = 0;
  v_prev = 0.0;
  v_top = 0.0;
  v_valley = 0.0;
  v_apd90 = 0.0;
  t_di = 0.0;
  print_freq = (1. / p_param->dt) * 2;
  apd_vec = new double[p_param->num_pace1];
  apd_vec_size = 0;
  alternant_counter = 0;
  is_alternant = false;
 
  fp_apdr = fopen( "apdr.plt", "w" );
  fp_anm = fopen( "anm.plt", "w" );
  fp_apd90 = fopen( "debug.plt", "w" );
  fp_vm = fopen( "vmcheck.plt", "w" );

  int iout = 0;

  /* Print the header of the output files */
  fprintf( fp_vm, "%s %s\n", "TIME", "Vm" );
  /* Extra flag for S1 stimulus */
  p_cell->isS1 = true;

  /* core simulation loop */
  for (int i = 0; ; i++) {
    v_prev = p_cell->STATES[0];
    /* 
      Control the flow of simulation.
      In case of APDR simulation, this part will also 
      decrease the bcl based on the bcl_decrement value.
    */
    if (bcl_curr <= 300.0) bcl_dec = 10;
    if (t_curr <= t_prev + bcl_curr * p_param->num_pace1) {
      if (t_curr > t_stim_end) {
        printf( "%s = %d %s = %lf %s = %lf %s = %lf\n",
                "depolarization_counter", depolarization_counter,
                "t_stim_start", t_stim_start,
                "t_stim_end", t_stim_end,
                "bcl", bcl_curr );
        t_stim_start += bcl_curr;
        t_stim_end = t_stim_start + p_cell->CONSTANTS[stim_duration];
      }
    }
    else {
        printf("%s: %lf %s: %d\n", "BCL", bcl_curr, "Size APD", apd_vec_size);
        if( apd_vec_size < 11 ) {
          anm = -1000.0;
          is_alternant = true;
        }
        else{
          anm = get_ANM(apd_vec, apd_vec_size);
          if(anm > 0.05) is_alternant = true;
          if(alternant_counter > 4) is_alternant = true;
        }
        fprintf(fp_anm, "BCL: %lf ANM: %lf\n", bcl_curr, anm);

        if( is_alternant == true ) break; 
        bcl_curr -= bcl_dec;
        if (bcl_curr >= 750.0) bcl_curr = 500.0;
        if (bcl_curr <= 300.0) bcl_dec = 10;
        p_cell->CONSTANTS[stim_period] = bcl_curr;
        t_prev = t_curr;
        apd_vec_size = 0;
        alternant_counter = 0;
        if (bcl_curr < p_param->bcl_end) break;
    }

    retval = CVode( cvode_mem, t_out, states_vec, &t_curr, CV_NORMAL  );
    if( retval == CV_SUCCESS ){
      iout++;
      t_out += p_param->dt;
    }

    /* find out whether the simulation can produce ACTION POTENTIAL or not */
    if (p_cell->STATES[0] > v_prev && !is_ap_increasing) {
      is_ap_increasing = true;
      v_valley = v_prev;
      is_action_potential = v_top - v_valley > p_param->v_apd90;
      if (is_action_potential) {
        //v_apd90 = v_top - (0.9 * (v_top - v_valley));
        v_apd90 = -70;
        fprintf( fp_apd90, "%lf %lf %lf %lf\n",
                 round(t_curr), v_top, v_prev, v_apd90 );
      }
    }
    else if (p_cell->STATES[0] < v_prev && is_ap_increasing) {
      is_ap_increasing = false;
      v_top = v_prev;
    }
    // repolarization, when the potential moves from HIGH to LOW
    if (v_prev > v_apd90 && p_cell->STATES[0] <= v_apd90) {
      t_repolarize = t_curr;
      apd = t_repolarize - t_depolarize;
      if(apd > bcl_curr) {
         printf("WARNING!! APD is larger than BCL!! Alternant detected!!\n");
         alternant_counter++;
      }
      apd_vec[apd_vec_size] = apd;
      printf("APD VEC %d: %lf\n", apd_vec_size,  apd_vec[apd_vec_size]);  
      apd_vec_size++;
    }
    // depolarization, when the potential moves from LOW to HIGH
    if (v_prev <= v_apd90 && p_cell->STATES[0] > v_apd90) {
      depolarization_counter++;
      t_depolarize = t_curr;
      if (depolarization_counter >= 2) {
        t_di = t_depolarize - t_repolarize;
        fprintf( fp_apdr, "%d %lf %lf %lf %lf\n",
                 depolarization_counter, bcl_curr, apd,
                 t_di, t_curr );
     }
    }
    /* print the result */
    if (i % print_freq == 0 ){ 
      fprintf( fp_vm, "%lf %lf\n", t_curr, p_cell->STATES[V] );
    }

  }


  fclose( fp_vm );
  fclose( fp_anm );
  fclose( fp_apd90 );
  fclose( fp_apdr );

  N_VDestroy(states_vec);
  CVodeFree(&cvode_mem);
  delete []apd_vec;
 
  /* generate APD vs DiastolicInterval and APD vs BCL result files */
  create_apdr_output(p_param->num_pace1, depolarization_counter + 1);

  // generate output
#ifdef IS_EDISON == 1
    char command[255];
    sprintf(command, "./convert-to-edison-format.sh vmcehck.plt %s > vmcheck_edison.plt", "ms");
    system(command);
    system("mv *.plt result");
#endif


  delete p_cell;
  free( p_param );
  return 0;
}


static int rhs_fn( realtype t, N_Vector y, N_Vector ydot, void *user_data )
{
  CellML *data = (CellML*)user_data;
  data->computeRates( t,
                      data->CONSTANTS,
                      N_VGetArrayPointer_Serial(ydot),
                      N_VGetArrayPointer_Serial(y),
                      data->ALGEBRAIC );
  return 0;
}

static double average( double *apd_arr, int start_idx, int end_idx )
{
  double sum = 0.0;
  for (int i = start_idx; i <= end_idx; i++ ) {
    sum += apd_arr[i];
  }
  return sum / (end_idx - start_idx + 1);
}

static double get_ANM( double *apd_arr, int length )
{
  double alternant_magnitude;
  double change_magnitude;
  double mean_apd;
  int start;

  start = length - 10 - 1;
  change_magnitude = 0.0;

  printf("Start: %d length: %d\n", start, length);
  for ( int i = start; i < length-1; i++ ) {
    change_magnitude += fabs(apd_arr[i] - apd_arr[i+1]);
  }

  alternant_magnitude = change_magnitude / 10.0;
  mean_apd = average(apd_arr, start, length-1);
  printf("Chg_mgntd: %lf Alt_mgntd: %lf mean_apd: %lf\n", change_magnitude, alternant_magnitude, mean_apd);
  return alternant_magnitude / mean_apd;
}


static void create_apdr_output(int num_pace1, int number_of_lines)
{
  long int CP = 0;

  int* icnt = (int*)malloc(number_of_lines * sizeof(int));

  printf("%s: %d\n", "Number of Lines", number_of_lines);

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
    if (bcl_arr[i] != bcl_arr[i - 1]) {
      CP += bcl_arr[i - 2] * num_pace1;
      fprintf(fp_apdr_last_2_cycles, "%lf %lf %lf %lf %lf %ld\n", 
              bcl_arr[i - 2], ap_duration_arr[i - 2], 
              ap_duration_arr[i - 1], di_arr[i - 2], di_arr[i - 1], CP);
      fprintf(fp_apdr_di, "%lf %lf %lf\n", 
              di_arr[i - 2], ap_duration_arr[i - 2], ap_duration_arr[i - 1]);
      fprintf(fp_apdr_bcl, "%lf %lf %lf\n", 
              bcl_arr[i - 2], ap_duration_arr[i - 2], ap_duration_arr[i - 1]);
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
