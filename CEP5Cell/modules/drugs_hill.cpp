#include "drugs_hill.hpp"
#include "rhs.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <cvode/cvode.h>
#include <cvode/cvode_dense.h>
#include <nvector/nvector_serial.h>

#if defined ORUDY2011_STATIC
int drugs_hill( int argc, char **argv, param_t *p_param, bool is_dutta )
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
  int iout;
  int imax;
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
      printf("Dosage: %lf\tCmax: %lf\tConcentration: %lf\n", dosage, cmax, conc);

      p_cell->initConsts(p_param->celltype, conc, hill, is_dutta, p_param->drug_name);
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

      iout = 0;
      imax = ( p_param->num_pace1 *  p_cell->CONSTANTS[stim_period] ) / p_param->dt;
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
      qnet_prev = 0.0;
      qnet_curr = 0.0;
      INaL_auc_prev = 0.0;
      ICaL_auc_prev = 0.0;
      INaL_auc_curr = 0.0;
      ICaL_auc_curr = 0.0;

      if( p_param->is_print_vm == 1 ){
        fprintf( fp_vm, "%s %s\n", 
                 "TIME", 
                 "Vm" );
      
      fprintf( fp_ca, "%s %s\n", 
               "TIME", 
               "cai" );

      fprintf( fp_ires, "%s %s %s %s %s %s %s %s %s %s\n",
               "TIME",
               "INa",
               "INaL",
               "ICaL",
               "Ito",
               "IKr",
               "IKs",
               "IK1",
               "INaCa",
               "INaK" );
      }
      fprintf( fp_apd, "%s %s %s\n",
               "AP_Dur",
               "T_Repol",
               "T_Depol");

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
          if( t >= (p_param->num_pace1 * p_param->bcl_init) - p_param->bcl_init && 
              is_apd_written == false){
            fprintf( fp_apd, "%lf %lf %lf\n",apd, t_repolarize, t_depolarize  );
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
         
           fprintf( fp_ca, "%lf %lf\n",
                  t,
                  p_cell->STATES[Ca_i] );
           fprintf( fp_ires, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
                    t,
                    p_cell->ALGEBRAIC[i_Na],
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

       // save qnet and auc(area under curve) for the last and before last cycles
       if( iout == 
           ((int)(p_cell->CONSTANTS[stim_period]/p_param->dt) * (p_param->num_pace1-2)) +
            (int)p_cell->CONSTANTS[stim_period]  ){
            qnet_prev = p_cell->STATES[qnet];
            INaL_auc_prev = p_cell->STATES[INaL_AUC];
            ICaL_auc_prev = p_cell->STATES[ICaL_AUC];
            printf("TEST1: %d %lf %lf %lf %lf\n", iout, t, qnet_prev, INaL_auc_prev, ICaL_auc_prev);
       }
       else if( iout ==
                ((int)(p_cell->CONSTANTS[stim_period]/p_param->dt) * (p_param->num_pace1-1)) +
                 (int)p_cell->CONSTANTS[stim_period]  ){
            qnet_curr = p_cell->STATES[qnet];
            INaL_auc_curr = p_cell->STATES[INaL_AUC];
            ICaL_auc_curr = p_cell->STATES[ICaL_AUC];
            printf("TEST2: %d %lf %lf %lf %lf\n", iout, t, qnet_curr, INaL_auc_curr, ICaL_auc_curr);
       }

    }
    /* end simulation loop */
    
    /* saves qnet result to the file */
    fprintf( fp_qnet, "%lf\n", (qnet_curr - qnet_prev)/1000.0 );

    /* calculate qinward for dosages > 0. */
    if((int)ceil(conc) == 0){
      INaL_auc_control = INaL_auc_curr - INaL_auc_prev;
      ICaL_auc_control = ICaL_auc_curr - ICaL_auc_prev;
    }
    else{
      INaL_auc_drug = INaL_auc_curr - INaL_auc_prev;
      ICaL_auc_drug = ICaL_auc_curr - ICaL_auc_prev;
      qinward = ( (INaL_auc_drug/INaL_auc_control) + (ICaL_auc_drug/ICaL_auc_control) ) * 0.5;
      fprintf( fp_qinward, "%lf\n", qinward );
    }

    if( p_param->is_print_vm == 1 ){
      fclose( fp_ires );
      fclose( fp_vm );
      fclose( fp_ca );
    }
    fclose( fp_qinward );
    fclose( fp_qnet );
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
  system("rm -rf result");
  system("mkdir result");
  system("mv *.plt result");
  system("mv output.dat result");

  delete p_cell;
  return 0;
}
#endif
