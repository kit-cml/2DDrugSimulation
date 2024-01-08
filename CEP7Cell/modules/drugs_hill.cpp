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
  FILE *fp_ap_profile;
  FILE *fp_ca_profile;
  double hill[14];
  double conc;
  char *pch;
  char concs_str[100];
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

  /* Calcium related variable */
  double ca_prev;
  double ca_peak;
  double ca_dia;
  double ca_amp;
  double t1_cad90;
  double t2_cad90;
  double cad90;
  double t1_cad50;
  double t2_cad50;
  double cad50;
  bool is_ca90_inc;
  bool is_ca50_inc;
  bool is_cad90_found;
  bool is_cad50_found;


  /* AP related variables */
  bool is_ap_increasing;
  bool is_action_potential;
  bool is_apd_written;
  double apd90;
  double apd50;
  double max_dvmdt;
  double vm_peak;
  double vm_dia;
  double v_prev;
  double v_apd90;
  double v_apd50;
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

  /* Skip the header in the sample file */
  fgets(hill_sample, sizeof(hill_sample), fp_hill);

  /* begin population loop */
  sim_id = 1;
  time_t begin = time(NULL); 
  while(fgets(hill_sample, sizeof(hill_sample), fp_hill) != NULL)
  {
    printf("Sample #%d: %s\n", sim_id, hill_sample);
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

    strncpy(concs_str, p_param->concs, sizeof(concs_str));
    pch = strtok(concs_str, ",");
    /* begin splitting the concs_str string loop */
    do{
      conc = strtod(pch, NULL);
      //printf("Concentration: %lf\n", conc);

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
        sprintf(result_file, "%s_%.1lf_vmcheck_#%d.plt", p_param->drug_name, conc, sim_id );
        result_file[sizeof(result_file) - 1] = '\0';
        fp_vm = fopen( result_file, "w" );
/*        sprintf(result_file, "%s_%.1lf_ca_i_#%d.plt", p_param->drug_name, conc, sim_id );
        result_file[sizeof(result_file) - 1] = '\0';
        fp_ca = fopen( result_file, "w" );
        sprintf(result_file, "%s_%.1lf_ires_#%d.plt", p_param->drug_name, conc, sim_id );
        result_file[sizeof(result_file) - 1] = '\0';
        fp_ires = fopen( result_file, "w" );
*/
      }
/*
      sprintf(result_file, "%s_%.1lf_apd#%d.plt", p_param->drug_name, conc, sim_id );
      result_file[sizeof(result_file) - 1] = '\0';
      fp_apd = fopen( result_file, "w" );
*/
      sprintf(result_file, "%s_%.1lf_ap_profile.plt", p_param->drug_name, conc );
      result_file[sizeof(result_file) - 1] = '\0';
      fp_ap_profile = fopen( result_file, "a" );
      sprintf(result_file, "%s_%.1lf_ca_profile.plt", p_param->drug_name, conc );
      result_file[sizeof(result_file) - 1] = '\0';
      fp_ca_profile = fopen( result_file, "a" );
      sprintf(result_file, "%s_%.1lf_qnet.plt", p_param->drug_name, conc );
      result_file[sizeof(result_file) - 1] = '\0';
      fp_qnet = fopen( result_file, "a" );
      sprintf(result_file, "%s_%.1lf_qinward.plt", p_param->drug_name, conc );
      result_file[sizeof(result_file) - 1] = '\0';
      fp_qinward = fopen( result_file, "a" );
      is_ap_increasing = false;
      is_action_potential = false;
      is_apd_written = false;
      v_apd90 = p_param->v_apd90;
      v_apd50 = p_param->v_apd50;
      qnet_prev = 0.0;
      qnet_curr = 0.0;
      INaL_auc_prev = 0.0;
      ICaL_auc_prev = 0.0;
      INaL_auc_curr = 0.0;
      ICaL_auc_curr = 0.0;
      max_dvmdt = p_cell->STATES[V];
      vm_peak = p_cell->STATES[V];

      ca_peak = p_cell->STATES[Ca_i];
      t1_cad90 = 0.0;
      t2_cad90 = 0.0;
      cad90 = 0.0;
      t1_cad50 = 0.0;
      t2_cad50 = 0.0;
      cad50 = 0.0;
      is_ca90_inc = false;
      is_ca50_inc = false;
      is_cad90_found = false;
      is_cad50_found = false;

      if( p_param->is_print_vm == 1 ){
        fprintf( fp_vm, "%s %s\n", 
                 "TIME", 
                 "Vm" );
/*      
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
*/
      }

      if(sim_id == 1){
/*        fprintf( fp_apd, "%s %s %s\n",
                 "AP_Durartion",
                 "T_Depolarize",
                 "T_Repolarize");
*/
        fprintf( fp_ap_profile, "%s %s %s %s %s %s\n",
                 "Max_Dvm/Dt",
                 "Vm_Peak",
                 "Vm_Resting",
                 "APD90",
                 "APD50",
                 "APDTri");
        fprintf( fp_ca_profile, "%s %s %s %s %s\n",
                 "Ca_Peak",
                 "Ca_Diastole",
                 "CaD90",
                 "CaD50",
                 "Catri");
     }


      /* begin simulation loop */
      while(1){
        v_prev = p_cell->STATES[V];
        ca_prev = p_cell->STATES[Ca_i];

        retval = CVode( cvode_mem, t_out, states_vec, &t, CV_NORMAL  );

        if( retval == CV_SUCCESS ){
          iout++;
          t_out += p_param->dt;
        }
        if (iout >= imax) break;

        // repolarization, when the potential moves from HIGH to LOW
        if (v_prev > v_apd90 && p_cell->STATES[V] <= v_apd90) {
          t_repolarize = t;
          apd90 = t_repolarize - t_depolarize;
//          fprintf( fp_apd, "%lf %lf %lf\n",apd90, t_depolarize, t_repolarize  );
          is_apd_written = true;
        }
        // depolarization, when the potential moves from LOW to HIGH
        if (v_prev <= v_apd90 && p_cell->STATES[V] > v_apd90) {
          t_depolarize = t;
       }
/*
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
*/
       // save qnet and auc(area under curve) for the last and before last cycles
       if( iout == 
           ((int)(p_cell->CONSTANTS[stim_period]/p_param->dt) * (p_param->num_pace1-2)) +
            (int)p_cell->CONSTANTS[stim_period]  ){
          qnet_prev = p_cell->STATES[qnet];
          INaL_auc_prev = p_cell->STATES[INaL_AUC];
          ICaL_auc_prev = p_cell->STATES[ICaL_AUC];
          //printf("TEST1: %d %lf %lf %lf %lf\n", iout, t, qnet_prev, INaL_auc_prev, ICaL_auc_prev);
       }
       else if( iout ==
                ((int)(p_cell->CONSTANTS[stim_period]/p_param->dt) * (p_param->num_pace1-1)) +
                 (int)p_cell->CONSTANTS[stim_period]  ){
          qnet_curr = p_cell->STATES[qnet];
          INaL_auc_curr = p_cell->STATES[INaL_AUC];
          ICaL_auc_curr = p_cell->STATES[ICaL_AUC];
          //printf("TEST2: %d %lf %lf %lf %lf\n", iout, t, qnet_curr, INaL_auc_curr, ICaL_auc_curr);
       }

       
       if( iout == ((int)(p_cell->CONSTANTS[stim_period]/p_param->dt) * (p_param->num_pace1-2)) ){
          // find diastolic calcium value and time
          ca_dia = p_cell->STATES[Ca_i];
       }

       if( t > p_cell->CONSTANTS[stim_period] * (p_param->num_pace1-2) &&
           t < p_cell->CONSTANTS[stim_period] * (p_param->num_pace1-1) ){
         // find out the peak calcium and time
         if ( ca_peak < p_cell->STATES[Ca_i]) { 
           ca_peak = p_cell->STATES[Ca_i];
           ca_amp = ca_peak - ca_dia;
         }
       }


       if( t > p_cell->CONSTANTS[stim_period] * (p_param->num_pace1-1) ){
         // get the cad90
         if( p_cell->STATES[Ca_i] < (0.1*ca_amp)+ca_dia && is_ca90_inc && !is_cad90_found ){
           t2_cad90 = t;
           cad90 = t2_cad90 - t1_cad90;
           //printf("CAD90 T2: %lf T1:%lf\n", t2_cad90, t1_cad90);
           is_cad90_found = true;
         }
         else if( p_cell->STATES[Ca_i] >= (0.1*ca_amp)+ca_dia && !is_ca90_inc ){
           t1_cad90 = t;
           is_ca90_inc = true;
         }

         // get the cad50
         if( p_cell->STATES[Ca_i] < (0.5*ca_amp)+ca_dia && is_ca50_inc && !is_cad50_found ){
           t2_cad50 = t;
           cad50 = t2_cad50 - t1_cad50;
           //printf("CAD50 T2: %lf T1:%lf\n", t2_cad50, t1_cad50);
           is_cad50_found = true;
         }
         else if( p_cell->STATES[Ca_i] >= (0.5*ca_amp)+ca_dia && !is_ca50_inc ){
           t1_cad50 = t;
           is_ca50_inc = true;
         }

         if(iout % print_freq == 0){
           if( p_param->is_print_vm == 1 ){
             fprintf( fp_vm, "%lf %lf\n",
                      t,
                      p_cell->STATES[V] );
           }
         }

         // get the maximum dv/dt at the last cycle
         p_cell->computeRates(t_out, p_cell->CONSTANTS, p_cell->RATES, p_cell->STATES, p_cell->ALGEBRAIC);
         if( p_cell->RATES[V] > max_dvmdt) max_dvmdt = p_cell->RATES[V];

         // get the peak Vm at the last cycle
         if( p_cell->STATES[V] > vm_peak) vm_peak = p_cell->STATES[V];

         // repolarization, when the potential moves from HIGH to LOW
         if (v_prev > v_apd50 && p_cell->STATES[V] <= v_apd50) {
           t_repolarize = t;
           apd50 = t_repolarize - t_depolarize;
         }
         // depolarization, when the potential moves from LOW to HIGH
         if (v_prev <= v_apd50 && p_cell->STATES[V] > v_apd50) {
           t_depolarize = t;
         }

       }


    }
    /* end simulation loop */
    
    /* saves results to the file */
    fprintf( fp_qnet, "%lf\n", (qnet_curr - qnet_prev)/1000.0 );
    fprintf( fp_ap_profile, "%lf %lf %lf %lf %lf %lf\n",
             max_dvmdt,
             vm_peak,
             p_cell->STATES[V],
             apd90,
             apd50,
             apd90-apd50);
    fprintf( fp_ca_profile, "%lf %lf %lf %lf %lf\n",
             ca_peak,
             ca_dia,
             cad90,
             cad50,
             cad90-cad50 );

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
      fclose( fp_vm );
/*
      fclose( fp_ires );
      fclose( fp_ca );
*/
    }
    fclose( fp_ap_profile );
    fclose( fp_ca_profile );
    fclose( fp_qinward );
    fclose( fp_qnet );
//    fclose( fp_apd );


    N_VDestroy(states_vec);
    CVodeFree(&cvode_mem);

    pch = strtok(NULL, ","); 
    }
    while( pch != NULL );
    /* end splitting the concs_str string loop */
    sim_id++;
  }
  /* end population loop */
  time_t end = time(NULL);
  printf("Time elapsed for simulation is: %ld sec", end - begin);

  // generate EDISON output
  system("rm -rf result");
  system("mkdir result");
  system("mv *.plt result");
  system("mv output.dat result");

  delete p_cell;
  return 0;
}
#endif
