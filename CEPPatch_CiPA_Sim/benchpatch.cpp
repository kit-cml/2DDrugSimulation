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

double glob_vm = 0.0;


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
  /* CellML object pointer */
  patch_clamp *p_patch;

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
  bool is_ca_peak;


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

  param_t *p_param;
  p_param = (param_t*)malloc( sizeof( param_t ) );
  edison_assign_params_single(argc,argv,p_param);

  p_patch = new Ohara_Rudy_2011();

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
      hill[i] = strtod(pch, NULL);
      pch = strtok(NULL, ",");
      i++;
    }
    while( pch != NULL );
    /* end splitting sample loop */

    strncpy(concs_str, p_param->concs, sizeof(concs_str));
    pch = strtok(concs_str, ",");

    /* begin splitting the concs_str string loop */
    do
    {
      conc = strtod(pch, NULL);
      //printf("Concentration: %lf\n", conc);

      p_patch->initConsts(p_param->celltype, conc, hill, true);
      p_patch->CONSTANTS[stim_period] = p_param->bcl_init;
      p_patch->CONSTANTS[amp] *= p_param->stim_amp;

      cvode_mem = CVodeCreate( CV_BDF,CV_NEWTON );
      CVodeSetUserData( cvode_mem, p_patch );
      states_vec = N_VMake_Serial( p_patch->states_size, p_patch->STATES );
      CVodeInit( cvode_mem, rhs_fn, 0.0, states_vec );
      t_out = p_param->dt;
      CVodeSetMaxStep( cvode_mem, p_param->dt );
      CVDense( cvode_mem, p_patch->states_size );
      CVodeSStolerances( cvode_mem, 1.0e-7, 1.0e-7 );

      print_freq = (1. / p_param->dt) * p_param->dt_write;
      if( p_param->is_print_vm == 1 ){
        sprintf(result_file, "%s_%.1lf_vmcheck_#%d.plt", p_param->drug_name, conc, sim_id );
        result_file[sizeof(result_file) - 1] = '\0';
        fp_vm = fopen( result_file, "w" );
        sprintf(result_file, "%s_%.1lf_ca_i_#%d.plt", p_param->drug_name, conc, sim_id );
        result_file[sizeof(result_file) - 1] = '\0';
        fp_ca = fopen( result_file, "w" );
/*        sprintf(result_file, "%s_%.1lf_ires_#%d.plt", p_param->drug_name, conc, sim_id );
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
      is_ca_peak = false;
      v_apd90 = p_param->v_apd90;
      v_apd50 = p_param->v_apd50;
      qnet_prev = 0.0;
      qnet_curr = 0.0;
      INaL_auc_prev = 0.0;
      ICaL_auc_prev = 0.0;
      INaL_auc_curr = 0.0;
      ICaL_auc_curr = 0.0;
      max_dvmdt = p_patch->STATES[V];
      vm_peak = p_patch->STATES[V];

      ca_peak = 0.0;
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
      fprintf( fp_ca, "%s %s\n",
               "TIME",
               "cai" );
/*      fprintf( fp_ires, "%s %s %s %s %s %s %s %s %s %s\n",
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

     iout = 0;
     imax = ( p_param->num_pace1 *  p_patch->CONSTANTS[stim_period] ) / p_param->dt;
     /* begin simulation loop */
      while(iout < imax){
        v_prev = p_patch->STATES[V];
        ca_prev = p_patch->STATES[cai];

        retval = CVode( cvode_mem, t_out, states_vec, &t, CV_NORMAL  );

        if( retval == CV_SUCCESS ){
          iout++;
          t_out += p_param->dt;
        }
        else{
          printf("CVode error at sample #%d and concentration %.1lf\n", sim_id, conc);
          break;
        }

        // repolarization, when the potential moves from HIGH to LOW
        if (v_prev > v_apd90 && p_patch->STATES[V] <= v_apd90) {
          t_repolarize = t;
          apd90 = t_repolarize - t_depolarize;
//          fprintf( fp_apd, "%lf %lf %lf\n",apd90, t_depolarize, t_repolarize  );
          is_apd_written = true;
        }
        // depolarization, when the potential moves from LOW to HIGH
        if (v_prev <= v_apd90 && p_patch->STATES[V] > v_apd90) {
          t_depolarize = t;
       }
/*
       if(iout % print_freq == 0){
         if( p_param->is_print_vm == 1 ){
           fprintf( fp_vm, "%lf %lf\n",
                    t,
                    p_patch->STATES[V] );
           fprintf( fp_ca, "%lf %lf\n",
                  t,
                  p_patch->STATES[cai] );
           fprintf( fp_ires, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
                    t,
                    p_patch->ALGEBRAIC[i_Na],
                    p_patch->ALGEBRAIC[i_NaL],
                    p_patch->ALGEBRAIC[i_CaL],
                    p_patch->ALGEBRAIC[i_to],
                    p_patch->ALGEBRAIC[i_Kr],
                    p_patch->ALGEBRAIC[i_Ks],
                    p_patch->ALGEBRAIC[i_K1],
                    p_patch->ALGEBRAIC[i_Nacai],
                    p_patch->ALGEBRAIC[i_NaK] );
           }
       }
*/
       // save qnet and auc(area under curve) for the last and before last cycles
       if( iout ==
           ((int)(p_patch->CONSTANTS[stim_period]/p_param->dt) * (p_param->num_pace1-2)) +
            (int)p_patch->CONSTANTS[stim_period]  ){
          qnet_prev = p_patch->STATES[qnet];
          INaL_auc_prev = p_patch->STATES[INaL_AUC];
          ICaL_auc_prev = p_patch->STATES[ICaL_AUC];
          //printf("TEST1: %d %lf %lf %lf %lf\n", iout, t, qnet_prev, INaL_auc_prev, ICaL_auc_prev);
       }
       else if( iout ==
                ((int)(p_patch->CONSTANTS[stim_period]/p_param->dt) * (p_param->num_pace1-1)) +
                 (int)p_patch->CONSTANTS[stim_period]  ){
          qnet_curr = p_patch->STATES[qnet];
          INaL_auc_curr = p_patch->STATES[INaL_AUC];
          ICaL_auc_curr = p_patch->STATES[ICaL_AUC];
          //printf("TEST2: %d %lf %lf %lf %lf\n", iout, t, qnet_curr, INaL_auc_curr, ICaL_auc_curr);
       }


       if( iout == ((int)(p_patch->CONSTANTS[stim_period]/p_param->dt) * (p_param->num_pace1-2)) ){
          // find diastolic calcium value and time
          ca_dia = p_patch->STATES[cai];
       }

       if( t > p_patch->CONSTANTS[stim_period] * (p_param->num_pace1-2) &&
           t < p_patch->CONSTANTS[stim_period] * (p_param->num_pace1-1) ){
         // find out the peak calcium and time
         if( ca_peak < p_patch->STATES[cai] ){
           ca_peak = p_patch->STATES[cai];
           ca_amp = ca_peak - ca_dia;
         }
       }


       if( t > p_patch->CONSTANTS[stim_period] * (p_param->num_pace1-1) ){
         // get the cad90
         if( p_patch->STATES[cai] < (0.1*ca_amp)+ca_dia && is_ca90_inc && !is_cad90_found ){
           t2_cad90 = t;
           cad90 = t2_cad90 - t1_cad90;
           printf("CAD90 T2: %lf T1:%lf\n", t2_cad90, t1_cad90);
           is_cad90_found = true;
         }
         
         else if( p_patch->STATES[cai] >= (0.1*ca_amp)+ca_dia && !is_ca90_inc ){
           t1_cad90 = t;
           is_ca90_inc = true;
         }

         // get the cad50
         if( (p_patch->STATES[cai])*10E6 < (0.5*ca_amp*10E6)+(ca_dia*10E6) && is_ca50_inc && !is_cad50_found ){
           t2_cad50 = t;
           cad50 = t2_cad50 - t1_cad50;
           printf("CAD50 T2: %lf T1:%lf\n", t2_cad50, t1_cad50);
           is_cad50_found = true;
         }
         else if( (p_patch->STATES[cai])*10E6 >= (0.5*ca_amp*10E6)+(ca_dia*10E6) && !is_ca50_inc ){
           t1_cad50 = t;
           is_ca50_inc = true;
         }

         if(iout % print_freq == 0){
           if( p_param->is_print_vm == 1 ){
             fprintf( fp_vm, "%lf %lf\n",
                      t,
                      p_patch->STATES[V] );
             fprintf( fp_ca, "%lf %lf\n",
                      t,
                      p_patch->STATES[cai] );

           }
         }

         // get the maximum dv/dt at the last cycle
         p_patch->computeRates(t_out, p_patch->CONSTANTS, p_patch->RATES, p_patch->STATES, p_patch->ALGEBRAIC);
         if( p_patch->RATES[V] > max_dvmdt) max_dvmdt = p_patch->RATES[V];

         // get the peak Vm at the last cycle
         if( p_patch->STATES[V] > vm_peak) vm_peak = p_patch->STATES[V];

         // repolarization, when the potential moves from HIGH to LOW
         if (v_prev > v_apd50 && p_patch->STATES[V] <= v_apd50) {
           t_repolarize = t;
           apd50 = t_repolarize - t_depolarize;
         }
         // depolarization, when the potential moves from LOW to HIGH
         if (v_prev <= v_apd50 && p_patch->STATES[V] > v_apd50) {
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
             p_patch->STATES[V],
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
        fclose( fp_ca );
        //fclose( fp_ires );
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
    sim_id++;
  }



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
    else if (strcasecmp(key, "GKr_Scale") == 0) {
      p_param->gkr_scale =  strtod(value, NULL) ;
      printf( "%s -- %lf\n", "GKr_Scale", p_param->gkr_scale);
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

