#include "apdr_bench.hpp"
#include "commons.hpp"
#include "../cellmodels/tentusscher_noble_noble_panfilov_2004_b.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <cvode/cvode.h>
#include <cvode/cvode_dense.h>

void create_apdr_output(int num_pace1, int number_of_lines);


void apdr_bench( int argc, char **argv, param_t *p_param )
{
  // I/O variables
  int print_freq;
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
  patch_clamp *p_cell;

  // SUNDIALs variables
  int retval;
  double tnext, tcurr;
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
  double bcl;
  double v_prev;
  double v_top;
  double v_valley;
  double v_apd90;

  p_cell = new tentusscher_noble_noble_panfilov_2004_b();
  p_cell->initConsts();
  p_cell->CONSTANTS[stim_period] = p_param->bcl_init;

  // Create CVODE solver
  cvode_main = CVodeCreate( CV_BDF,CV_NEWTON );
  // Give p_cell as CVode User Data
  CVodeSetUserData( cvode_main, p_cell );
  // Create the states vector based on the STATES array
  states_vec = N_VMake_Serial( p_cell->states_size, p_cell->STATES );
  // Initalize CVODE solver
  CVodeInit( cvode_main, rhs_fn, 0.0, states_vec );
  // Set up the maximum step size (dt)
  CVodeSetMaxStep( cvode_main, p_param->dt );
  // Set up the linear solver type
  CVDense( cvode_main, p_cell->states_size );
  // Set up the numerical error tolerances
  CVodeSStolerances( cvode_main, 1.0e-7, 1.0e-7 );

  bcl = p_param->bcl_init;
  bcl_decrement = p_param->bcl_decrement;
  t_prev = 0.0;
  t_stim_start = 0.0;
  t_stim_end = t_stim_start + p_cell->CONSTANTS[stim_duration];
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
  int imax = ( p_param->num_pace1 * p_cell->CONSTANTS[stim_period] ) / p_param->dt;
  tnext = p_param->dt;
  tcurr = 0.0;
  print_freq = (1. / p_param->dt) * p_param->dt_write;
  time_t begin = time(NULL);
  printf("BCL: %lf\n", bcl);
  while( bcl >= p_param->bcl_end ){
    v_prev = p_cell->STATES[0];
    if (bcl <= 300.0) bcl_decrement = 10;
    if (tcurr <= t_prev + bcl * p_param->num_pace1) {
      if (tcurr > t_stim_end) {
/*        printf( "%s = %d %s = %lf %s = %lf %s = %lf\n",
                "depolarization_counter", depolarization_counter,
                "t_stim_start", t_stim_start,
                "t_stim_end", t_stim_end,
                "bcl", bcl );
*/
        t_stim_start += bcl;
        t_stim_end = t_stim_start + p_cell->CONSTANTS[stim_duration];
      }
    }
    else {
      //printf("BCL: %lf LAST APD: %lf\n", bcl, apd);
      bcl -= p_param->bcl_decrement;
      p_cell->CONSTANTS[stim_period] = bcl;

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

      t_prev = tcurr;
      depolarization_counter = 0;
    }

    /* solve ODE */
    retval = CVode( cvode_main, tnext, states_vec, &tcurr, CV_NORMAL  );
    if( retval == CV_SUCCESS ){
      iout++;
      tnext += p_param->dt;
    }
    else{
      printf("CVode calculation error at %lf msec!!!\n", tcurr);
      break;
    }

    /* find out whether the simulation can produce ACTION POTENTIAL or not */
    if (p_cell->STATES[0] > v_prev && !is_ap_increasing) {
      is_ap_increasing = true;
      v_valley = v_prev;
      is_action_potential = v_top - v_valley > v_apd90;
      if (is_action_potential && v_prev < 0) {
        v_apd90 = v_top - (0.9 * (v_top - v_valley));
        fprintf( fp_apd90, "%lf %lf %lf %lf %lf\n",
                 tcurr, v_top, v_valley, v_apd90, bcl );
      }
    }
    else if (p_cell->STATES[0] < v_prev && is_ap_increasing) {
      is_ap_increasing = false;
      v_top = v_prev;
    }


    // repolarization, when the potential moves from HIGH to LOW
    if ( v_prev > v_apd90 && p_cell->STATES[0] <= v_apd90 ) {
      t_repolarize = tcurr;
      apd = t_repolarize - t_depolarize;
      fprintf(fp_apdur, "%d %lf %lf %lf %lf\n", depolarization_counter, t_depolarize, apd, v_apd90, bcl);
    }
    // depolarization, when the potential moves from LOW to HIGH
    if ( v_prev <= v_apd90 && p_cell->STATES[0] > v_apd90 ) {
      depolarization_counter++;
      apdr_file_lines++;
      t_depolarize_prev = t_depolarize;
      t_depolarize = tcurr;
      if (depolarization_counter >= 3) {
        t_di_prev = t_di;
        t_di = t_depolarize - t_repolarize;
        if(depolarization_counter > 3 && t_di-t_di_prev > 100.){
          apd = tcurr - t_depolarize_prev;
        }
        fprintf( fp_apdr, "%d %lf %lf %lf %lf\n",
                 depolarization_counter, bcl, apd, t_di, tcurr);
      }
    }

    if(iout % print_freq == 0){
      fprintf(fp_vm, "%lf %lf\n", tcurr, p_cell->STATES[V]);
      fprintf(fp_vm_all, "%lf %lf\n", tcurr, p_cell->STATES[V]);
      fprintf(fp_curr, "%lf %lf %lf %lf %lf %lf %lf %lf\n",
             tcurr, p_cell->ALGEBRAIC[i_Na], p_cell->ALGEBRAIC[i_CaL], 
             p_cell->ALGEBRAIC[i_NaCa], p_cell->ALGEBRAIC[i_to], p_cell->ALGEBRAIC[i_Kr], 
             p_cell->ALGEBRAIC[i_Ks], p_cell->ALGEBRAIC[i_K1]);
      fprintf(fp_curr_all, "%lf %lf %lf %lf %lf %lf %lf %lf\n",
           tcurr, p_cell->ALGEBRAIC[i_Na], p_cell->ALGEBRAIC[i_CaL], 
           p_cell->ALGEBRAIC[i_NaCa], p_cell->ALGEBRAIC[i_to], p_cell->ALGEBRAIC[i_Kr], 
           p_cell->ALGEBRAIC[i_Ks], p_cell->ALGEBRAIC[i_K1]);
      fprintf(fp_conc, "%lf %lf %f\n", tcurr, p_cell->STATES[Na_i], p_cell->STATES[Ca_i]);
      fprintf(fp_conc_all, "%lf %lf %f\n", tcurr, p_cell->STATES[Na_i], p_cell->STATES[Ca_i]);
    }

  }

  time_t end = time(NULL);
  // end simulation loop
  printf("Time elapsed for simulation is: %ld sec\n", end - begin);

  /* generate APD vs DiastolicInterval and APD vs BCL result files */
  create_apdr_output(p_param->num_pace1, apdr_file_lines + 1);
 

  // Memory Cleanup
  N_VDestroy(states_vec);
  CVodeFree(&cvode_main);
  fclose(fp_vm_all);
  fclose(fp_curr_all);
  fclose(fp_conc_all);
  fclose( fp_apd90 );
  fclose( fp_apdr );
  fclose( fp_apdur );
  delete p_cell;
  free(p_param);

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
