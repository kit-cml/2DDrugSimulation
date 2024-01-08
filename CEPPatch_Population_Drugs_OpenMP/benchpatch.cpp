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

void population_drug(int argc, char **argv);
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
  population_drug(argc, argv);

  return 0;
}



void population_drug(int argc, char **argv)
{
  int a,b,c,d,e,f,g,h;
  int start_val, end_val, inc;
  double scalars[8];
  char scenario_name[255] = "";


  /* patch_clamp object pointer */
  patch_clamp *p_patch;
  /* time related variables */
  int print_freq;
  int err_code;
  char file_name[100];
  FILE *fp_apd;
  FILE *fp_vm;
  FILE *fp_ca;
  FILE *fp_ires;
  FILE *fp_log;

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
  double qnet_tot;
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
  double cadtri;
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
  double apdtri;
  double t_depolarize;
  double t_repolarize;

  /* SUNDIALs variables */
  int retval;
  double t_out, t;
  void *cvode_mem;
  N_Vector states_vec;
  int iout;
  int imax;

  /* Matrices variables */
  double qnet4[4 + 1][4 + 1];
  double qnet16[16 + 1][16 + 1];
  double qnet64[64 + 1][64 + 1];
  double qnet256[256 + 1][256 + 1];

  double qinward4[4 + 1][4 + 1];
  double qinward16[16 + 1][16 + 1];
  double qinward64[64 + 1][64 + 1];
  double qinward256[256 + 1][256 + 1];

  double max_dvmdt4[4 + 1][4 + 1];
  double max_dvmdt16[16 + 1][16 + 1];
  double max_dvmdt64[64 + 1][64 + 1];
  double max_dvmdt256[256 + 1][256 + 1];

  double vm_peak4[4 + 1][4 + 1];
  double vm_peak16[16 + 1][16 + 1];
  double vm_peak64[64 + 1][64 + 1];
  double vm_peak256[256 + 1][256 + 1];

  double vm_dia4[4 + 1][4 + 1];
  double vm_dia16[16 + 1][16 + 1];
  double vm_dia64[64 + 1][64 + 1];
  double vm_dia256[256 + 1][256 + 1];

  double v_apd90_4[4 + 1][4 + 1];
  double v_apd90_16[16 + 1][16 + 1];
  double v_apd90_64[64 + 1][64 + 1];
  double v_apd90_256[256 + 1][256 + 1];

  double v_apd50_4[4 + 1][4 + 1];
  double v_apd50_16[16 + 1][16 + 1];
  double v_apd50_64[64 + 1][64 + 1];
  double v_apd50_256[256 + 1][256 + 1];

  double apdtri4[4 + 1][4 + 1];
  double apdtri16[16 + 1][16 + 1];
  double apdtri64[64 + 1][64 + 1];
  double apdtri256[256 + 1][256 + 1];


  double ca_peak4[4 + 1][4 + 1];
  double ca_peak16[16 + 1][16 + 1];
  double ca_peak64[64 + 1][64 + 1];
  double ca_peak256[256 + 1][256 + 1];

  double ca_dia4[4 + 1][4 + 1];
  double ca_dia16[16 + 1][16 + 1];
  double ca_dia64[64 + 1][64 + 1];
  double ca_dia256[256 + 1][256 + 1];

  double cad90_4[4 + 1][4 + 1];
  double cad90_16[16 + 1][16 + 1];
  double cad90_64[64 + 1][64 + 1];
  double cad90_256[256 + 1][256 + 1];

  double cad50_4[4 + 1][4 + 1];
  double cad50_16[16 + 1][16 + 1];
  double cad50_64[64 + 1][64 + 1];
  double cad50_256[256 + 1][256 + 1];

  double cadtri4[4 + 1][4 + 1];
  double cadtri16[16 + 1][16 + 1];
  double cadtri64[64 + 1][64 + 1];
  double cadtri256[256 + 1][256 + 1];

  //--------------------------------------------------------------------------------------------------------

  start_val = 100;
  end_val = 100;
  inc = -25;

  param_t *p_param;
  p_param = (param_t*)malloc( sizeof( param_t ) );
  edison_assign_params_single(argc,argv,p_param);

  // PLEASE DON'T PUT THIS CALL INSIDE THE DEEPEST LOOP!!!
  p_patch = new Ohara_Rudy_2011();

  a = p_param->gna_scale;
  b = p_param->gnal_scale;
  c = p_param->gcal_scale;
  d = p_param->gkr_scale;
          for(e = 100; e >= 25; e += inc){
            for(f = 100; f >= 25; f += inc){
              for(g = 100; g >= 25; g += inc){
                for(h = 100; h >= 25; h += inc){

                  scalars[0] = (double)a/(double)100;
                  scalars[1] = (double)b/(double)100;
                  scalars[2] = (double)c/(double)100;
                  scalars[3] = (double)d/(double)100;
                  scalars[4] = (double)e/(double)100;
                  scalars[5] = (double)f/(double)100;
                  scalars[6] = (double)g/(double)100;
                  scalars[7] = (double)h/(double)100;

                  p_patch->initConsts();
                  p_patch->CONSTANTS[stim_period] = p_param->bcl_init;
                  p_patch->CONSTANTS[amp] *= p_param->stim_amp;
                  p_patch->CONSTANTS[duration] = p_param->stim_dur;

                  p_patch->CONSTANTS[GNa] *= scalars[0];
                  p_patch->CONSTANTS[GNaL] *= scalars[1];
                  p_patch->CONSTANTS[PCa] *= scalars[2];
                  p_patch->CONSTANTS[GKr] *= scalars[3];
                  p_patch->CONSTANTS[GKs] *= scalars[4];
                  p_patch->CONSTANTS[GK1] *= scalars[5];
                  p_patch->CONSTANTS[GKb] *= scalars[6];
                  p_patch->CONSTANTS[Gto] *= scalars[7];


                  sprintf( scenario_name,
                    "GNa%.2fGNaL%.2fPCa%.2fGKr%.2fGKs%.2fGK1%.2fGKb%.2fGto%.2f",
                    scalars[0],
                    scalars[1],
                    scalars[2],
                    scalars[3],
                    scalars[4],
                    scalars[5],
                    scalars[6],
                    scalars[7] );
                  sprintf(result_file, "vmcheck_%s.plt", scenario_name );
                  result_file[sizeof(result_file) - 1] = '\0';
                  fp_vm = fopen( result_file, "w" );
                  sprintf(result_file, "ap_profile_%s.plt", scenario_name );
                  result_file[sizeof(result_file) - 1] = '\0';
                  fp_ap_profile = fopen( result_file, "a" );
                  sprintf(result_file, "ca_profile_%s.plt", scenario_name );
                  result_file[sizeof(result_file) - 1] = '\0';
                  fp_ca_profile = fopen( result_file, "a" );
                  sprintf(result_file, "qnet_%s.plt", scenario_name );
                  result_file[sizeof(result_file) - 1] = '\0';
                  fp_qnet = fopen( result_file, "a" );

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
                  qnet_tot = 0.0;
                  qinward = 0.0;


                  //printf("LOOP %s\n", scenario_name);
                  t_out = 0.01;
                  cvode_mem = CVodeCreate( CV_BDF,CV_NEWTON );
                  CVodeSetUserData( cvode_mem, p_patch );
                  states_vec = N_VMake_Serial( p_patch->states_size, p_patch->STATES );
                  CVodeInit( cvode_mem, rhs_fn, 0.0, states_vec );
                  CVodeSetMaxStep( cvode_mem, p_param->dt );
                  CVDense( cvode_mem, p_patch->states_size );
                  CVodeSStolerances( cvode_mem, 1.0e-7, 1.0e-7 );

                  iout = 0;
                  imax = ( p_param->num_pace1 *  p_patch->CONSTANTS[stim_period] ) / p_param->dt;
                  print_freq = (1. / p_param->dt) * p_param->dt_write;
                  fprintf( fp_vm, "%s %s\n", "TIME", "Vm" );
                  fprintf( fp_ap_profile, "%s %s %s %s %s %s\n","Max_Dvm/Dt","Vm_Peak","Vm_Resting","APD90","APD50","APDTri");
                  fprintf( fp_ca_profile, "%s %s %s %s %s\n","Ca_Peak","Ca_Diastole","CaD90","CaD50","Catri");


                  /* begin simulation loop */
                  while(iout < imax){
                    v_prev = p_patch->STATES[V];
                    ca_prev = p_patch->STATES[cai];
                    retval = CVode( cvode_mem, t_out, states_vec, &t, CV_NORMAL  );
                    if( retval == CV_SUCCESS ){
                      t_out += p_param->dt;
                    }
                    else{
                      printf("CVode error at scneario %s\n", scenario_name);
                      break;
                    }


                    // repolarization, when the potential moves from HIGH to LOW
                    if (v_prev > v_apd90 && p_patch->STATES[V] <= v_apd90) {
                      t_repolarize = t;
                      apd90 = t_repolarize - t_depolarize;
                     // fprintf( fp_apd, "%lf %lf %lf\n",apd90, t_depolarize, t_repolarize  );
                      is_apd_written = true;
                    }
                    // depolarization, when the potential moves from LOW to HIGH
                    if (v_prev <= v_apd90 && p_patch->STATES[V] > v_apd90) {
                     t_depolarize = t;
                    }

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
                      if ( ca_peak < p_patch->STATES[cai]) { 
                        ca_peak = p_patch->STATES[cai];
                        ca_amp = ca_peak - ca_dia;
                      }
                    }


                    if( t > p_patch->CONSTANTS[stim_period] * (p_param->num_pace1-1) ){
                      // get the cad90
                      if( p_patch->STATES[cai] < (0.1*ca_amp)+ca_dia && is_ca90_inc && !is_cad90_found ){
                        t2_cad90 = t;
                        cad90 = t2_cad90 - t1_cad90;
                        //printf("CAD90 T2: %lf T1:%lf CAD90:%lf\n", t2_cad90, t1_cad90, t2_cad90-t1_cad90);
                        is_cad90_found = true;
                      }
                      else if( p_patch->STATES[cai] >= (0.1*ca_amp)+ca_dia && !is_ca90_inc ){
                        t1_cad90 = t;
                        is_ca90_inc = true;
                      }
            
                      // get the cad50
                      if( p_patch->STATES[cai] < (0.5*ca_amp)+ca_dia && is_ca50_inc && !is_cad50_found ){
                        t2_cad50 = t;
                        cad50 = t2_cad50 - t1_cad50;
                        //printf("CAD50 T2: %lf T1:%lf CAD50: %lf\n", t2_cad50, t1_cad50, t2_cad50-t1_cad50);
                        is_cad50_found = true;
                      }
                      else if( p_patch->STATES[cai] >= (0.5*ca_amp)+ca_dia && !is_ca50_inc ){
                        t1_cad50 = t;
                        is_ca50_inc = true;
                      }
             
                      if(iout % print_freq == 0)fprintf( fp_vm, "%lf %lf\n",t,p_patch->STATES[V] );
           

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
                    iout++;
                  }
                  /* end simulation loop */
                  /* saves results to the file */
                  qnet_tot = (qnet_curr - qnet_prev)/1000.0;
                  cadtri = cad90-cad50;
                  apdtri = apd90-apd50;
                  vm_dia = p_patch->STATES[V];
                  fprintf( fp_qnet, "%lf\n", qnet_tot );
                  fprintf( fp_ap_profile, "%lf %lf %lf %lf %lf %lf\n",
                     max_dvmdt,vm_peak,vm_dia,apd90,apd50,apd90-apd50);
                  fprintf( fp_ca_profile, "%lf %lf %lf %lf %lf\n",
                     ca_peak,ca_dia,cad90,cad50,cad90-cad50 );

                  fp_log = fopen( "temp.plt", "a" );
                  if( ftell(fp_log) == 0 ){
                    fprintf( fp_log, "%s %s %s %s %s %s %s %s %s %s %s %s %s\n","scenario_name", "qnet", "max_dvmdt","vm_peak","vm_dia","v_apd90","v_apd50","apdtri","ca_peak","ca_dia","cad90","cad50","cadtri" );
                  }
                  fprintf( fp_log, "%s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",scenario_name, qnet_tot, max_dvmdt,vm_peak,vm_dia,apd90,apd50,apdtri,ca_peak,ca_dia,cad90,cad50,cadtri );


                  N_VDestroy(states_vec);
                  CVodeFree(&cvode_mem);
                  fclose( fp_ap_profile );
                  fclose( fp_ca_profile );
                  fclose( fp_qnet );
                  fclose( fp_vm );
                  fclose( fp_log );

                  qnet4[g/25][h/25] = qnet_tot;
                  qinward4[g/25][h/25] = qinward;
                  max_dvmdt4[g/25][h/25] = max_dvmdt;
                  vm_peak4[g/25][h/25] = vm_peak;
                  vm_dia4[g/25][h/25] = vm_dia;
                  v_apd90_4[g/25][h/25] = v_apd90;
                  v_apd50_4[g/25][h/25] = v_apd50;
                  apdtri4[g/25][h/25] = apdtri;
                  ca_peak4[g/25][h/25] = ca_peak;
                  ca_dia4[g/25][h/25] = ca_dia;
                  cad90_4[g/25][h/25] = cad90;
                  cad50_4[g/25][h/25] = cad50;
                  cadtri4[g/25][h/25] = cadtri;

                  printf("kimoo: qnet4\n");
                }
              }

              for(int x = 1; x < 4+1; x++){
                for (int y = 1; y < 4+1; y++) {
                  qnet16[(4 - (e/25))*4 + x][(4 - (f/25))*4 + y] = qnet4[x][y];
                  qinward16[(4 - (e/25))*4 + x][(4 - (f/25))*4 + y] = qinward4[x][y];
                  max_dvmdt16[(4 - (e/25))*4 + x][(4 - (f/25))*4 + y] = max_dvmdt4[x][y];
                  vm_peak16[(4 - (e/25))*4 + x][(4 - (f/25))*4 + y] = vm_peak4[x][y];
                  vm_dia16[(4 - (e/25))*4 + x][(4 - (f/25))*4 + y] = vm_dia4[x][y];
                  v_apd90_16[(4 - (e/25))*4 + x][(4 - (f/25))*4 + y] = v_apd90_4[x][y];
                  v_apd50_16[(4 - (e/25))*4 + x][(4 - (f/25))*4 + y] = v_apd50_4[x][y];
                  apdtri16[(4 - (e/25))*4 + x][(4 - (f/25))*4 + y] = apdtri4[x][y];
                  ca_peak16[(4 - (e/25))*4 + x][(4 - (f/25))*4 + y] = ca_peak4[x][y];
                  ca_dia16[(4 - (e/25))*4 + x][(4 - (f/25))*4 + y] = ca_dia4[x][y];
                  cad90_16[(4 - (e/25))*4 + x][(4 - (f/25))*4 + y] = cad90_4[x][y];
                  cad50_16[(4 - (e/25))*4 + x][(4 - (f/25))*4 + y] = cad50_4[x][y];
                  cadtri16[(4 - (e/25))*4 + x][(4 - (f/25))*4 + y] = cadtri4[x][y];
                }
              }
              printf("kimoo: qnet16\n");

            }
          }

          for(int x = 1; x < 16+1; x++){
            for (int y = 1; y < 16+1; y++) {
              qnet64[(4 - (c/25))*4 + x][(4 - (d/25))*4 + y] = qnet16[x][y];
              qinward64[(4 - (c/25))*4 + x][(4 - (d/25))*4 + y] = qinward16[x][y];
              max_dvmdt64[(4 - (c/25))*4 + x][(4 - (d/25))*4 + y] = max_dvmdt16[x][y];
              vm_peak64[(4 - (c/25))*4 + x][(4 - (d/25))*4 + y] = vm_peak16[x][y];
              vm_dia64[(4 - (c/25))*4 + x][(4 - (d/25))*4 + y] = vm_dia16[x][y];
              v_apd90_64[(4 - (c/25))*4 + x][(4 - (d/25))*4 + y] = v_apd90_16[x][y];
              v_apd50_64[(4 - (c/25))*4 + x][(4 - (d/25))*4 + y] = v_apd50_16[x][y];
              apdtri64[(4 - (c/25))*4 + x][(4 - (d/25))*4 + y] = apdtri16[x][y];
              ca_peak64[(4 - (c/25))*4 + x][(4 - (d/25))*4 + y] = ca_peak16[x][y];
              ca_dia64[(4 - (c/25))*4 + x][(4 - (d/25))*4 + y] = ca_dia16[x][y];
              cad90_64[(4 - (c/25))*4 + x][(4 - (d/25))*4 + y] = cad90_16[x][y];
              cad50_64[(4 - (c/25))*4 + x][(4 - (d/25))*4 + y] = cad50_16[x][y];
              cadtri64[(4 - (c/25))*4 + x][(4 - (d/25))*4 + y] = cadtri16[x][y];
            }
          }
          printf("kimoo: qnet64\n");

      for(int x = 1; x < 64+1; x++){
         for (int y = 1; y < 64+1; y++) {
           qnet256[(4 - (a/25))*4 + x][(4 - (b/25))*4 + y] = qnet64[x][y];
           qinward256[(4 - (a/25))*4 + x][(4 - (b/25))*4 + y] = qinward64[x][y];
           max_dvmdt256[(4 - (a/25))*4 + x][(4 - (b/25))*4 + y] = max_dvmdt64[x][y];
           vm_peak256[(4 - (a/25))*4 + x][(4 - (b/25))*4 + y] = vm_peak64[x][y];
           vm_dia256[(4 - (a/25))*4 + x][(4 - (b/25))*4 + y] = vm_dia64[x][y];
           v_apd90_256[(4 - (a/25))*4 + x][(4 - (b/25))*4 + y] = v_apd90_64[x][y];
           v_apd50_256[(4 - (a/25))*4 + x][(4 - (b/25))*4 + y] = v_apd50_64[x][y];
           apdtri256[(4 - (a/25))*4 + x][(4 - (b/25))*4 + y] = apdtri64[x][y];
           ca_peak256[(4 - (a/25))*4 + x][(4 - (b/25))*4 + y] = ca_peak64[x][y];
           ca_dia256[(4 - (a/25))*4 + x][(4 - (b/25))*4 + y] = ca_dia64[x][y];
           cad90_256[(4 - (a/25))*4 + x][(4 - (b/25))*4 + y] = cad90_64[x][y];
           cad50_256[(4 - (a/25))*4 + x][(4 - (b/25))*4 + y] = cad50_64[x][y];
           cadtri256[(4 - (a/25))*4 + x][(4 - (b/25))*4 + y] = cadtri64[x][y];
         }
      }
      printf("kimoo: qnet256\n");

  //--------------------------------------------------------------------------------------------------------
  free(p_param);
  delete p_patch;

  printf("Create matrices\n");
  FILE *fp_matrix;
  fp_matrix = fopen( "qnet_matrix.mat", "w" );
  for (int y = 1; y < 256 + 1; y++) {
    for (int x = 1; x < 256 + 1; x++) {
      fprintf( fp_matrix, "%lf\n", qnet256[y][x] );
    }
  }
  fclose(fp_matrix);
/*  fp_matrix = fopen( "qinward_matrix.mat", "w" );
  for (int y = 1; y < 256 + 1; y++) {
    for (int x = 1; x < 256 + 1; x++) {
      fprintf( fp_matrix, "%lf\n", qinward256[y][x] );
    }
  }
  fclose(fp_matrix);
*/
  fp_matrix = fopen( "max_dvmdt_matrix.mat", "w" );
  for (int y = 1; y < 256 + 1; y++) {
    for (int x = 1; x < 256 + 1; x++) {
      fprintf( fp_matrix, "%lf\n", max_dvmdt256[y][x] );
    }
  }
  fclose(fp_matrix);
  fp_matrix = fopen( "vm_peak_matrix.mat", "w" );
  for (int y = 1; y < 256 + 1; y++) {
    for (int x = 1; x < 256 + 1; x++) {
      fprintf( fp_matrix, "%lf\n", vm_peak256[y][x] );
    }
  }
  fclose(fp_matrix);
  fp_matrix = fopen( "vm_dia_matrix.mat", "w" );
  for (int y = 1; y < 256 + 1; y++) {
    for (int x = 1; x < 256 + 1; x++) {
      fprintf( fp_matrix, "%lf\n", vm_dia256[y][x] );
    }
  }
  fclose(fp_matrix);
  fp_matrix = fopen( "v_apd90_matrix.mat", "w" );
  for (int y = 1; y < 256 + 1; y++) {
    for (int x = 1; x < 256 + 1; x++) {
      fprintf( fp_matrix, "%lf\n", v_apd90_256[y][x] );
    }
  }
  fclose(fp_matrix);
  fp_matrix = fopen( "v_apd50_matrix.mat", "w" );
  for (int y = 1; y < 256 + 1; y++) {
    for (int x = 1; x < 256 + 1; x++) {
      fprintf( fp_matrix, "%lf\n", v_apd50_256[y][x] );
    }
  }
  fclose(fp_matrix);
  fp_matrix = fopen( "apdtri_matrix.mat", "w" );
  for (int y = 1; y < 256 + 1; y++) {
    for (int x = 1; x < 256 + 1; x++) {
      fprintf( fp_matrix, "%lf\n", apdtri256[y][x] );
    }
  }
  fclose(fp_matrix);
  fp_matrix = fopen( "ca_peak_matrix.mat", "w" );
  for (int y = 1; y < 256 + 1; y++) {
    for (int x = 1; x < 256 + 1; x++) {
      fprintf( fp_matrix, "%lf\n", ca_peak256[y][x] );
    }
  }
  fclose(fp_matrix);
  fp_matrix = fopen( "ca_dia_matrix.mat", "w" );
  for (int y = 1; y < 256 + 1; y++) {
    for (int x = 1; x < 256 + 1; x++) {
      fprintf( fp_matrix, "%lf\n", ca_dia256[y][x] );
    }
  }
  fclose(fp_matrix);
  fp_matrix = fopen( "cad90_matrix.mat", "w" );
  for (int y = 1; y < 256 + 1; y++) {
    for (int x = 1; x < 256 + 1; x++) {
      fprintf( fp_matrix, "%lf\n", cad90_256[y][x] );
    }
  }
  fclose(fp_matrix);
  fp_matrix = fopen( "cad50_matrix.mat", "w" );
  for (int y = 1; y < 256 + 1; y++) {
    for (int x = 1; x < 256 + 1; x++) {
      fprintf( fp_matrix, "%lf\n", cad50_256[y][x] );
    }
  }
  fclose(fp_matrix);
  fp_matrix = fopen( "cadtri_matrix.mat", "w" );
  for (int y = 1; y < 256 + 1; y++) {
    for (int x = 1; x < 256 + 1; x++) {
      fprintf( fp_matrix, "%lf\n", cadtri256[y][x] );
    }
  }
  fclose(fp_matrix);

  system("mkdir result");
  system("mv *.plt result/");
  system("mv *.mat result/");


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

