#include "population_orudy.hpp"
#include "rhs.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <cvode/cvode.h>
#include <cvode/cvode_dense.h>
#include <nvector/nvector_serial.h>

#if defined ORUDY2011_STATIC34234
int population_ORudy(int argc, char **argv, param_t *p_param)
{
  /* CellML object pointer */
  CellML *p_cell;
  /* time related variables */
  double t_depolarize;
  double t_repolarize;
  double ap_duration;
  /* volatge related variables */
  double v_prev;
  double v_top;
  double v_valley;
  double v_apd90;
  /* supporting variables */
  bool is_ap_increasing;
  int print_freq;
  int err_code;
  FILE *fp_temp;
  FILE *fp_res;

  double apd3[3 + 1][3 + 1];
  double apd9[9 + 1][9 + 1];
  double apd27[27 + 1][27 + 1];
  double apd81[81 + 1][81 + 1];

  double mean_ca_i3[3 + 1][3 + 1];
  double mean_ca_i9[9 + 1][9 + 1];
  double mean_ca_i27[27 + 1][27 + 1];
  double mean_ca_i81[81 + 1][81 + 1];

  double sum_ca_i;
  double mean_ca_i;
  double scalars[10]; 

  /* SUNDIALs variables */
  int retval;
  int iout, imax;
  double tnext, tcurr;
  void *cvode_main;
  N_Vector states_vec;


 
  /* Initialize p_cell based on the ionmodel value */
  p_cell = init_cell( );
  if( p_cell == NULL ){
    return 20;
  }

  char scenario_name[255] = "";
  char result_file[100] = "";
  char command[255] = "";
  int a,b,c,d,e,f,g,h,m;

  for (int y = 1; y < 81 + 1; y++) {
    for (int x = 1; x < 81 + 1; x++) {
      apd81[y][x] = 0.0;
      mean_ca_i81[y][x] = 0.0;
    }
  }


#ifdef IS_EDISON == 1
  system("rm -rf result");
  system("mkdir result");
#endif
  a = 100;
  b = 100;
      for( c = 0; c <= 100; c += 50) {
        for( d = 0; d <= 100; d += 50) {
          for( e = 0; e <= 100; e += 50) {
            for( f = 0; f <= 100; f += 50) {
              for( g = 0; g <= 100; g += 50) {
                for( h = 0; h <= 100; h += 50) {

                  p_cell->initConsts();

                  scalars[0] = (double)a/(double)100;
                  scalars[1] = (double)b/(double)100;
                  scalars[2] = (double)c/(double)100;
                  scalars[3] = (double)d/(double)100;
                  scalars[4] = (double)e/(double)100;
                  scalars[5] = (double)f/(double)100;
                  scalars[6] = (double)g/(double)100;
                  scalars[7] = (double)h/(double)100;

                  p_cell->CONSTANTS[GNa] *= scalars[0];
                  p_cell->CONSTANTS[GNaL] *= scalars[1];      
                  p_cell->CONSTANTS[PCa] *= scalars[2];
                  p_cell->CONSTANTS[GKr] *= scalars[3];
                  p_cell->CONSTANTS[GKs] *= scalars[4];
                  p_cell->CONSTANTS[GK1] *= scalars[5];
                  p_cell->CONSTANTS[GKb] *= scalars[6];
                  p_cell->CONSTANTS[Gto] *= scalars[7];

                  p_cell->CONSTANTS[stim_period] = p_param->bcl_init;
                  p_cell->CONSTANTS[stim_amplitude] *= p_param->stim_amp;
                  p_cell->CONSTANTS[stim_duration] = p_param->stim_dur;
                  is_ap_increasing = false;
                  ap_duration = 0.0;
                  v_prev = 0.0;
                  v_top = 0.0;
                  v_valley = 0.0;
                  v_apd90 = p_param->v_apd90;
                  print_freq = (1. / p_param->dt) * p_param->dt_write;

                  // Create CVODE solver
                  cvode_main = CVodeCreate( CV_BDF,CV_NEWTON );
                  // Give p_cell as CVode User Data
                  CVodeSetUserData( cvode_main, p_cell );
                  // Create the states vector based on the STATES array
                  states_vec = N_VMake_Serial( p_cell->states_size, p_cell->STATES );
                  // Initalize CVODE solver
                  CVodeInit( cvode_main, rhs_fn, 0.0, states_vec );
                  // Set up the future time (must be more than 0.0)
                  tnext = p_param->dt;
                  // Set up the maximum step size (dt)
                  CVodeSetMaxStep( cvode_main, p_param->dt );
                  // Set up the linear solver type
                  CVDense( cvode_main, p_cell->states_size );
                  // Set up the numerical error tolerances
                  CVodeSStolerances( cvode_main, 1.0e-7, 1.0e-7 );

                  iout = 0;
                  imax = ( p_param->num_pace1 * p_cell->CONSTANTS[stim_period] ) / p_param->dt;

                  sum_ca_i = 0.0;
                  m = 0;
      
                  sprintf(scenario_name, 
                  "GNa%.2fGNaL%.2fPCa%.2fGKr%.2fGKs%.2fGK1%.2fGKb%.2fGto%.2f", 
                  scalars[0],
                  scalars[1],
                  scalars[2],
                  scalars[3],
                  scalars[4],
                  scalars[5],
                  scalars[6],
                  scalars[7]
                  );
                  sprintf( result_file, "%s-%s.plt", "elec_res", scenario_name );	
                  fp_res = fopen( result_file, "w" );
                  fprintf( fp_res, "%s %s %s\n",
                        "TIME","Vm","Ca_i" );

                  /* START SIMULATION CORE */
                  while(1) {
                    v_prev = p_cell->STATES[0];

                    // solve ODE
                    retval = CVode( cvode_main, tnext, states_vec, &tcurr, CV_NORMAL  );
                    if( retval == CV_SUCCESS ){
                      iout++;
                      tnext += p_param->dt;
                    }

                    /* find out whether the simulation can produce ACTION POTENTIAL or not */
                    if (p_cell->STATES[0] > v_prev && !is_ap_increasing) {
                      is_ap_increasing = true;
                      v_valley = v_prev;
                    }
                    else if (p_cell->STATES[0] < v_prev && is_ap_increasing) {
                      is_ap_increasing = false;
                      v_top = v_prev;
                    }
                    // repolarization, when the potential moves from HIGH to LOW
                    if (v_prev > v_apd90 && p_cell->STATES[0] <= v_apd90) {
                      t_repolarize = tcurr;
                      ap_duration = t_repolarize - t_depolarize;
                    }
                    // depolarization, when the potential moves from LOW to HIGH
                    if (v_prev <= v_apd90 && p_cell->STATES[0] > v_apd90) {
                      t_depolarize = tcurr;
                    }

                    // retrieve some informations from the last action potential
                    if( tnext >= (p_param->num_pace1-1) * p_cell->CONSTANTS[stim_period] ) {
                      fprintf( fp_res, "%lf %lf %lf\n",
                        tcurr,
                        p_cell->STATES[0],
                        p_cell->STATES[Ca_i] );
                        sum_ca_i += p_cell->STATES[Ca_i];
                        m++;
                    }
                  
                    if (iout >= imax) break;
                  }
                  /* END SIMULATION CORE */

                  // calculate both of the average Ca_i and i_NaK
                  m--;
                  mean_ca_i = sum_ca_i / m;

                  apd4[(g+50)/50][(h+50)/50] = ap_duration;
                  mean_ca_i4[g/25][h/25] = mean_ca_i;

                  fp_temp = fopen( "temp.plt", "a" );
                  if( ftell(fp_temp) == 0 ){
                    fprintf( fp_temp, "%s %s %s %s\n",
                              "scenario_name", "ap_duration", "mean_ca_i" );
                  }

                  fprintf( fp_temp, "%s %lf %lf\n",
                            scenario_name, 
                            ap_duration, 
                            mean_ca_i );
                  fclose(fp_temp);
                  fclose(fp_res);

                  p_cell->CONSTANTS[GNa] /= scalars[0];
                  p_cell->CONSTANTS[GNaL] /= scalars[1];
                  p_cell->CONSTANTS[PCa] /= scalars[2];
                  p_cell->CONSTANTS[GKr] /= scalars[3];
                  p_cell->CONSTANTS[GKs] /= scalars[4];
                  p_cell->CONSTANTS[GK1] /= scalars[5];
                  p_cell->CONSTANTS[GKb] /= scalars[6];
                  p_cell->CONSTANTS[Gto] /= scalars[7];

                }
              }
              for(int x = 1; x < 4+1; x++){
                for (int y = 1; y < 4+1; y++) {
                  apd16[(e/25 - 1) * 4 + x][(f/25 - 1) * 4 + y] = apd4[x][y];
                  mean_ca_i16[(e/25 - 1) * 4 + x][(f/25 - 1) * 4 + y] = mean_ca_i4[x][y];
                }
              }
              printf("kimoo: mapd16\n");
            }
          }
          for(int x=1; x< 16+1; x++){
            for(int y = 1; y < 16+1; y++) {
              apd64[(c/25 - 1) * 16 + x][(d/25 - 1) * 16 + y] = apd16[x][y];
              mean_ca_i64[(c/25 - 1) * 16 + x][(d/25 - 1) * 16 + y] = mean_ca_i16[x][y];
            }
          }
          printf("kimoo: mapd64\n");
        }
      }
      for(int x=1; x< 64+1; x++){ 
        for (int y = 1; y < 64+1; y++) {
          apd256[(a/25 - 1) * 64 + x][(b/25 - 1) * 64 + y] = apd64[x][y];
          mean_ca_i256[(a/25 - 1) * 64 + x][(b/25 - 1) * 64 + y] = mean_ca_i64[x][y];
        }
      }
      printf("kimoo: mapd256\n");


  printf("Create matrices\n");
  FILE *fp_matrix;
  fp_matrix = fopen( "apd_matrix.txt", "w" );
  for (int y = 1; y < 256 + 1; y++) {
    for (int x = 1; x < 256 + 1; x++) {
      fprintf( fp_matrix, "%lf\n", apd256[y][x] );
    }
  }
  fclose( fp_matrix );
  fp_matrix = fopen( "mean_ca_i_matrix.txt", "w" );
  for (int y = 1; y < 256+1; y++) {
    for (int x = 1; x < 256+1; x++) {
      fprintf( fp_matrix, "%lf\n", mean_ca_i256[y][x] );
    }
  }
  fclose( fp_matrix );

  printf("Finished creating matrices\n");


  delete p_cell;
  free( p_param );

  return 0;
}
#endif
