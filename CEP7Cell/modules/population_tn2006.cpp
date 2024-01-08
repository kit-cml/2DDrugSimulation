#include "population_tn2006.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>


#if defined TN2006ENDO || defined TN2006M || defined TN2006EPI
int population_TN2006(int argc, char **argv, param_t *p_param)
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
  FILE *fp_res;
  FILE *fp_log;

  double apd4[4 + 1][4 + 1];
  double apd16[16 + 1][16 + 1];
  double apd64[64 + 1][64 + 1];
  double apd256[256 + 1][256 + 1];
  double apd1024[1024 + 1][1024 + 1];

  double mean_ca_i4[4 + 1][4 + 1];
  double mean_ca_i16[16 + 1][16 + 1];
  double mean_ca_i64[64 + 1][64 + 1];
  double mean_ca_i256[256 + 1][256 + 1];
  double mean_ca_i1024[1024 + 1][1024 + 1];


  double sum_ca_i;
  double mean_ca_i;
  double scalars[10];
 
  /* SUNDIALs variables */
  int retval;
  int iout, imax;
  double tnext, tcurr;


 
  /* Initialize p_cell based on the ionmodel value */
  p_cell = init_cell( );
  if( p_cell == NULL ){
    return 20;
  }

  char scenario_name[255] = "";
  char result_file[100] = "";
  char command[255] = "";
  int a,b,c,d,e,f,g,h,i,j,m;

#ifdef IS_EDISON == 1
  system("rm -rf result");
  system("mkdir result");
#endif
  a = 100;
  b = 100;
      for( c = 25; c <= 100; c += 25) {
        for( d = 25; d <= 100; d += 25) {
          for( e = 25; e <= 100; e += 25) {
            for( f = 25; f <= 100; f += 25) {
              for( g = 25; g <= 100; g += 25) {
                for( h = 25; h <= 100; h += 25) {
                  for( i = 25; i <= 100; i += 25) {
                    for( j = 25; j <= 100; j += 25) { 

                      p_cell->initConsts();

                      scalars[0] = (double)a/(double)100;
                      scalars[1] = (double)b/(double)100;
                      scalars[2] = (double)c/(double)100;
                      scalars[3] = (double)d/(double)100;
                      scalars[4] = (double)e/(double)100;
                      scalars[5] = (double)f/(double)100;
                      scalars[6] = (double)g/(double)100;
                      scalars[7] = (double)h/(double)100;
                      scalars[8] = (double)i/(double)100;
                      scalars[9] = (double)j/(double)100;

                      p_cell->CONSTANTS[g_Ks] *= scalars[0];
                      p_cell->CONSTANTS[g_Kr] *= scalars[1];      
                      p_cell->CONSTANTS[g_K1] *= scalars[2];
                      p_cell->CONSTANTS[g_Na] *= scalars[3];
                      p_cell->CONSTANTS[g_bna] *= scalars[4];
                      p_cell->CONSTANTS[g_CaL] *= scalars[5];
                      p_cell->CONSTANTS[g_bca] *= scalars[6];
                      p_cell->CONSTANTS[g_to] *= scalars[7];
                      p_cell->CONSTANTS[g_pCa] *= scalars[8];
                      p_cell->CONSTANTS[g_pK] *= scalars[9];
	
                      p_cell->CONSTANTS[stim_period] = p_param->bcl_init;
                      p_cell->CONSTANTS[stim_amplitude] *= p_param->stim_amp;
                      p_cell->CONSTANTS[stim_duration] = p_param->stim_dur;
                      is_ap_increasing = false;
                      ap_duration = 0.0;
                      v_prev = 0.0;
                      v_top = 0.0;
                      v_prev = 0.0;
                      v_top = 0.0;
                      v_valley = 0.0;
                      v_apd90 = p_param->v_apd90;
                      print_freq = (1. / p_param->dt) * p_param->dt_write;

                      tcurr = 0.;
                      iout = 0;
                      imax = ( p_param->num_pace1 * p_cell->CONSTANTS[stim_period] ) / p_param->dt;

		      sum_ca_i = 0.0;
                      m = 0;
					
                      sprintf(scenario_name, 
                      "g_Ks%.2fg_Kr%.2fg_K1%.2fg_Na%.2fg_bna%.2fg_CaL%.2fg_bca%.2fg_to%.2fg_pCa%.2fg_pK%.2f", 
                      scalars[0],
                      scalars[1],
                      scalars[2],
                      scalars[3],
                      scalars[4],
                      scalars[5],
                      scalars[6],
                      scalars[7],
                      scalars[8],
                      scalars[9]);
                      sprintf( result_file, "%s-%s.plt", "elec_res", scenario_name );	
                      fp_res = fopen( result_file, "w" );
                      fprintf( fp_res, "%s %s %s\n",
                           "TIME","Vm","Ca_i" );

                      /* START SIMULATION CORE */
                      while(iout < imax) {
                       v_prev = p_cell->STATES[0];
                       
                       /* calculate ALGEBRAIC and RATES */
                       p_cell->computeRates( tcurr,
                                            p_cell->CONSTANTS,
                                            p_cell->RATES,
                                            p_cell->STATES,
                                            p_cell->ALGEBRAIC );
                       /* solve the RATES using Analytical Method */
                       p_cell->solveAnalytical(p_param->dt);
                       tcurr += p_param->dt;
 
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
                       if( tcurr >= (p_param->num_pace1-1) * p_cell->CONSTANTS[stim_period] ) {
                         fprintf( fp_res, "%lf %lf %lf\n",
                           tcurr,
                           p_cell->STATES[0],
                           p_cell->STATES[Ca_i] );
                         sum_ca_i += p_cell->STATES[Ca_i];
                         m++;
                       }
                    
                       iout++;
                     }
                      /* END SIMULATION CORE */

                      // calculate both of the average Ca_i and i_NaK
                      m--;
                      mean_ca_i = sum_ca_i / m;

                      apd4[i/25][j/25] = ap_duration;
                      mean_ca_i4[i/25][j/25] = mean_ca_i;

                      fp_log = fopen( "temp.plt", "a" );
                      if( ftell(fp_log) == 0 ){
                        fprintf( fp_log, "%s %s %s\n",
                                 "scenario_name", "ap_duration", "mean_ca_i" );
                      }
                      fprintf( fp_log, "%s %lf %lf\n",
                               scenario_name, 
                               ap_duration, 
                               mean_ca_i );
                      fclose(fp_log);
                      fclose(fp_res);

                      p_cell->CONSTANTS[g_Ks] /= scalars[0];
                      p_cell->CONSTANTS[g_Kr] /= scalars[1];
                      p_cell->CONSTANTS[g_K1] /= scalars[2];
                      p_cell->CONSTANTS[g_Na] /= scalars[3];
                      p_cell->CONSTANTS[g_bna] /= scalars[4];
                      p_cell->CONSTANTS[g_CaL] /= scalars[5];
                      p_cell->CONSTANTS[g_bca] /= scalars[6];
                      p_cell->CONSTANTS[g_to] /= scalars[7];
                      p_cell->CONSTANTS[g_pCa] /= scalars[8];
                      p_cell->CONSTANTS[g_pK] /= scalars[9];

                    }
                  }
                  for(int x = 1; x < 4+1; x++){
                    for (int y = 1; y < 4+1; y++) {
                      apd16[(g/25 - 1) * 4 + x][(h/25 - 1) * 4 + y] = apd4[x][y];
                      mean_ca_i16[(g/25 - 1) * 4 + x][(h/25 - 1) * 4 + y] = mean_ca_i4[x][y];
                    }
                  }
                  printf("kimoo: mapd16\n");
                }
              }
              for(int x=1; x < 16+1; x++){ 
                for(int y = 1; y < 16+1; y++) {
                  apd64[(e/25 - 1) * 16 + x][(f/25 - 1) * 16 + y] = apd16[x][y];
                  mean_ca_i64[(e/25 - 1) * 16 + x][(f/25 - 1) * 16 + y] = mean_ca_i16[x][y];
                }
              }
              printf("kimoo: mapd64\n");
            }
          }
          for(int x=1; x < 64+1; x++){ 
            for (int y = 1; y < 64+1; y++) {
              apd256[(c/25 - 1) * 64 + x][(d/25 - 1) * 64 + y] = apd64[x][y];
              mean_ca_i256[(c/25 - 1) * 64 + x][(d/25 - 1) * 64 + y] = mean_ca_i64[x][y];
            }  
          }
          printf("kimoo: mapd256\n");
        }
      }
      for(int x=1; x<256+1; x++){ 
        for (int y = 1; y < 256 + 1; y++) {
          apd1024[(a/25 - 1) * 256 + x][(b/25 - 1) * 256 + y] = apd256[x][y];
          mean_ca_i1024[(a/25 - 1) * 156 + x][(b/25 - 1) * 256 + y] = mean_ca_i256[x][y];
        }
      }
      printf("kimoo: mapd1024\n");


  printf("Create matrices\n");
  FILE *fp_matrix;
  fp_matrix = fopen( "apd_matrix.txt", "w" );
  for (int y = 1; y < 1024 + 1; y++) {
    for (int x = 1; x < 1024 + 1; x++) {
      fprintf( fp_matrix, "%lf\n", apd1024[y][x] );
    }
  }
  fclose(fp_matrix);
  fp_matrix = fopen( "mean_ca_i_matrix.txt", "w" );
  for (int y = 1; y < 1024 + 1; y++) {
    for (int x = 1; x < 1024 + 1; x++) {
      fprintf( fp_matrix, "%lf\n", mean_ca_i1024[y][x] );
    }
  }
  fclose( fp_matrix );
  printf("Finished creating matrices\n");


  delete p_cell;
  free( p_param );

  return 0;
}
#endif


