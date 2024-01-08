#include "population_ep.hpp"

#if defined TN2006ENDO || defined TN2006M || defined TN2006EPI
int population_EP(int argc, char **argv, param_t *p_param)
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


