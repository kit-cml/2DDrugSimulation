#include "patches/tentusscher_noble_noble_panfilov_2004_b.hpp"

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

void population_ANM(int argc, char **argv);
void edison_assign_params_single(int argc, char *args[], param_t *p_param);
double get_ANM( double *apd_arr, int length );
double average( double *apd_arr, int start_idx, int end_idx );

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
  population_ANM(argc, argv);

  return 0;
}



void population_ANM(int argc, char **argv)
{
  int a,b,c,d,e,f,g,h,i,j;
  int start_val, end_val, inc;
  double scalars[10];
  char scenario_name[255] = "";
  char result_file[100] = "";
  FILE *fp_log;


  double aocl4[4+1][4+1];
  double aocl16[16+1][16+1];
  double aocl64[64+1][64+1];
  double aocl256[256+1][256+1];
  double aocl1024[1024+1][1024+1];

  double apdb4[4+1][4+1];
  double apdb16[16+1][16+1];
  double apdb64[64+1][64+1];
  double apdb256[256+1][256+1];
  double apdb1024[1024+1][1024+1];

  double apde4[4+1][4+1];
  double apde16[16+1][16+1];
  double apde64[64+1][64+1];
  double apde256[256+1][256+1];
  double apde1024[1024+1][1024+1];

  double apdb_aocl4[4+1][4+1];
  double apdb_aocl16[16+1][16+1];
  double apdb_aocl64[64+1][64+1];
  double apdb_aocl256[256+1][256+1];
  double apdb_aocl1024[1024+1][1024+1];

  double apde_aocl4[4+1][4+1];
  double apde_aocl16[16+1][16+1];
  double apde_aocl64[64+1][64+1];
  double apde_aocl256[256+1][256+1];
  double apde_aocl1024[1024+1][1024+1];

  /* patch_clamp object pointer */
  patch_clamp *p_patch;
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
  double apdb;
  double apde;
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
  bool is_beginning;
  /* SUNDIALs variables */
  int retval;
  double t_out, t_curr;
  void *cvode_mem;
  N_Vector states_vec;

 

  //--------------------------------------------------------------------------------------------------------

  param_t *p_param;
  p_param = (param_t*)malloc( sizeof( param_t ) );
  edison_assign_params_single(argc,argv,p_param);

  // PLEASE DON'T PUT THIS CALL INSIDE THE DEEPEST LOOP!!!
  p_patch = new tentusscher_noble_noble_panfilov_2004_b();

  a = p_param->gks_scale;
  b = p_param->gkr_scale;
  c = p_param->gk1_scale;
  d = p_param->gna_scale;
          for(e = 100; e >= 25; e -= 25){ // 16 hour
            for(f = 100; f >= 25; f -= 25){ //4 hour
              for(g = 100; g >= 25; g -= 25){ // 1hour
                for(h = 100; h >= 25; h -= 25){
                  for(i = 100; i >= 25; i -= 25){
                    for(j = 100; j >= 25; j -= 25){


                      p_patch->initConsts();
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

                      p_patch->CONSTANTS[g_Ks] *= scalars[0];
                      p_patch->CONSTANTS[g_Kr] *= scalars[1];
                      p_patch->CONSTANTS[g_K1] *= scalars[2];
                      p_patch->CONSTANTS[g_Na] *= scalars[3];
                      p_patch->CONSTANTS[g_bna] *= scalars[4];
                      p_patch->CONSTANTS[g_CaL] *= scalars[5];
                      p_patch->CONSTANTS[g_bca] *= scalars[6];
                      p_patch->CONSTANTS[g_to] *= scalars[7];
                      p_patch->CONSTANTS[g_pCa] *= scalars[8];
                      p_patch->CONSTANTS[g_pK] *= scalars[9];



                      sprintf(scenario_name,
                      "g_Ks%.2fg_Kr%.2fg_K1%.2fg_Na%.2fg_bna%.2fg_CaL%.2fg_bca%.2fg_to%.2fg_pCa%.2fg_pK%.2f",
                      scalars[0],scalars[1],scalars[2],scalars[3],scalars[4],scalars[5],scalars[6],scalars[7],scalars[8],scalars[9]);

                      cvode_mem = CVodeCreate( CV_BDF,CV_NEWTON );
                      CVodeSetUserData( cvode_mem, p_patch );
                      states_vec = N_VMake_Serial( p_patch->states_size, p_patch->STATES );
                      CVodeInit( cvode_mem, rhs_fn, 0., states_vec );
                      t_out = 0.1;
                      CVodeSetMaxStep( cvode_mem, p_param->dt );
                      CVDense( cvode_mem, p_patch->states_size );
                      CVodeSStolerances( cvode_mem, 1.0e-7, 1.0e-7 );

                      bcl_curr = p_param->bcl_init;
                      bcl_dec = p_param->bcl_decrement;
                      p_patch->CONSTANTS[stim_period] = p_param->bcl_init;
                      p_patch->CONSTANTS[stim_amplitude] *= p_param->stim_amp;
                      p_patch->CONSTANTS[stim_duration] = p_param->stim_dur;
                      is_ap_increasing = false;
                      is_action_potential = false;
                      t_curr = 0.0;
                      t_prev = 0.0;
                      t_stim_start = 0.0;
                      t_stim_end = t_stim_start + p_patch->CONSTANTS[stim_duration];
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
                      is_beginning = true;

                      fp_apdr = fopen( "apdr.plt", "w" );
                      sprintf( result_file, "anm_%s.plt", scenario_name );
                      fp_anm = fopen( result_file, "w" );
                      fp_apd90 = fopen( "debug.plt", "w" );
                      sprintf( result_file, "vmcheck_%s.plt", scenario_name );
                      fp_vm = fopen( result_file, "w" );

                      int iout = 0;
                      /* Print the header of the output files */
                      fprintf( fp_vm, "%s %s\n", "TIME", "Vm" );


                     /* core simulation loop */
                     while (bcl_curr >= p_param->bcl_end) {
                       v_prev = p_patch->STATES[0];
                       // Control the flow of simulation.
                       // In case of APDR simulation, this part will also
                       // decrease the bcl based on the bcl_decrement value.
                       if (bcl_curr <= 300.0) bcl_dec = 10;
                       if (t_curr <= t_prev + bcl_curr * p_param->num_pace1) {
                         if (t_curr > t_stim_end) {
                           t_stim_start += bcl_curr;
                           t_stim_end = t_stim_start + p_patch->CONSTANTS[stim_duration];
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
                           if(anm > 0.05){
                             printf("Onsite cycle length: %lf\n", bcl_curr);
                             is_alternant = true;
                           }
                           if(alternant_counter > 4) is_alternant = true;
                         }
                         fprintf(fp_anm, "BCL: %lf ANM: %lf\n", bcl_curr, anm);

                         if( is_beginning == true ){
                           printf("BCL: %lf APDB: %lf SHOULD ONLY ONCE\n", bcl_curr, apd );
                           apdb = apd;
                           is_beginning = false;
                         }
                         if( is_alternant == true ){
                           printf("BCL: %lf APDE: %lf SHOULD ONLY ONCE\n", bcl_curr, apd );
                           apde = apd;
                           break;
                         }
                         bcl_curr -= bcl_dec;
                         if (bcl_curr >= 750.0) bcl_curr = 500.0;
                         if (bcl_curr <= 300.0) bcl_dec = 10;
                         p_patch->CONSTANTS[stim_period] = bcl_curr;
                         t_prev = t_curr;
                         apd_vec_size = 0;
                         alternant_counter = 0;
                       }

                       retval = CVode( cvode_mem, t_out, states_vec, &t_curr, CV_NORMAL  );
                       if( retval == CV_SUCCESS ){
                         iout++;
                         t_out += p_param->dt;
                       }

                       /* find out whether the simulation can produce ACTION POTENTIAL or not */
                       if (p_patch->STATES[0] > v_prev && !is_ap_increasing) {
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
                       else if (p_patch->STATES[0] < v_prev && is_ap_increasing) {
                         is_ap_increasing = false;
                         v_top = v_prev;
                       }
                       // repolarization, when the potential moves from HIGH to LOW
                       if (v_prev > v_apd90 && p_patch->STATES[0] <= v_apd90) {
                         t_repolarize = t_curr;
                         apd = t_repolarize - t_depolarize;
                         if(apd > bcl_curr) {
                           printf("WARNING!! APD is larger than BCL!! Alternant detected!!\n");
                           alternant_counter++;
                         }
                         apd_vec[apd_vec_size] = apd;
                         //printf("APD VEC %d: %lf\n", apd_vec_size,  apd_vec[apd_vec_size]);
                         apd_vec_size++;
                       }
                       // depolarization, when the potential moves from LOW to HIGH
                       if (v_prev <= v_apd90 && p_patch->STATES[0] > v_apd90) {
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
                         fprintf( fp_vm, "%lf %lf\n", t_curr, p_patch->STATES[V] );
                       }



                     }

                      fp_log = fopen( "temp.plt", "a" );
                      if( ftell(fp_log) == 0 )fprintf( fp_log, "%s %s %s %s %s %s\n","scenario_name", "AOCL", "APD_b", "APD_e", "AOCL/APD_b", "AOCL/APD_e" );
                      fprintf( fp_log, "%s %lf %lf %f %lf %lf\n",scenario_name, bcl_curr, apdb, apde, bcl_curr/apdb, bcl_curr/apde );


                      fclose( fp_vm );
                      fclose( fp_anm );
                      fclose( fp_apd90 );
                      fclose( fp_apdr );
                      fclose( fp_log );

                      N_VDestroy(states_vec);
                      CVodeFree(&cvode_mem);
                      delete []apd_vec;

                      aocl4[i/25][j/25] = bcl_curr;
                      apdb4[i/25][j/25] = apdb;
                      apde4[i/25][j/25] = apde;
                      apdb_aocl4[i/25][j/25] = bcl_curr/apdb;
                      apde_aocl4[i/25][j/25] = bcl_curr/apde;
                      printf("kimoo: aocl4\n");
                    }
                  }

                  for(int x = 1; x < 4+1; x++){
                    for (int y = 1; y < 4+1; y++) {
                      aocl16[(4 - (g/25))*4 + x][(4 - (h/25))*4 + y] = aocl4[x][y];
                      apdb16[(4 - (g/25))*4 + x][(4 - (h/25))*4 + y] = apdb4[x][y];
                      apde16[(4 - (g/25))*4 + x][(4 - (h/25))*4 + y] = apde4[x][y];
                      apdb_aocl16[(4 - (g/25))*4 + x][(4 - (h/25))*4 + y] = apdb_aocl4[x][y];
                      apde_aocl16[(4 - (g/25))*4 + x][(4 - (h/25))*4 + y] = apde_aocl4[x][y];
                    }
                  }
                  printf("kimoo: aocl16\n");

                }
              }

              for(int x = 1; x < 16+1; x++){
                for (int y = 1; y < 16+1; y++) {
                  aocl64[(4 - (e/25))*4 + x][(4 - (f/25))*4 + y] = aocl16[x][y];
                  apdb64[(4 - (e/25))*4 + x][(4 - (f/25))*4 + y] = apdb16[x][y];
                  apde64[(4 - (e/25))*4 + x][(4 - (f/25))*4 + y] = apde16[x][y];
                  apdb_aocl64[(4 - (e/25))*4 + x][(4 - (f/25))*4 + y] = apdb_aocl16[x][y];
                  apde_aocl64[(4 - (e/25))*4 + x][(4 - (f/25))*4 + y] = apde_aocl16[x][y];
                }
              }
              printf("kimoo: aocl64\n");

            }
          }

          for(int x = 1; x < 64+1; x++){
            for (int y = 1; y < 64+1; y++) {
              aocl256[(4 - (c/25))*4 + x][(4 - (d/25))*4 + y] = aocl64[x][y];
              apdb256[(4 - (c/25))*4 + x][(4 - (d/25))*4 + y] = apdb64[x][y];
              apde256[(4 - (c/25))*4 + x][(4 - (d/25))*4 + y] = apde64[x][y];
              apdb_aocl256[(4 - (c/25))*4 + x][(4 - (d/25))*4 + y] = apdb_aocl64[x][y];
              apde_aocl256[(4 - (c/25))*4 + x][(4 - (d/25))*4 + y] = apde_aocl64[x][y];
            }
          }
          printf("kimoo: aocl256\n");

       
     

      for(int x = 1; x < 256+1; x++){
         for (int y = 1; y < 256+1; y++) {
           aocl1024[(4 - (a/25))*4 + x][(4 - (b/25))*4 + y] = aocl256[x][y];
           apdb1024[(4 - (a/25))*4 + x][(4 - (b/25))*4 + y] = apdb256[x][y];
           apde1024[(4 - (a/25))*4 + x][(4 - (b/25))*4 + y] = apde256[x][y];
           apdb_aocl1024[(4 - (a/25))*4 + x][(4 - (b/25))*4 + y] = apdb_aocl256[x][y];
           apde_aocl1024[(4 - (a/25))*4 + x][(4 - (b/25))*4 + y] = apde_aocl256[x][y];
         }
      }
      printf("kimoo: aocl1024\n");

  //--------------------------------------------------------------------------------------------------------
  free(p_param);
  delete p_patch;

  printf("Create matrices\n");
  FILE *fp_matrix;
  fp_matrix = fopen( "aocl_matrix.mat", "w" );
  for (int y = 1; y < 1024 + 1; y++) {
    for (int x = 1; x < 1024 + 1; x++) {
      fprintf( fp_matrix, "%lf\n", aocl1024[y][x] );
    }
  }
  fclose(fp_matrix);
/*
  fp_matrix = fopen( "apdb_matrix.mat", "w" );
  for (int y = 1; y < 1024 + 1; y++) {
    for (int x = 1; x < 1024 + 1; x++) {
      fprintf( fp_matrix, "%lf\n", apdb1024[y][x] );
    }
  }
  fclose(fp_matrix);

  fp_matrix = fopen( "apde_matrix.mat", "w" );
  for (int y = 1; y < 1024 + 1; y++) {
    for (int x = 1; x < 1024 + 1; x++) {
      fprintf( fp_matrix, "%lf\n", apde1024[y][x] );
    }
  }
  fclose(fp_matrix);

  fp_matrix = fopen( "apdb_aocl_matrix.mat", "w" );
  for (int y = 1; y < 1024 + 1; y++) {
    for (int x = 1; x < 1024 + 1; x++) {
      fprintf( fp_matrix, "%lf\n", apdb_aocl1024[y][x] );
    }
  }
  fclose(fp_matrix);

  fp_matrix = fopen( "apde_aocl_matrix.mat", "w" );
  for (int y = 1; y < 1024 + 1; y++) {
    for (int x = 1; x < 1024 + 1; x++) {
      fprintf( fp_matrix, "%lf\n", apde_aocl1024[y][x] );
    }
  }
  fclose(fp_matrix);
*/
  system("mkdir result");
  system("mv *.plt result");
  system("mv *.mat result");

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
    else if (strcasecmp(key, "GKs_Scale") == 0) {
      p_param->gks_scale =  strtod(value, NULL) ;
      printf( "%s -- %lf\n", "GKs_Scale", p_param->gks_scale);
    }
    else if (strcasecmp(key, "GKr_Scale") == 0) {
      p_param->gkr_scale =  strtod(value, NULL) ;
      printf( "%s -- %lf\n", "GKr_Scale", p_param->gkr_scale);
    }
    else if (strcasecmp(key, "GK1_Scale") == 0) {
      p_param->gk1_scale =  strtod(value, NULL) ;
      printf( "%s -- %lf\n", "GK1_Scale", p_param->gk1_scale);
    }
    else if (strcasecmp(key, "GNa_Scale") == 0) {
      p_param->gna_scale =  strtod(value, NULL) ;
      printf( "%s -- %lf\n", "GNa_Scale", p_param->gna_scale);
    }
    else if (strcasecmp(key, "Is_Print_Vm") == 0) {
      p_param->is_print_vm = strtol( value, NULL, 10 );
      printf( "%s -- %d\n", "Is_Print_Vmcheck", p_param->is_print_vm );
    }

  }
  fclose( fp_inputdeck );
}

double get_ANM( double *apd_arr, int length )
{
  double alternant_magnitude;
  double change_magnitude;
  double mean_apd;
  int start;

  start = length - 10 - 1;
  change_magnitude = 0.0;

  printf("Start_idx: %d length_idx: %d\n", start, length);
  for ( int i = start; i < length-1; i++ ) {
    change_magnitude += fabs(apd_arr[i] - apd_arr[i+1]);
  }

  alternant_magnitude = change_magnitude / 10.0;
  mean_apd = average(apd_arr, start, length-1);
  //printf("Chg_mgntd: %lf Alt_mgntd: %lf mean_apd: %lf\n", change_magnitude, alternant_magnitude, mean_apd);
  return alternant_magnitude / mean_apd;
}

double average( double *apd_arr, int start_idx, int end_idx )
{
  double sum = 0.0;
  for (int i = start_idx; i <= end_idx; i++ ) {
    sum += apd_arr[i];
  }
  return sum / (end_idx - start_idx + 1);
}

