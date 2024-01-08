#include "commons.hpp"
#include "globals.hpp"
#include "../cellmodels/patch_clamp.hpp"
#include "../libs/zip.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <sstream>
#include <vector>

#include <dirent.h>
#include <fcntl.h>
#include <magic.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

using std::string;
using std::stringstream;
using std::vector;


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

int is_file_exists (char *filename) {
  struct stat   buffer;
  return (stat (filename, &buffer) == 0);
}

void set_default_values( param_t *p_param )
{
  p_param->celltype = 0;
  p_param->simulation_mode = 0;
  p_param->t_write_vtk = 0.0;
  p_param->dt = 0.01;
  p_param->dt_write = 1.0;
  p_param->stim_amp = 1.0;
  p_param->stim_dur = 0.5;
  p_param->num_pace1 = 3;
  p_param->bcl_decrement = 50.0;
  p_param->bcl_end = 2000000.0;
  p_param->bcl_init = 1000.0;
  p_param->is_print_graph = 1;

  p_param->gks_scale =  1. ;
  p_param->gkr_scale =  1. ;
  p_param->gk1_scale =  1. ;
  p_param->gcal_scale = 1. ;
  p_param->gna_scale =  1. ;
  p_param->gnal_scale =  1. ;
  p_param->gbna_scale = 1. ;
  p_param->gbca_scale = 1. ;
  p_param->gto_scale = 1. ;
  p_param->gpca_scale = 1. ;
  p_param->gpk_scale = 1. ;
  
  p_param->is_s1s2 = 0;
  p_param->is_using_output = 0;
  p_param->is_dutta = 0;
  sprintf( p_param->variant, "%s", "ori" );
  sprintf( p_param->drug_name, "%s", "bepridil" );
  p_param->cmax = 33;
  sprintf( p_param->concs, "%s", "0.0" );

  glob_var::A1656D_mode = 0;
  glob_var::is_hrv = false;

}


void edison_assign_params_single(int argc, char *args[], param_t *p_param)
{
  char buffer[100];
  char key[100];
  char value[100];
  char file_name[255];
  FILE *fp_inputdeck;

  for (int i = 1; i < argc; i += 2) {
    if (!strcasecmp(args[i], "-input_deck"))
      strcpy(file_name, args[i + 1]);
    else if (!strcasecmp(args[i], "-hill_file"))
      strcpy(p_param->hill_file, args[i + 1]);
    else if (!strcasecmp(args[i], "-herg_file"))
      strcpy(p_param->herg_file, args[i + 1]);
    else if (!strcasecmp(args[i], "-hrv_files"))
      strcpy(p_param->hrv_files, args[i + 1]);
  }

  if(strlen(p_param->hill_file) != 0){
    if(mympi::rank == 0)  printf( "Hill File: %s\n", p_param->hill_file );
  }
  if(strlen(p_param->herg_file) != 0){
    if(mympi::rank == 0) printf( "hERG File: %s\n", p_param->herg_file );
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
      if(mympi::rank == 0) printf( "%s -- %lf\n", "Celltype", p_param->celltype);
    }
    else if (strcasecmp(key, "Variant") == 0) {
      strcpy( p_param->variant, value );
      if(mympi::rank == 0) printf( "%s -- %s\n", "Variant", p_param->variant);
    }
    else if (strcasecmp(key, "Simulation_Mode") == 0) {
      p_param->simulation_mode = strtol( value, NULL, 10 );
      if(mympi::rank == 0) printf( "%s -- %d\n", "Simulation_Mode", p_param->simulation_mode);
    }
    else if (strcasecmp(key, "Write_Time") == 0) {
      p_param->t_write_vtk = strtod( value, NULL );
      if(mympi::rank == 0) printf( "%s -- %lf\n", "Write_Time", p_param->t_write_vtk);
    }
    else if (strcasecmp(key, "Number_of_Pacing") == 0) {
      p_param->num_pace1 = strtol( value, NULL, 10 );
      if(mympi::rank == 0) printf( "%s -- %d\n", "Number_of_Pacing", p_param->num_pace1);
    }
    else if (strcasecmp(key, "Basic_Cycle_Length") == 0) {
      p_param->bcl_init = strtod( value, NULL );
      if(mympi::rank == 0) printf( "%s -- %lf\n", "Basic_Cycle_Length", p_param->bcl_init);
    }
    else if (strcasecmp(key, "Time_Step") == 0) {
      p_param->dt = strtod( value, NULL );
      if(mympi::rank == 0) printf( "%s -- %lf\n", "Time_Step", p_param->dt);
    }
    else if (strcasecmp(key, "Write_Step") == 0) {
      p_param->dt_write = strtod( value, NULL );
      if(mympi::rank == 0) printf( "%s -- %lf\n", "Writing_Step", p_param->dt_write);
    }
    else if (strcasecmp(key, "Stim_Amplitude") == 0) {
      p_param->stim_amp = strtod( value, NULL );
      if(mympi::rank == 0) printf( "%s -- %lf\n", "Stim_Amplitude", p_param->stim_amp);
    }
    else if (strcasecmp(key, "Stim_Duration") == 0) {
      p_param->stim_dur = strtod( value, NULL );
      if(mympi::rank == 0) printf( "%s -- %lf\n", "Stim_Duration", p_param->stim_dur);
    }
    else if (strcasecmp(key, "BCL_Step") == 0) {
      p_param->bcl_decrement = strtod( value, NULL );
      if(mympi::rank == 0) printf( "%s -- %lf\n", "BCL_Step", p_param->bcl_decrement);
    }
    else if (strcasecmp(key, "BCL_End") == 0) {
      p_param->bcl_end = strtod( value, NULL );
      if(mympi::rank == 0) printf( "%s -- %lf\n", "BCL_End", p_param->bcl_end);
    }
    else if (strcasecmp(key, "BCL_Print_Start") == 0) {
      p_param->bcl_print_start = strtod( value, NULL );
      if(mympi::rank == 0) printf( "%s -- %lf\n", "BCL_Print_Start", p_param->bcl_print_start);
    }
    else if (strcasecmp(key, "BCL_Print_End") == 0) {
      p_param->bcl_print_end = strtod( value, NULL );
      if(mympi::rank == 0) printf( "%s -- %lf\n", "BCL_Print_End", p_param->bcl_print_end);
    }

    else if (strcasecmp(key, "Scaled_Gate") == 0) {
      p_param->scaled_gate = strtol( value, NULL, 10 );
      if(mympi::rank == 0) printf( "%s -- %d\n", "Scaled_Gate", p_param->scaled_gate);
    }
    else if (strcasecmp(key, "Drug_Name") == 0) {
      strcpy( p_param->drug_name, value );
      if(mympi::rank == 0) printf( "%s -- %s\n", "Drug_Name", p_param->drug_name);
    }
    else if (strcasecmp(key, "Concentrations") == 0) {
      strcpy( p_param->concs, "0," );
      strcat( p_param->concs, value );
      if(mympi::rank == 0) printf( "%s -- %s\n", "Concentrations", p_param->concs);
    }

    else if (strcasecmp(key, "GKs_Scale") == 0) {
      p_param->gks_scale =  strtod(value, NULL) ;
      if(mympi::rank == 0) printf( "%s -- %lf\n", "GKs_Scale", p_param->gks_scale);
    }
    else if (strcasecmp(key, "GKr_Scale") == 0) {
      p_param->gkr_scale =  strtod(value, NULL) ;
      if(mympi::rank == 0) printf( "%s -- %lf\n", "GKr_Scale", p_param->gkr_scale);
    }
    else if (strcasecmp(key, "GCaL_Scale") == 0) {
      p_param->gcal_scale =  strtod(value, NULL) ;
      if(mympi::rank == 0) printf( "%s -- %lf\n", "GCaL_Scale", p_param->gcal_scale);
    }
    else if (strcasecmp(key, "GNa_Scale") == 0) {
      p_param->gna_scale =  strtod(value, NULL) ;
      if(mympi::rank == 0) printf( "%s -- %lf\n", "GNa_Scale", p_param->gna_scale);
    }
    else if (strcasecmp(key, "GK1_Scale") == 0) {
      p_param->gk1_scale =  strtod(value, NULL) ;
      if(mympi::rank == 0) printf( "%s -- %lf\n", "GK1_Scale", p_param->gk1_scale);
    }
    else if (strcasecmp(key, "GbNa_Scale") == 0) {
      p_param->gbna_scale =  strtod(value, NULL) ;
      if(mympi::rank == 0) printf( "%s -- %lf\n", "GbNa_Scale", p_param->gbna_scale);
    }
    else if (strcasecmp(key, "GbCa_Scale") == 0) {
      p_param->gbca_scale =  strtod(value, NULL) ;
      if(mympi::rank == 0) printf( "%s -- %lf\n", "GbCa_Scale", p_param->gbca_scale);
    }
    else if (strcasecmp(key, "Gto_Scale") == 0) {
      p_param->gto_scale =  strtod(value, NULL) ;
      if(mympi::rank == 0) printf( "%s -- %lf\n", "Gto_Scale", p_param->gto_scale);
    }
    else if (strcasecmp(key, "GpCa_Scale") == 0) {
      p_param->gpca_scale =  strtod(value, NULL) ;
      if(mympi::rank == 0) printf( "%s -- %lf\n", "GpCa_Scale", p_param->gpca_scale);
    }
    else if (strcasecmp(key, "GpK_Scale") == 0) {
      p_param->gpk_scale =  strtod(value, NULL) ;
      if(mympi::rank == 0) printf( "%s -- %lf\n", "GpK_Scale", p_param->gpk_scale);
    }
    else if (strcasecmp(key, "GNaL_Scale") == 0) {
      p_param->gnal_scale =  strtod(value, NULL) ;
      if(mympi::rank == 0) printf( "%s -- %lf\n", "GNaL_Scale", p_param->gnal_scale);
    }

    else if (strcasecmp(key, "Is_Print_Graph") == 0) {
      p_param->is_print_graph = strtol( value, NULL, 10 );
      if(mympi::rank == 0) printf( "%s -- %d\n", "Is_Print_Graph", p_param->is_print_graph );
    }
    else if (strcasecmp(key, "Is_Dutta") == 0) {
      p_param->is_dutta = strtol( value, NULL, 10 );
      if(mympi::rank == 0) printf( "%s -- %d\n", "Is_Dutta", p_param->is_dutta );
    }
    else if (strcasecmp(key, "Is_Using_Output") == 0) {
      p_param->is_using_output = strtol( value, NULL, 10 );
      if(mympi::rank == 0) printf( "%s -- %d\n", "Is_Using_Output", p_param->is_using_output );
    }
    else if (strcasecmp(key, "A1656D_Type") == 0) {
      glob_var::A1656D_mode = strtol( value, NULL, 10 );
      if(mympi::rank == 0) printf( "%s -- %d\n", "A1656D_Type", glob_var::A1656D_mode );
    }

  }
  fclose( fp_inputdeck );
}

vector<double> get_vector_from_string( string delim_str )
{
  vector<double> result;
  string substr;
  stringstream ss(delim_str);

  while( ss.good() )
  {
    getline(ss, substr, ',');
    result.push_back(strtod(substr.c_str(), NULL));
  }
  return result;
}

int extract_zip(const char *zip_filename, const char *target_folder)
{
  return zip_extract( zip_filename, target_folder, NULL, NULL);
}

void create_zip(const char *zip_filename, const char *target_folder)
{
  struct zip_t *zp_result = zip_open(zip_filename, ZIP_DEFAULT_COMPRESSION_LEVEL, 'w');
  if( zp_result != NULL ){
    zip_walk(zp_result, target_folder);
  }
  else{
    fprintf(stderr, "Error when zipping result folder!!\n");
    return;
  }
  if( zp_result != NULL ) zip_close(zp_result);
}

void zip_walk(struct zip_t *zip, const char *path) {
    DIR *dir;
    struct dirent *entry;
    char fullpath[MAX_PATH];
    struct stat s;

    memset(fullpath, 0, MAX_PATH);
    dir = opendir(path);
    if(dir == NULL){
      fprintf(stderr, "Error opening the directory %s\n", path);
      return;
    }


    while ((entry = readdir(dir))) {
      // skip "." and ".."
      if (!strcasecmp(entry->d_name, ".\0") || !strcasecmp(entry->d_name, "..\0"))
        continue;

      snprintf(fullpath, sizeof(fullpath), "%s/%s", path, entry->d_name);
      stat(fullpath, &s);
      if (S_ISDIR(s.st_mode))
        zip_walk(zip, fullpath);
      else {
        zip_entry_open(zip, fullpath);
        zip_entry_fwrite(zip, fullpath);
        zip_entry_close(zip);
      }
    }

    closedir(dir);
}

int convert_base_n(int number,int base){
    if(number == 0 || base==10)
        return number;
    return (number % base) + 10*convert_base_n(number / base, base);
}


void remove_dir_content(const char *path)
{
    struct dirent *de;
    char fname[300];
    DIR *dr = opendir(path);
    if(dr == NULL){
      printf("No file or directory found\n");
      return;
    }
    while((de = readdir(dr)) != NULL){
      int ret = -1;
      struct stat statbuf;
      sprintf(fname,"%s/%s",path,de->d_name);
      if(!strcasecmp(de->d_name, ".") || !strcasecmp(de->d_name, "..")) continue;
      if(!stat(fname, &statbuf)){
        if(S_ISDIR(statbuf.st_mode)){
          //printf("Is dir: %s\n",fname);
          if(ret != 0){
            remove_dir_content(fname);
            rmdir(fname);
          }
        }
        else{
          unlink(fname);
          //printf("Is file: %s\n",fname);
          //printf("Err: %d\n",unlink(fname));
        }
      }
    }
    closedir(dr);
}

double get_ANM( vector<double> vec )
{
  double alternant_magnitude;
  double change_magnitude;
  double mean_apd;
  int start;

  start = vec.size()-10-1;
  change_magnitude = 0.0;

  for ( int i = start; i < vec.size()-1; i++ ) {
    change_magnitude += fabs(vec[i] - vec[i+1]);
  }

  alternant_magnitude = change_magnitude / 10.0;
  mean_apd = average(vec, start, vec.size()-1);
  return alternant_magnitude / mean_apd;
}


double average( vector<double> vec, int start_idx, int end_idx )
{
  if( end_idx > vec.size() ){
    fprintf( stderr, "end_indx %d exceeded vector size %d\n", end_idx, vec.size() );
    return 0.;
  }
  double sum = 0.0;

  for (int idx = start_idx; idx <= end_idx; idx++ ) {
    sum += vec[idx];
  }

  return sum / (end_idx - start_idx + 1);
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

  char buffer[150] = {'\0'};

  fp_apdr = fopen( "result/apdr.plt", "r" );
  fp_apdr_last_2_cycles = fopen( "result/apdr2.plt", "w" );
  fp_apdr_di = fopen( "result/apdr_di.plt", "w" );
  fp_apdr_bcl = fopen( "result/apdr_bcl.plt", "w" );

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

  printf("APDR processing is successful!!!\n");
}
