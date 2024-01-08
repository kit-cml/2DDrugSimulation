// Code located in this project directory
#include "hrv_bench.hpp"
#include "commons.hpp"
#include "globals.hpp"
#include "../cellmodels/Ohara_Rudy_2011.hpp"

// Standard libs
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <map>
#include <string>
#include <vector>
#include <mpi.h>

// SUNDIALs libs
#include <cvode/cvode.h>
#include <cvode/cvode_dense.h>
#include <nvector/nvector_serial.h>

// POSIX libs
#include <dirent.h>
#include <magic.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

using std::copy;
using std::map;
using std::multimap;
using std::pair;
using std::string;
using std::vector;


void hrv_bench_optimal( int argc, char **argv, param_t *p_param )
{
  // CellML object pointer
  patch_clamp *p_cell;

  // Supporting variables
  char buffer[150] = {'\0'};
  char *pch;
  int idx, jdx, sample_id;
  int print_freq;
  string foldername;
  vector<double> temp_data;

  // HRV-related variables
  FILE *fp_hrv;
  map< string, vector<double> > hrv_data;

  // CiPA related variables
  cipa_t temp_result;
  double inal_auc_ctl, ical_auc_ctl, inal_auc_drug, ical_auc_drug, steepest_time, qnet, qinward, inal_auc, ical_auc;
  double hill_data_arr[14];
  FILE *fp_vm = NULL, *fp_ca = NULL;
  FILE *fp_hill, *fp_qni, *fp_debug;
  int pace, steepest_pace;
  vector< vector<double> > hill_data;  
  vector<double> concs;

  // AP related variables
  bool is_vm_peak_found;
  bool is_eligible_AP;
  double t_depol, vm_repol30, vm_repol90, vm_repol50;
  FILE *fp_ap_profile;

  // Ca related variables
  double ca_amp50, ca_amp90;
  double cad50_prev, cad90_prev, cad50_curr, cad90_curr;
  double t_ca_peak;
  FILE *fp_ca_profile, *fp_ca_debug;
  multimap<double, double> cai_data;

  // POSIX variable
  DIR *dir;
  struct magic_set *magic_cookie;
  struct stat st = {0};
  struct dirent *dent;
  
  // SUNDIALs variables
  bool cvode_firsttime;
  int cvode_retval, iout, imax;
  double tnext, tcurr = 0.0;
  void *cvode_mem;
  N_Vector states_vec;


  // libmagic setup for detecting filetype
  magic_cookie = magic_open(MAGIC_MIME|MAGIC_CHECK);
  magic_load(magic_cookie,NULL);

  // detect whether the hrv_files is a folder or zip
  string filetype( magic_file( magic_cookie, p_param->hrv_files ) );
  if( filetype.find("zip") != std::string::npos ){
    int zip_retval;
    printf("The HRV input is a zip file!!!\n");
    string zipname( p_param->hrv_files );
    foldername = zipname.substr( 0, zipname.size()-4 );
    zip_retval = extract_zip( p_param->hrv_files, foldername.c_str());
    if( zip_retval < 0 ){
      printf("Problem occured when extracting zipfile %s!!\n", p_param->hrv_files);
      return;
    }
    else{
      printf("Zipfile %s successfully extracted to folder %s!!\n", p_param->hrv_files, foldername.c_str());
    }
  }
  else{
    printf("The HRV input is not a zip file!!! Assume it as a folder....\n");
    foldername = p_param->hrv_files;
  }

  magic_close(magic_cookie);

  // Open the folder and process the data inside all of the files.
  dir = opendir(foldername.c_str());
  if(dir == NULL){
    fprintf(stderr, "Error opening the directory %s\n", foldername.c_str());
    return;
  }
  while( (dent = readdir(dir)) != NULL){ // begin hrv_files loop
    string hrv_filename(dent->d_name);
    if( hrv_filename.size() > 6 && hrv_filename.substr(hrv_filename.size()-3, 3) == "csv" )
    {
      // Open file
      snprintf( buffer, sizeof(buffer), "%s/%s", foldername.c_str(), dent->d_name );
      printf("File name: %s\t\t", buffer);
      fp_hrv = fopen(buffer, "r" );
      if( fp_hrv == NULL ) {fprintf(stderr, "Cannot open file %s\n", buffer ); continue;}

      // Fill hrv_data vector with the data inside the file
      temp_data.clear();
      while( fgets(buffer, sizeof(buffer), fp_hrv) != NULL ) temp_data.push_back(strtod(buffer, NULL));
      temp_data.push_back(1000.);
      printf("Size: %d\n", temp_data.size());

      fclose(fp_hrv);
      hrv_data.insert( std::pair< string, vector<double> > ( hrv_filename.substr(0, hrv_filename.size()-4), temp_data) );
    }
  } // end hrv_files loop
  closedir(dir);

#ifdef DEBUG
  for( std::map< string, vector<double> >::iterator itrmap = hrv_data.begin(); itrmap != hrv_data.end() ; itrmap++ )
  { // begin map debug
    printf( "File Name: %s\n", (itrmap->first).c_str() );
    printf( "Data: \n{\n" );
    for( std::vector<double>::iterator itrvec = (itrmap->second).begin(); itrvec != (itrmap->second).end(); itrvec++ )
    {
      printf("  %lf\n", *itrvec);
    }
    printf( "}\n\n" );
  } // end map debug
#endif

  // setup concentration vector
  concs = get_vector_from_string( std::string(p_param->concs) );
  if( concs.size() > 4+1 ){  // +1 because we include concentration 0
    fprintf(stderr, "Only 4 concentrations is allowed in this program! Exiting...\n");
    return;
  }
#ifdef DEBUG
  printf("Concentration: ");
  for( std::vector<double>::iterator itrvec = concs.begin(); itrvec != concs.end(); itrvec++ )
  {
    printf("%lf ", *itrvec);
  }
  printf("\n\n");
#endif

  // make folders from the concentration
  for( idx = 0; idx < concs.size(); idx++ ){
    snprintf( buffer, sizeof(buffer), "result/%.2lf", concs[idx] );
    if(mympi::rank == 0 && stat("result/", &st) == 0) mkdir(buffer, 0775);
  }

  // setup hill data vector
  fp_hill = fopen(p_param->hill_file, "r");
  if (fp_hill == NULL) {
      fprintf(stderr, "Cannot open Hill file %s\n", p_param->hill_file);
      return;
  }
  fgets(buffer, sizeof(buffer), fp_hill); // skip header
  temp_data.clear();
  while( fgets(buffer, sizeof(buffer), fp_hill) != NULL )
  { // begin data loop

    temp_data = get_vector_from_string( std::string(buffer) );
    hill_data.push_back( temp_data );

  } // end data loop
  fclose(fp_hill);
  if( hill_data.size() > 2000 ){  // sample cannot exceed 2000
    fprintf(stderr, "Hill samples are more than 2000! Exiting...\n");
    return;
  }

#ifdef DEBUG
  for( std::vector< vector<double> >::iterator itr = hill_data.begin(); itr != hill_data.end() ; itr++ )
  { 
    printf( "Hill Data: \n{\n" );
    for( std::vector<double>::iterator itrvec = itr->begin(); itrvec != itr->end(); itrvec++ )
    {
      printf("  %lf\n", *itrvec);
    }
    printf( "}\n\n" );
  } 
#endif

  // variable initialization
  p_cell = new Ohara_Rudy_2011();
  cvode_mem = CVodeCreate( CV_BDF,CV_NEWTON );
  cvode_firsttime = true;
  print_freq = (1. / p_param->dt) * p_param->dt_write;

  // main loop
  double begin = MPI_Wtime();
  for( sample_id = 0; sample_id < hill_data.size(); sample_id++ ){ // begin hill_data loop

    // convert one hill_data vector into an array
    copy( hill_data[sample_id].begin(), hill_data[sample_id].end(), hill_data_arr );

    for( idx = 0; idx < concs.size(); idx++ ){ // begin concentration loop
      printf("Concentration: %lf\n", concs[idx]);

      for( std::map< string, vector<double> >::iterator itrmap = hrv_data.begin(); 
            itrmap != hrv_data.end() ; itrmap++ )
      { // begin HRV file loop

        // init cell model
        p_cell->initConsts(p_param->celltype, concs[idx], hill_data_arr, true);
        p_cell->CONSTANTS[amp] *= p_param->stim_amp;

        // setup CVode setting
        tcurr = 0.0;
        tnext = p_param->dt;
        states_vec = N_VMake_Serial( p_cell->states_size, p_cell->STATES );
        if( cvode_firsttime == true ) {
          CVodeInit( cvode_mem, rhs_fn, tcurr, states_vec );
          CVodeSetUserData( cvode_mem, p_cell );
          CVodeSStolerances( cvode_mem, 1.0e-6, 1.0e-6 );
          CVodeSetMaxStep( cvode_mem, p_param->dt );
          CVDense( cvode_mem, p_cell->states_size );
          cvode_firsttime = false;
        }
        else CVodeReInit( cvode_mem, tcurr, states_vec );

        if( p_param->is_print_graph == 1 ){
          snprintf(buffer, sizeof(buffer), "./result/%.2lf/%s_%s_%.2lf_vmcheck.plt", concs[idx], p_param->drug_name, (itrmap->first).c_str(), concs[idx]);
          fp_vm = fopen( buffer, "w");
          snprintf(buffer, sizeof(buffer), "./result/%.2lf/%s_%s_%.2lf_cai.plt", concs[idx], p_param->drug_name, (itrmap->first).c_str(), concs[idx]);
          fp_ca = fopen( buffer, "w");
        }
        snprintf(buffer, sizeof(buffer), "./result/%.2lf/%s_%s_%.2lf_qni.plt", concs[idx], p_param->drug_name, (itrmap->first).c_str(), concs[idx]);
        fp_qni = fopen( buffer, "w");
        snprintf(buffer, sizeof(buffer), "./result/%.2lf/%s_%s_%.2lf_ap_profile.plt", concs[idx], p_param->drug_name, (itrmap->first).c_str(), concs[idx]);
        fp_ap_profile = fopen( buffer, "w");
        snprintf(buffer, sizeof(buffer), "./result/%.2lf/%s_%s_%.2lf_ca_profile.plt", concs[idx], p_param->drug_name, (itrmap->first).c_str(), concs[idx]);
        fp_ca_profile = fopen( buffer, "w");

        // Header Printing
        if( p_param->is_print_graph == 1 ){
          fprintf(fp_vm, "%s %s\n", "TIME", "Vm");
          fprintf(fp_ca, "%s %s\n", "TIME", "Cai");
        } 
        fprintf(fp_ap_profile, "%s %s %s %s %s %s %s\n", "Pace", "DVm/Dt_Max", "Vm_Peak", "Vm_Dia", "APD50", "APD90", "APDTri");
        fprintf( fp_ca_profile, "%s %s %s %s %s %s\n", "Pace", "Ca_Peak", "Ca_Diastole", "CaD90", "CaD50","Catri");
        fprintf(fp_qni, "%s %s %s\n", "Pace", "Qnet", "Qinward");

        // Initialization when entering new HRV data
        pace = 0;

        for( std::vector<double>::iterator itrvec = (itrmap->second).begin(); itrvec != (itrmap->second).end(); itrvec++ ){ // begin BCL loop

          // set the BCL value
          p_cell->CONSTANTS[BCL] = *itrvec;
          if(p_cell->CONSTANTS[BCL] <= 0.){printf("%s IS WRONG BCL VALUE!!!\n", buffer);continue;}
          //printf("BCL: %.2lf Line: %d\n", p_cell->CONSTANTS[BCL], pace+1);

          // mandatory initializations between paces (without this, expect some problems will happen)
          // THIS IS IMPORTANT!!!! PLEASE DON'T REMOVE THESE!!!
          qnet = inal_auc = ical_auc = 0.0;
          iout = 0;
          imax = p_cell->CONSTANTS[BCL] / p_param->dt;
          is_vm_peak_found = false;
          is_eligible_AP = false;
          t_ca_peak = tcurr;
          if(pace > 0){
            tcurr -= p_param->dt;
            tnext -= p_param->dt;
          }
          temp_result.init( p_cell->STATES[V], p_cell->STATES[cai] );

          while( iout <= imax ){ // begin one pacing loop
            // solve ODE
            cvode_retval = CVode( cvode_mem, tnext, states_vec, &tcurr, CV_NORMAL  );
            // get the RATES array
            p_cell->computeRates(tnext, p_cell->CONSTANTS, p_cell->RATES, N_VGetArrayPointer_Serial(states_vec), p_cell->ALGEBRAIC);
            if( cvode_retval == CV_SUCCESS ){
              iout++;
              tnext += p_param->dt;
            }
            else{
              printf("CVode calculation error at file %s line %d BCL %lf concentration %.2lf msec!!!\n", dent->d_name, pace+1, p_cell->CONSTANTS[BCL], concs[idx]);
              break;
            }
            // calculation for CiPA parameters
            qnet += (p_cell->ALGEBRAIC[INaL]+p_cell->ALGEBRAIC[ICaL]+p_cell->ALGEBRAIC[Ito]+p_cell->ALGEBRAIC[IKr]+p_cell->ALGEBRAIC[IKs]+p_cell->ALGEBRAIC[IK1])*p_param->dt;
            inal_auc += (p_cell->ALGEBRAIC[INaL])*p_param->dt;
            ical_auc += (p_cell->ALGEBRAIC[ICaL])*p_param->dt;

            // save output for each pace
            if( iout % print_freq == 0 ) {
              temp_result.cai_data.insert( std::pair<double, double> (tcurr, p_cell->STATES[cai]) );
              temp_result.vm_data.insert( std::pair<double, double> (tcurr, p_cell->STATES[V]) );
            }

            // get maximum dvmdt
            if( temp_result.dvmdt_max < p_cell->RATES[V] ) temp_result.dvmdt_max = p_cell->RATES[V];

            // get peak vm at the 6 seconds after stimulus start time
            if( iout % (int)((p_cell->CONSTANTS[stim_start]+6.) / p_param->dt) == 0 && is_vm_peak_found == false ) {
              temp_result.vm_peak = p_cell->STATES[V];
              t_depol = tcurr;
              // only finding the vm_repol30, vm_repol_50 and vm_repol90 if the peak Vm more than 0 mV
              if( temp_result.vm_peak > 0. ){
                vm_repol30 = temp_result.vm_peak - (0.3 * (temp_result.vm_peak - temp_result.vm_valley));
                vm_repol50 = temp_result.vm_peak - (0.5 * (temp_result.vm_peak - temp_result.vm_valley));
                vm_repol90 = temp_result.vm_peak - (0.9 * (temp_result.vm_peak - temp_result.vm_valley));
                is_eligible_AP = true;
              }
              is_vm_peak_found = true;
            }

            // if the pace is eligible AP, find the rest of paramters such as APD90, APD50, CAD90, CAD50
            if( is_eligible_AP == true ){
              // get the steepest dvmdt_repol around vm_repol30 and vm_repol90 if AP shape is eligible
              if( vm_repol30 >= p_cell->STATES[V] && p_cell->STATES[V] >= vm_repol90 && temp_result.dvmdt_repol < p_cell->RATES[V] ){
                temp_result.dvmdt_repol = p_cell->RATES[V];
                steepest_time = tcurr;
              }
              // get the APD50
              if( vm_repol50 > p_cell->STATES[V] && p_cell->STATES[V] > vm_repol50-p_param->dt ){
                temp_result.apd50 = tcurr - t_depol;
              }
              // get the APD90
              if( vm_repol90 > p_cell->STATES[V] && p_cell->STATES[V] > vm_repol90-(p_param->dt*2) ){
                temp_result.apd90 = tcurr - t_depol;
              }
            }

            // finding peak calcium
            if( temp_result.ca_peak < p_cell->STATES[cai] ){
              temp_result.ca_peak = p_cell->STATES[cai];
              ca_amp50 = temp_result.ca_peak - ( 0.5 * (temp_result.ca_peak - temp_result.ca_valley) );
              ca_amp90 = temp_result.ca_peak - ( 0.9 * (temp_result.ca_peak - temp_result.ca_valley) );
              t_ca_peak = tcurr;
            }

          } // end one pacing loop

          temp_result.vm_dia = p_cell->STATES[V];
          temp_result.ca_dia = p_cell->STATES[cai];
          for(std::multimap<double, double>::iterator itrmap = cai_data.begin(); itrmap != cai_data.end() ; itrmap++ ){
            // before the peak calcium
            if( itrmap->first < t_ca_peak ){
              if( itrmap->second < ca_amp50 ) cad50_prev = itrmap->first;
              if( itrmap->second < ca_amp90 ) cad90_prev = itrmap->first;
            }
            // after the peak calcium
            else{
              if( itrmap->second > ca_amp50 ) cad50_curr = itrmap->first;
              if( itrmap->second > ca_amp90 ) cad90_curr = itrmap->first;
            }
          }
          temp_result.cad50 = cad50_curr - cad50_prev;
          temp_result.cad90 = cad90_curr - cad90_prev;
          temp_result.qnet = qnet;
          temp_result.inal_auc = inal_auc;
          temp_result.ical_auc = ical_auc;

          // writing final results
          if( cvode_retval == CV_SUCCESS ){

            if(p_param->is_print_graph == 1 ) {
              for(std::multimap<double, double>::iterator itrmap = temp_result.cai_data.begin(); itrmap != temp_result.cai_data.end() ; itrmap++ ){
                fprintf(fp_ca, "%lf %lf\n", itrmap->first, itrmap->second);
              }
              for(std::multimap<double, double>::iterator itrmap = temp_result.vm_data.begin(); itrmap != temp_result.vm_data.end() ; itrmap++ ){
                fprintf(fp_vm, "%lf %lf\n", itrmap->first, itrmap->second);
              }
            }

            fprintf(fp_ap_profile, "%d %lf %lf %lf %lf %lf %lf\n", pace, temp_result.dvmdt_max, temp_result.vm_peak, temp_result.vm_dia, temp_result.apd50, temp_result.apd90, temp_result.apd90-temp_result.apd50);
            fprintf( fp_ca_profile, "%d %lf %lf %lf %lf %lf\n", pace, temp_result.ca_peak, temp_result.ca_dia, temp_result.cad90, temp_result.cad50, temp_result.cad90-temp_result.cad50 );
            if((int)ceil(concs[idx]) == 0){
              inal_auc_ctl = temp_result.inal_auc;
              ical_auc_ctl = temp_result.ical_auc;
              qinward = 0.0;
            }
            else{
              inal_auc_drug = temp_result.inal_auc;
              ical_auc_drug = temp_result.ical_auc;
              qinward = ((inal_auc_drug/inal_auc_ctl) + (ical_auc_drug/ical_auc_ctl)) * 0.5;
            }
            fprintf(fp_qni, "%d %lf %lf\n", pace, temp_result.qnet, qinward  );
          }
          else{
            fprintf(fp_ap_profile, "%d ERR ERR ERR ERR ERR ERR\n", pace );
            fprintf(fp_ca_profile, "%d ERR ERR ERR ERR ERR\n", pace );
            fprintf(fp_qni, "%d ERR ERR\n", pace );
          }


          // increase pace and accumulated BCL
          pace++;
          p_cell->CONSTANTS[cumm_BCL] += p_cell->CONSTANTS[BCL];

        } // end BCL loop

        printf("File %s with concentration %.4lf and total pace %d using sample %d finished!!!\n", (itrmap->first).c_str(), concs[idx], (itrmap->second).size(), sample_id);

        // memory cleanup
        N_VDestroy_Serial(states_vec);
        if( p_param->is_print_graph == 1 ){
          fclose(fp_vm);
          fclose(fp_ca);
        }
        fclose(fp_qni);
        fclose(fp_ap_profile);
        fclose(fp_ca_profile);

      } // end HRV file debug

    } // end concentration loop

  } // end hill_data loop
  double end = MPI_Wtime();

  // Memory cleanup
  CVodeFree(&cvode_mem);
  free(p_param);
  delete p_cell;

  // setup timestamp
  char buff_time[15] = {'\0'};
  FILE *fp_perf_log;
  struct tm* tm_info;
  time_t timer;
  timer = time(NULL);
  tm_info = localtime(&timer);
  strftime( buff_time, 15, "%Y%m%d%H%M%S", tm_info );

  // zip result folder
  snprintf( buffer, sizeof(buffer), "result_%s.zip", buff_time );
  create_zip(buffer, "result");

  // delete the content inside result folder
  remove_dir_content("result");

  // copy the zip into the result folder
  char buff_new[150] = {'\0'};
  snprintf( buffer, sizeof(buffer), "result_%s.zip", buff_time );
  snprintf( buff_new, sizeof(buff_new), "result/result_%s.zip", buff_time );
  if( rename( buffer, buff_new ) == 0 ){
    printf("Zip file %s moved succesfully to %s!!\n", buffer, buff_new);
  }
  else{
    printf("Error moving file %s to %s!!\n", buffer, buff_new);
    return;
  }

  // make a timestamp log
  snprintf( buffer, sizeof(buffer), "result/performance_%s.log", buff_time );
  fp_perf_log = fopen( buffer, "w" );
  fprintf( fp_perf_log, "File_Count CPU Concentration Time(Minutes)\n" );
  fprintf( fp_perf_log, "%d %d %d %.5lf \n", hrv_data.size(), mympi::size, concs.size(), (end-begin)/60. );
  fclose( fp_perf_log );
  printf("Computation Time: %lf minutes\n", (end-begin)/60.);

}
