// Code located in this project directory
#include "hrv_bench_combined.hpp"
#include "commons.hpp"
#include "globals.hpp"
#include "../cellmodels/Ohara_Rudy_2011.hpp"

// Standard libs
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <map>
#include <string>
#include <sstream>
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
using std::getline;
using std::istringstream;
using std::map;
using std::multimap;
using std::pair;
using std::string;
using std::vector;

void hrv_bench_combined( int argc, char **argv, param_t *p_param )
{
  // CellML object pointer
  patch_clamp *p_cell;

  // Supporting variables
  char buffer[150] = {'\0'};
  char *pch;
  int idx, jdx, group_id, sample_id, sim_counter;
  int print_freq;
  istringstream iss;
  pair< string, double > pr;
  string temp;
  vector<double> temp_data;

  // HRV-related variables
  FILE *fp_hrv;
  vector<double> hrv_data;

  // CiPA related variables
  cipa_t cipa_result, temp_result;
  double inal_auc_ctl, ical_auc_ctl, inal_auc_drug, ical_auc_drug, qnet, qinward, inal_auc, ical_auc;
  double hill_data_arr[14];
  FILE *fp_vm = NULL, *fp_ca = NULL;
  FILE *fp_hill, *fp_qni, *fp_debug, *fp_cmax;
  int pace, steepest_pace;
  multimap< string, vector< vector<double> > > hill_data;  
  multimap< string, double > cmax_data;  
  vector<double> conc_scalars;

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
  struct dirent *dent;
  struct stat st = {0};
  
  // SUNDIALs variables
  bool cvode_firsttime;
  int cvode_retval, iout, imax;
  double tnext, tcurr = 0.0;
  void *cvode_mem;
  N_Vector states_vec;


  // Open the file and process the data inside it.
  fp_hrv = fopen(p_param->hrv_files, "r" );
  if( fp_hrv == NULL ) {fprintf(stderr, "Cannot open file %s\n", buffer ); return;}
  while( fgets(buffer, sizeof(buffer), fp_hrv) != NULL ) hrv_data.push_back(strtod(buffer, NULL));
  hrv_data.push_back(1000.);
  fclose(fp_hrv);

  // setup concentration vector
  conc_scalars = get_vector_from_string( std::string(p_param->concs) );

  // read cmax data
  fp_cmax = fopen("test_drug_data/cmax_tables.csv", "r");
  fgets(buffer, sizeof(buffer), fp_cmax); // skip header
  while( fgets(buffer, sizeof(buffer), fp_cmax) != NULL ) {
    temp = buffer;
    iss.str( temp );
    getline( iss, temp, ',' );
    pr.first = temp;
    getline( iss, temp, ',' );
    pr.second = strtod(temp.c_str(), NULL);
    cmax_data.insert(pr);
    iss.clear();
  }
  fclose(fp_cmax);

  for( std::multimap< string, double >::iterator itrmap = cmax_data.begin(); itrmap != cmax_data.end() ; itrmap++ ) printf("%s %lf\n", (itrmap->first).c_str(), itrmap->second);

  // Open the folder and process the data inside all of the files.
  dir = opendir( "test_drug_data" );
  assert(dir != NULL && "directory is null!!");
  
  while( (dent = readdir(dir)) != NULL){ // begin hrv_files loop

    string hill_filename(dent->d_name);
    if( hill_filename.size() > 6 && hill_filename.substr(hill_filename.size()-3, 3) == "csv" && hill_filename.find("IC50") != std::string::npos )
    {
      // Open file
      snprintf( buffer, sizeof(buffer), "%s/%s", "test_drug_data", dent->d_name );
      printf("File name: %s\t\t\n", buffer);

      fp_hill = fopen( "buffer", "r" );
      
      fclose(fp_hill);

    }

  } // end hrv_files loop  
  
  closedir(dir);

}
