#ifndef COMMONS_HPP
#define COMMONS_HPP

#include "param.hpp"

#include <cstdio>
#include <map>
#include <string>
#include <vector>

#include <nvector/nvector_serial.h>

using std::multimap;
using std::string;
using std::vector;

// struct for storing the result of CiPA simulation
struct cipa_t{
  double qnet;
  double inal_auc;
  double ical_auc;
  double dvmdt_repol;
  double dvmdt_max;
  double vm_peak;
  double vm_valley;
  double vm_dia;
  double apd90;
  double apd50;
  double apd_tri;
  double ca_peak;
  double ca_valley;
  double ca_dia;
  double cad90;
  double cad50;
  double cad_tri;
  int bcl;
  multimap<double, double> vm_data;
  multimap<double, double> cai_data;
  multimap<double, string> ires_data;

  cipa_t()
  {

  }

  void copy(const cipa_t &source)
  {
    qnet = source.qnet;
    inal_auc = source.inal_auc;
    ical_auc = source.ical_auc;
    dvmdt_repol = source.dvmdt_repol;
    dvmdt_max = source.dvmdt_max;
    vm_peak = source.vm_peak;
    vm_valley = source.vm_valley;
    vm_dia = source.vm_dia;
    apd90 = source.apd90;
    apd50 = source.apd50;
    ca_peak = source.ca_peak;
    ca_valley = source.ca_valley;
    ca_dia = source.ca_dia;
    cad90 = source.cad90;
    cad50 = source.cad50;
    bcl = source.bcl;
    vm_data.clear();
    cai_data.clear();
    ires_data.clear();
    vm_data.insert( (source.vm_data).begin(), (source.vm_data).end() );
    cai_data.insert( (source.cai_data).begin(), (source.cai_data).end() );  
    ires_data.insert( (source.ires_data).begin(), (source.ires_data).end() );  
  }

  cipa_t( const cipa_t &source )
  {
    copy(source);
  }


  cipa_t& operator=(const cipa_t & source)
  {
    if( this != &source )
    {
      copy(source);
    }

    return *this;
  }

  void init(const double vm_val, const double ca_val)
  {
    qnet = 0.;
    inal_auc = 0.;
    ical_auc = 0.;
    dvmdt_repol = -999;
    dvmdt_max = -999;
    vm_peak = -999;
    vm_valley = vm_val;
    vm_dia = -999;
    apd90 = 0.;
    apd50 = 0.;
    ca_peak = -999;
    ca_valley = ca_val;
    ca_dia = -999;
    cad90 = 0.;
    cad50 = 0.;
    vm_data.clear();
    cai_data.clear();
    ires_data.clear();
  }

};

struct anm_result_t{
  double anm_aocl;
  double anm_alternant;
  double aocl;
  double mean_apd;
  double cl_alternant;
  double mean_apd_alternant;
};

int rhs_fn( realtype t, N_Vector y, N_Vector ydot, void *user_data );
int is_file_exists( char *filename );
void set_default_values( param_t *p_param );
void edison_assign_params_single(int argc, char *args[], param_t *p_param);

// convert comma-delimited string into a vector
vector<double> get_vector_from_string( string delim_str );

// Main interface for zipping a folder
void create_zip(const char *zip_filename, const char *target_folder);

// Main interface for extracting a zipfile
int extract_zip(const char *zip_filename, const char *target_folder);

// Zip folder recursively
// Source:
// https://github.com/kuba--/zip
void zip_walk(struct zip_t *zip, const char *path);

// delete folder that have contents
// Source:
// https://stackoverflow.com/a/54956690/981481
void remove_dir_content(const char *path);

// convert decimal number into a base number
int convert_base_n(int number,int base);

// get the ANM from APD vector
double get_ANM( vector<double> vec );

// get the average from a vector within start_idx and end_idx
double average( vector<double> vec, int start_idx, int end_idx );


void create_apdr_output(int num_pace1, int number_of_lines);


#endif
