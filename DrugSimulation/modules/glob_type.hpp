#ifndef GLOB_TYPE_HPP
#define GLOB_TYPE_HPP

#include <cstddef>
#include <vector>


// global variable for MPI.
struct mympi
{
  static char host_name[255];
  static int host_name_len;
  static int rank;
  static int size;
};

// data structure for IC50
typedef struct row_data { double data[14]; } row_data;
typedef std::vector< row_data > drug_t;

// data structure for hERG bootstrap
typedef struct row_data2 { double data[6]; } row_data2;
typedef std::vector< row_data2 > boot_t;

// data structure to store
// ICaL/INaL control value
// for calculating qinward
// control means 0 concentration
// otherwise, drugs
typedef struct{
  double ical_auc_control;
  double inal_auc_control;
  double ical_auc_drug;
  double inal_auc_drug;
} qinward_t;

// holder for curl fetch
typedef struct tuple {
    char *payload;
    size_t size;
}curl_fetch_t;

// structure of lock data
typedef struct file_content{
	char ipaddress[20];
	unsigned long long timestamp;
}lock_t;

#endif
