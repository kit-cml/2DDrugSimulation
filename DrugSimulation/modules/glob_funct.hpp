#ifndef GLOB_FUNCT_HPP
#define GLOB_FUNCT_HPP

#include <cstdio>
#include <curl/curl.h>

#include "glob_type.hpp"
#include "param.hpp"
#include "../cellmodels/cellmodel.hpp"

// custom printf for MPI
// to avoid duplicate printing.
// In Windows, this function will be same
// as your usual printf() or fprintf().
// Will be defined for the sake of portability.
void mpi_printf(unsigned short node_id, const char *fmt, ...);
void mpi_fprintf(unsigned short node_id, FILE *stream, const char *fmt, ...);

// parameter setup function
void edison_assign_params(int argc, char *argv[], param_t *p_param);

#ifndef _WIN32
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
#endif

// create a directory.
// supporting different OS.
int make_directory(const char* dirname );

// checking file availability
// supporting different OS.
int is_file_existed(const char* pathname);

/* Fetch and return URL body */
CURLcode curl_fetch_url(CURL *curlptr, const char *url, curl_fetch_t *fetch);

/* Callback function for WRITEFUNCTION option */
size_t curl_callback (void *contents, size_t size, size_t nmemb, void *userptr);

/* Function to get real-time data from internet */
/* modified from:
 *  * https://gist.github.com/leprechau/e6b8fef41a153218e1f4
 *   */
int get_credentials(lock_t *lockptr);

/* Read the contents of the lock file. */
/* This function shouldn't be exposed carelessly. */
int read_lock_file(lock_t *lockptr);

/* Return true if the trial license still active. */
/* Timestamps use UNIX time. */
bool is_trial_active(unsigned long long current_timestamp,unsigned long long first_timestamp, unsigned short *usage_days);

#endif
