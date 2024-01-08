#include "glob_funct.hpp"
#include "../libs/zip.h"

#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>

// to make it more "portable" between OSes.
#if defined _WIN32
  #include <direct.h>
  #define snprintf _snprintf
  #define vsnprintf _vsnprintf
  #define strcasecmp _stricmp
  #define strncasecmp _strnicmp
#else
  #include <dirent.h>
#endif
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <curl/curl.h>
#include <json-c/json.h>

void mpi_printf(unsigned short node_id, const char *fmt, ...)
{
#ifndef _WIN32
  if(mympi::rank == node_id){
    va_list args;
    va_start(args, fmt);
    vprintf(fmt, args);
    va_end(args);
  }
#else
  va_list args;
  va_start(args, fmt);
  vprintf(fmt, args);
  va_end(args);	
#endif
}

void mpi_fprintf(unsigned short node_id, FILE *stream, const char *fmt, ...)
{
#ifndef _WIN32
  if(mympi::rank == node_id){
    va_list args;
    va_start(args, fmt);
    vfprintf(stream, fmt, args);
    va_end(args);
  }
#else
  va_list args;
  va_start(args, fmt);
  vfprintf(stream, fmt, args);
  va_end(args);	
#endif
}

void edison_assign_params(int argc, char *argv[], param_t *p_param)
{
  bool is_default;
  char buffer[100];
  char key[100];
  char value[100];
  char file_name[150];
  FILE *fp_inputdeck;

  // parameters from arguments
  for (int idx = 1; idx < argc; idx += 2) {
    if (!strcasecmp(argv[idx], "-input_deck"))
      strncpy(file_name, argv[idx + 1], sizeof(file_name));
    else if (!strcasecmp(argv[idx], "-hill_file"))
      strncpy(p_param->hill_file, argv[idx + 1], sizeof(p_param->hill_file));
	else if (!strcasecmp(argv[idx], "-herg_file"))
      strncpy(p_param->herg_file, argv[idx + 1], sizeof(p_param->herg_file));
  }  

  is_default = false;
  fp_inputdeck = fopen( file_name, "r");
  if(fp_inputdeck == NULL){
    fprintf(stderr, "Cannot open input deck file %s!!!\nUse default value as the failsafe.\n", file_name);
    is_default = true;
  }

  // read input_deck line by line
  // and store each line to the buffer
  while ( is_default == false && fgets( buffer, 100, fp_inputdeck ) != NULL ) {
    sscanf( buffer, "%s %*s %s", key, value );
    if (strcasecmp(key, "Simulation_Mode") == 0) {
      p_param->simulation_mode = strtod( value, NULL );
    }
    else if (strcasecmp(key, "Celltype") == 0) {
      p_param->celltype = strtod( value, NULL );
    }
    else if (strcasecmp(key, "Is_Dutta") == 0) {
      p_param->is_dutta = strtol( value, NULL, 10 );
    }
    else if (strcasecmp(key, "Is_Print_Graph") == 0) {
      p_param->is_print_graph = strtol( value, NULL, 10 );
    }
    else if (strcasecmp(key, "Is_Using_Output") == 0) {
      p_param->is_using_output = strtol( value, NULL, 10 );
    }
    else if (strcasecmp(key, "Basic_Cycle_Length") == 0) {
      p_param->bcl = strtod( value, NULL );
    }
    else if (strcasecmp(key, "Number_of_Pacing") == 0) {
      p_param->pace_max = strtol( value, NULL, 10 );
    }
	else if (strcasecmp(key, "Last_Drug_Check_Pace") == 0) {
      p_param->last_drug_check_pace = strtol( value, NULL, 10 );
    }
    else if (strcasecmp(key, "Time_Step") == 0) {
      p_param->dt = strtod( value, NULL );
    }
    else if (strcasecmp(key, "Write_Step") == 0) {
      p_param->dt_write = strtod( value, NULL );
    }
    else if (strcasecmp(key, "Drug_Name") == 0) {
      strncpy( p_param->drug_name, value, sizeof(p_param->concs) );
    }
    else if (strcasecmp(key, "Inet_Vm_Threshold") == 0) {
      p_param->inet_vm_threshold = strtod( value, NULL );
    }
    else if (strcasecmp(key, "Concentrations") == 0) {
      strncpy( p_param->concs, value, sizeof(p_param->concs) );
    }

  }

  if( is_default == false ) fclose( fp_inputdeck );
}

#ifndef _WIN32
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

void zip_walk(struct zip_t *zip, const char *path) 
{
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
#endif

int make_directory(const char* dirname )
{
#if defined _WIN32
  return _mkdir(dirname);
#else
  return mkdir(dirname, 0775);
#endif	
}

int is_file_existed(const char* pathname)
{
#if defined _WIN32
  struct _stat buf;
  return _stat( pathname, &buf );
#else
  struct stat st = {0};
  return stat(pathname, &st);
#endif
}

CURLcode curl_fetch_url(CURL *curlptr, const char *url, curl_fetch_t *fetch)
{
	CURLcode curl_response;

	fetch->size = 0;
	fetch->payload = (char*)calloc(1, sizeof(fetch->payload));
	if(fetch->payload == NULL){
		fprintf(stderr, "ERROR: Failed to allocate payload in curl_fetch_url\n");
		return CURLE_FAILED_INIT;
	}

	/* Refer to the docs for those options
 	* https://curl.se/libcurl/c/curl_easy_setopt.html
  */
	curl_easy_setopt(curlptr, CURLOPT_URL, url);
	curl_easy_setopt(curlptr, CURLOPT_WRITEFUNCTION, curl_callback);
	curl_easy_setopt(curlptr, CURLOPT_WRITEDATA, (void*) fetch);
	curl_easy_setopt(curlptr, CURLOPT_USERAGENT, "libcurl-agent/1.0");
	curl_easy_setopt(curlptr, CURLOPT_TIMEOUT, 5);
	curl_easy_setopt(curlptr, CURLOPT_FOLLOWLOCATION, 1L);
	curl_easy_setopt(curlptr, CURLOPT_MAXREDIRS, 1);
	curl_easy_setopt(curlptr, CURLOPT_TIMEOUT, 5);

	curl_response = curl_easy_perform(curlptr);

	return curl_response;
}

size_t curl_callback (void *contents, size_t size, size_t nmemb, void *userptr)
{
	size_t realsize = size * nmemb;               /* calculate buffer size */
	curl_fetch_t *ptr = (curl_fetch_t*)userptr;   /* cast pointer to fetch struct */

	ptr->payload = (char *) realloc(ptr->payload, ptr->size + realsize + 1);
	if(ptr->payload == NULL){
		fprintf(stderr, "ERROR: Failed to expand buffer in curl_callback\n");
		return -1;
	}

	memcpy(&(ptr->payload[ptr->size]), contents, realsize);

	ptr->size += realsize;
	ptr->payload[ptr->size] = 0;


	return realsize;
}

int get_credentials(lock_t *lockptr)
{
	CURL *curl;
	CURLcode curl_response;

	json_object *json_master;                               /* json post body */
	json_object *json_ipaddress;                            /* json for IPAddress */      
	json_object *json_unixtime;                             /* json for UNIX Time */
	enum json_tokener_error jerr = json_tokener_success;    /* json parse error */

	curl_fetch_t cf;
  curl_fetch_t *curl_fetch = &cf;
  const char url[100] = "http://worldtimeapi.org/api/timezone/Asia/Seoul";
 
	curl = curl_easy_init();
	if(curl == NULL){
		fprintf(stderr, "ERROR: Cannot initialize curl object!!\n");
		return 1;
	}
	if(curl_fetch == NULL){
		fprintf(stderr, "ERROR: Cannot initialize fetch object!!\n");
		return 2;
	}
	json_master = json_object_new_object();
	if( json_master == NULL ){
		fprintf(stderr, "ERROR: Cannot initialize json_master object!!\n");
		return 3;
	}
	json_ipaddress = json_object_new_object();
	if( json_master == NULL ){
		fprintf(stderr, "ERROR: Cannot initialize json_ipaddress object!!\n");
		return 4;
	}
	json_unixtime = json_object_new_object();
	if( json_unixtime == NULL ){
		fprintf(stderr, "ERROR: Cannot initialize json_unixtime object!!\n");
		return 5;
	}


	curl_response = curl_fetch_url(curl, url, curl_fetch);
	if(curl_response != CURLE_OK){
		fprintf(stderr, "ERROR: curl_fetch_url() failed: %s\n",
			curl_easy_strerror(curl_response));
		return 6;
	}
	if( curl_fetch->payload == NULL ){
		fprintf(stderr, "ERROR: failed to populate payload");
		return 7;
	}

	//printf("CURL Fetch: %s\n", curl_fetch->payload);
	json_master = json_tokener_parse_verbose(curl_fetch->payload, &jerr);
	json_object_object_get_ex(json_master,"unixtime",&json_unixtime);
	json_object_object_get_ex(json_master,"client_ip",&json_ipaddress);
	
	const char *ip = json_object_to_json_string(json_ipaddress);
  strncpy(lockptr->ipaddress, ip, sizeof(lockptr->ipaddress));
	lockptr->timestamp = strtol(json_object_to_json_string(json_unixtime),NULL,10);


	/* debugging */
	//printf("Parsed JSON: %s\n", json_object_to_json_string(json_master));
	//printf("Parsed JSON IP Address: %s\n", lockptr->ipaddress);
	//printf("Parsed JSON UNIX Time: %llu\n",lockptr->timestamp);

	curl_easy_cleanup(curl);
	json_object_put(json_master);
	json_object_put(json_ipaddress);
	json_object_put(json_unixtime);
	free(curl_fetch->payload);
	return 0;

}

int read_lock_file(lock_t *lockptr)
{
  FILE *fptr;
	const char *LOCK_FILE_NAME = "cardio.lock";

  fptr = fopen(LOCK_FILE_NAME,"rb");
  if(fptr == NULL){
    fprintf(stderr, "ERROR: Cannot read lock file!!\n");
    return 1;
  }
  while( fread(lockptr, sizeof(lock_t), 1, fptr) == 1 )
    printf("IP Address:%s\nRegister Date:%llu\n", 
		lockptr->ipaddress, lockptr->timestamp);

  fclose(fptr);
  return 0;
}


bool is_trial_active(unsigned long long current_timestamp,unsigned long long first_timestamp, unsigned short *usage_days)
{
	const unsigned long long DAYS_TO_SECOND = 86400;
	const unsigned short TRIAL_MAX = 30;
	*usage_days = abs(current_timestamp-first_timestamp)/DAYS_TO_SECOND;
	if(*usage_days <= TRIAL_MAX) return true;
	else return false; 
}
