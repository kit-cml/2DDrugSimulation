#include "population_orudy.hpp"
#include "rhs.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <cvode/cvode.h>
#include <cvode/cvode_dense.h>
#include <nvector/nvector_serial.h>

#if defined ORUDY2011_STATIC123456
int population_ORudy(int argc, char **argv, param_t *p_param)
{
  /* CellML object pointer */
  CellML *p_cell;
  /* time related variables */
  int print_freq;
  int err_code;
  int i;
  char file_name[100];
  FILE *fp_apd;
  FILE *fp_vm;
  FILE *fp_ca;
  FILE *fp_ires;

  /* CiPA related variables */
  FILE *fp_hill;
  FILE *fp_qnet;
  FILE *fp_qinward;
  FILE *fp_ap_profile;
  FILE *fp_ca_profile;
  double hill[14];
  double conc;
  char *pch;
  char concs_str[100];
  char hill_sample[128];
  char result_file[128];
  int sim_id;
  double qnet_prev;
  double qnet_curr;
  double INaL_auc_prev;
  double INaL_auc_curr;
  double ICaL_auc_prev;
  double ICaL_auc_curr;
  double INaL_auc_control;
  double ICaL_auc_control;
  double INaL_auc_drug;
  double ICaL_auc_drug;
  double qinward;

  /* Calcium related variable */
  double ca_prev;
  double ca_peak;
  double ca_dia;
  double ca_amp;
  double t1_cad90;
  double t2_cad90;
  double cad90;
  double t1_cad50;
  double t2_cad50;
  double cad50;
  bool is_ca90_inc;
  bool is_ca50_inc;
  bool is_cad90_found;
  bool is_cad50_found;


  /* AP related variables */
  bool is_ap_increasing;
  bool is_action_potential;
  bool is_apd_written;
  double apd90;
  double apd50;
  double max_dvmdt;
  double vm_peak;
  double vm_dia;
  double v_prev;
  double v_apd90;
  double v_apd50;
  double t_depolarize;
  double t_repolarize;

  /* SUNDIALs variables */
  int retval;
  double t_out, t;
  void *cvode_mem;
  N_Vector states_vec;
  int iout;
  int imax;


  int a,b,c,d,e,f,g,h;

  for(a = 20; a <= 100; a += 20){
    for(b = 20; a <= 100; b += 20){
      for(c = 20; a <= 100; c += 20){
        for(d = 20; a <= 100; d += 20){
          for(e = 20; a <= 100; e += 20){
            for(f = 20; a <= 100; f += 20){
              for(g = 20; a <= 100; g += 20){
                for(h = 20; a <= 100; h += 20){
                  /* Initialize p_cell based on the ionmodel value */
                  p_cell = init_cell( );
                  if( p_cell == NULL ){
                    fprintf(stderr, "Problem when initializing cellmodel\n");
                    return EXIT_FAILURE;
                  }

                }
              }
            }
          }
        }
      }
    }
  }

      
  delete p_cell;
  free( p_param );

  return 0;
}
#endif
