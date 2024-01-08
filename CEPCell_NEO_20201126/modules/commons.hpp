#ifndef COMMONS_HPP
#define COMMONS_HPP

#include "param.hpp"
#include <nvector/nvector_serial.h>

typedef struct result{
  double qnet_prev;
  double qnet_curr;
  double INaL_auc_prev;
  double ICaL_auc_prev;
  double INaL_auc_curr;
  double ICaL_auc_curr;
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
}cipa_t;


class mympi
{
 public:
    static int rank, size;
};


int rhs_fn( realtype t, N_Vector y, N_Vector ydot, void *user_data );
void set_default_values( param_t *p_param );
void edison_assign_params_single(int argc, char *args[], param_t *p_param);


#endif
