#ifndef HELPER_HPP
#define HELPER_HPP

#include "../commons/usermacro.hpp"
#include "param.hpp"

// if you just make another cell models, add those to these line
#if defined CRN1998
  #include "../cellmodels/courtemanche_ramirez_nattel_1998.hpp"
#elif defined HHUXLEY1952
  #include "../cellmodels/hodgkin_huxley_squid_axon_model_1952.hpp"
#elif defined OHARA_RUDY2011
  #include "../cellmodels/Ohara_Rudy_2011.hpp"
#elif defined ORUDY_CIPA2017
  #include "../cellmodels/ohara_rudy_cipa_v1_2017.hpp"
#elif defined TN2004ENDO
  #include "../cellmodels/tentusscher_noble_noble_panfilov_2004_a.hpp"
#elif defined TN2004EPI
  #include "../cellmodels/tentusscher_noble_noble_panfilov_2004_b.hpp"
#elif defined TN2004M
  #include "../cellmodels/tentusscher_noble_noble_panfilov_2004_c.hpp"
#elif defined TN2006ENDO
  #include "../cellmodels/ten_tusscher_model_2006_IK1Ko_endo_units.hpp"
#elif defined TN2006EPI
  #include "../cellmodels/ten_tusscher_model_2006_IK1Ko_epi_units.hpp"
#elif defined TN2006M
  #include "../cellmodels/ten_tusscher_model_2006_IK1Ko_M_units.hpp"
#endif
// end of cell model definition

// mutation definition
#if defined CRN1998
  #include "../mutations/hERG.hpp"
  #include "../mutations/S140G.hpp"
  #include "../mutations/Fibrosis.hpp"
#elif defined TN2004ENDO || defined TN2004M || defined TN2004EPI
  #include "../mutations/A1656D.hpp"
#elif defined OHARA_RUDY2011 || defined ORUDY_CIPA2017
  #include "../mutations/G229D.hpp"
#endif
// end of cell model definition


double get_cmax( char *drug_name );
void set_default_values( param_t *p_param );
void assign_params_single(int argc, char **argv, param_t *p_param);
void edison_assign_params_single(int argc, char *args[], param_t *p_param);
int load_last_state(CellML *p_cell);
int save_last_state(CellML *p_cell, double t_laststate);

CellML* init_cell();
CellML* init_mutation(const char variant_ID[], double scale, CellML *p_cell);

#endif
