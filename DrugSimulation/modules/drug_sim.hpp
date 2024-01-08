#ifndef DRUG_SIM_HPP
#define DRUG_SIM_HPP

#include "cipa_t.hpp"
#include "param.hpp"
#include "glob_type.hpp"

// main simulation function for drug simulation
// return TRUE if EAD happened, otherwise false
bool do_drug_sim(const double conc, row_data ic50, row_data2 herg_pars,
const param_t* p_param, const unsigned short sample_id, const unsigned short group_id, qinward_t *p_qin);

#endif
