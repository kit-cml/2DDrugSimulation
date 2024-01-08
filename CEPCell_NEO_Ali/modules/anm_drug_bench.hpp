#ifndef ANM_BENCH_HPP
#define ANM_BENCH_HPP

#include "commons.hpp"
#include "../cellmodels/patch_clamp.hpp"

anm_result_t anm_drug_bench(int argc, char **argv, param_t *p_param, patch_clamp* p_cell, const int sample_id, const double concentration);

#endif
