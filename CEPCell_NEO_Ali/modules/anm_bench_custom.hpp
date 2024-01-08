#ifndef ANM_BENCH_CUSTOM_HPP
#define ANM_BENCH_CUSTOM_HPP

#include "commons.hpp"
#include "../cellmodels/patch_clamp.hpp"

anm_result_t anm_bench_custom(int argc, char **argv, param_t *p_param, patch_clamp* p_cell, const char* scenario_name);

#endif
