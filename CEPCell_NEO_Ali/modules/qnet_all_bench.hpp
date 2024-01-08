#ifndef QNET_ALL_BENCH_HPP
#define QNET_ALL_BENCH_HPP

#include "commons.hpp"
#include "../cellmodels/patch_clamp.hpp"

cipa_t qnet_all_bench(int argc, char **argv, param_t *p_param, patch_clamp* p_cell, const char* scenario_name);

#endif
