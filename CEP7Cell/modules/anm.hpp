#ifndef ANM_HPP
#define ANM_HPP

#include "../commons/helper.hpp"

#if  defined TN2004ENDO || defined TN2004M || defined TN2004EPI || defined TN2006ENDO || defined TN2006M || defined TN2006EPI
int ANM(int argc, char **argv, param_t *p_param);
double get_ANM( double *apd_arr, int length );
double average( double *apd_arr, int start_idx, int end_idx );
#endif

#endif
