#ifndef APDR_HPP
#define APDR_HPP

#include "../commons/helper.hpp"

#if  defined TN2004ENDO || defined TN2004M || defined TN2004EPI || defined TN2006ENDO || defined TN2006M || defined TN2006EPI
int APDR(int argc, char **argv, param_t *p_param);
#endif

#endif
