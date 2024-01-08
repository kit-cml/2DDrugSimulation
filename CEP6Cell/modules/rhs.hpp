#ifndef RHS_HPP
#define RHS_HPP

#include "../commons/helper.hpp"
#include <nvector/nvector_serial.h>

int rhs_fn( realtype t, N_Vector y, N_Vector ydot, void *user_data );

#endif
