#include "rhs.hpp"

int rhs_fn( realtype t, N_Vector y, N_Vector ydot, void *user_data )
{
  CellML *data = (CellML*)user_data;
  data->computeRates( t,
                      data->CONSTANTS,
                      N_VGetArrayPointer_Serial(ydot),
                      N_VGetArrayPointer_Serial(y),
                      data->ALGEBRAIC );
  return 0;
}
