#include "S140G.hpp"
#include "../enums/enum_courtemanche_ramirez_nattel_1998.hpp"

S140G::S140G( double phi )
{
  this->phi = phi;
}


void S140G::mutate(double* ALGEBRAIC, double* STATES, double* CONSTANTS)
{
  ALGEBRAIC[i_Ks] += CONSTANTS[Cm] * phi * CONSTANTS[g_Ks] * ( STATES[0] - (-75.3) );
}
