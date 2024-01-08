#include "Fibrosis.hpp"
#include "../enums/enum_courtemanche_ramirez_nattel_1998.hpp"

Fibrosis::Fibrosis()
{}


void Fibrosis::mutate( double *ALGEBRAIC, double *STATES, double *CONSTANTS )
{
    ALGEBRAIC[i_K1]   *= 0.5;
    ALGEBRAIC[i_CaL] *= 0.25;
    ALGEBRAIC[i_Na]   *= 0.4;
}
