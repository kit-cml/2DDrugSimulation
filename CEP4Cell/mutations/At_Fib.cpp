#include "At_Fib.hpp"
#include "../enums/enum_courtemanche_ramirez_nattel_1998.hpp"

At_Fib::At_Fib( )
{}


void At_Fib::mutate( double *ALGEBRAIC, double *STATES, double *CONSTANTS)
{
    ALGEBRAIC[i_Kur]   *= 0.7;
    ALGEBRAIC[i_to]    *= 0.5;
    ALGEBRAIC[i_CaL]  *= 0.5;
}
