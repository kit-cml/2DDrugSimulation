#include "hERG.hpp"
#include "../enums/enum_courtemanche_ramirez_nattel_1998.hpp"

#include <cmath>
#include <cstring>

hERG::hERG(const char* str)
{
  if(strncmp(str, "L532P", sizeof(str)) == 0) type = 1;
  else if(strncmp(str, "N588K", sizeof(str)) == 0) type = 2;
  else type = 0;
}

void hERG::mutate(double* ALGEBRAIC, double* STATES, double* CONSTANTS)
{
  double CTL[12]         = {0., 0.029412, -15.0, 22.4, 14.1, 6.5, 0.0003, 14.1, -5, -3.3328, 5.1237, 1};
  const double L532P[12] = {0., 0.091720, -15.54, 24.37, -9.88, 22.31, 0.00025, -196.86, -131.36, -40.00, 3.79e-6, 1};
  const double N588K[12] = {0., 0.029412, -38.65, 19.46, 16.49, 6.76, 0.0003, 14.1, -5, -3.3328, 5.1237, 2};

  switch (type){
    case 1:
      memcpy(CTL, L532P, sizeof(CTL));
      break;
    case 2:
      memcpy(CTL, N588K, sizeof(CTL));
      break;
    default:
      break;
  }

  ALGEBRAIC[alpha_xr]    = ( CTL[6]*(STATES[0]+CTL[7]) ) / ( 1 - exp((STATES[0]+CTL[7])/CTL[8]) );
  ALGEBRAIC[beta_xr]     = ( 7.3898e-5*(STATES[0]+CTL[9]) ) / ( exp((STATES[0] + CTL[9])/CTL[10]) - 1 );
  ALGEBRAIC[tau_xr]      = pow( (ALGEBRAIC[alpha_xr]+ALGEBRAIC[beta_xr])*CTL[11],-1 );
  ALGEBRAIC[xr_infinity] = pow( 1+exp((STATES[0]+CTL[4])/-CTL[5]), -1 );
  ALGEBRAIC[i_Kr]        = ( CONSTANTS[Cm]*CTL[1]*STATES[xr]*(STATES[0] - ALGEBRAIC[E_K]) ) / 
                          ( 1+exp((STATES[0]+CTL[2])/CTL[3]) );
}
