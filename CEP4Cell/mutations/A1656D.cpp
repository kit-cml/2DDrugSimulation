#include "A1656D.hpp"
#include "../enums/enum_tentusscher_noble_noble_panfilov_2004.hpp"

#include <cmath>
#include <cstring>

A1656D::A1656D( const char* str )
{
  if(strncmp(str, "WT", sizeof(str)) == 0) type = 1;
  else if(strncmp(str, "A1656D", sizeof(str)) == 0) type = 2;
  else if(strncmp(str, "MEX", sizeof(str)) == 0) type = 3;
  else if(strncmp(str, "FLE", sizeof(str)) == 0) type = 4;
  else if(strncmp(str, "RAN", sizeof(str)) == 0) type = 5;
  else type = 0;
}


void A1656D::mutate(double* ALGEBRAIC, double* STATES, double* CONSTANTS)
{
  double wt[21]        = {1125.0,
                          -36.445, -4.519, -52.000, 1.206, -46.861, -126.579, 0.226,
                          -78.929, 4.409, -50.0, 0.980, 117.620, 96.289, 0.314,
                          -3.0, 57.314, 104.658, 59.032, 0.0,
                          0.990};
  const double mt[21]  = {213.66,
                          -44.694, -5.369, -68.764, 5.078, 196.176, 85.086, 0.,
                          -69.476, 4.898, -36.690, 13.568, -181.416, 2.772, 0.,
                          -3.423, 449.708, 128.855, 165.944, 17.551,
                          0.863};
  const double mex[21] = {231.18,
                          -47.173, -4.998, -62.604, 6.202, 256.571, 126.690, 0.106,
                          -71.958, 5.486, -17.961, 12.620, 13.7, 82.331, 0.0,
                          -10.432, 89.477, 48.360, 160.984, 6.627,
                          0.902};
  const double fle[21] = {259.03,
                          -54.609, -3.615, -61.536, 1.347, 122.632, 53.127, 0.394,
                          -73.533, 8.220, 53.838, 11.021, 168.053, 35.413, 0.0,
                          -10.138, 9.198, -15.989, 362.571, 17.084,
                          0.5}; 
  const double ran[21] = {261.05,
                          -50.207, -4.950, -63.290, 3.394, 161.294, 90.247, 0.135,
                          -70.530, 4.886, 14.483, 19.607, 30.882, 84.456, 1.745,
                          6.126, 91.069, -103.071, -84.140, 5.697,
                          0.880}; 
  const double RT = CONSTANTS[R]*CONSTANTS[T];
  double mtemp, htemp, jtemp;

  switch (type){
    case 2: 
      memcpy(wt, mt, sizeof(wt));
      break;
    case 3:
      memcpy(wt, mex, sizeof(wt));
      break;
    case 4:
      memcpy(wt, fle, sizeof(wt));
      break;
    case 5:
      memcpy(wt, ran, sizeof(wt));
      break;
    default:
      break;
  }

  ALGEBRAIC[tau_m] = wt[4] / (exp(wt[5] * (STATES[V] - wt[3]) / RT) + exp(-wt[6] * (STATES[V] - wt[3]) / RT)) + wt[7];
  ALGEBRAIC[m_inf] = pow( 1.0 / (1.0 + exp((STATES[V]-wt[1]) / wt[2])), 1.0/3.0 );
  ALGEBRAIC[tau_h] = wt[11] / (exp(wt[12] * (STATES[V] - wt[10]) / RT) + exp(-wt[13] * (STATES[V] - wt[10]) / RT)) + wt[14];
  ALGEBRAIC[h_inf] = 1.0 / (1.0 + exp((STATES[V] - wt[8]) / wt[9]));
  ALGEBRAIC[tau_j] = wt[16] / (exp(wt[17] * (STATES[V] - wt[15]) / RT) + exp(-wt[18] * (STATES[V] - wt[15]) / RT)) + wt[19];
  ALGEBRAIC[j_inf] = ALGEBRAIC[h_inf];

//  ALGEBRAIC[E_Na] = 50.0;

  mtemp = ALGEBRAIC[m_inf]-( ( ALGEBRAIC[m_inf]-STATES[m] )*exp( -0.01/ALGEBRAIC[tau_m] ) ); 
  htemp = ALGEBRAIC[h_inf]-( ( ALGEBRAIC[h_inf]-STATES[h] )*exp( -0.01/ALGEBRAIC[tau_h] ) ); 
  jtemp = ALGEBRAIC[j_inf]-( ( ALGEBRAIC[j_inf]-STATES[j] )*exp( -0.01/ALGEBRAIC[tau_j] ) ); 

  ALGEBRAIC[i_Na] = CONSTANTS[Cm]*wt[0]*
                    pow(mtemp, 3.0)*
                    ((wt[20]*htemp) + ((1.0-wt[20])*jtemp))*
                    (STATES[0] - ALGEBRAIC[E_Na]);
}
