#include "G229D.hpp"
#if defined(ORUDY_CIPA2017)
  #include "../enums/enum_ohara_rudy_cipa_v1_2017.hpp"
#else
  #include "../enums/enum_Ohara_Rudy_2011.hpp"
#endif


#include <cmath>
#include <cstring>

G229D::G229D( const char* str )
{
  if(strncmp(str, "G229D", sizeof(str)) == 0) is_WT = false;
  else is_WT = true;
}



void G229D::mutate(double* ALGEBRAIC, double* STATES, double* CONSTANTS)
{
  double wt[21]        = {1, 28.8, 15.45, 
                          65.5, 227.2};
  const double mt[21]  = {0.85, 82.8, 41.72, 
                          119.5, 281.2};
  double xs1temp,xs2temp;


  if(!is_WT) memcpy(wt, mt, sizeof(wt));

  ALGEBRAIC[xs1ss] = wt[0]/( 1+exp( -( STATES[V]+wt[1] )/wt[2] ) );
  ALGEBRAIC[xs2ss] = ALGEBRAIC[xs1ss];
  ALGEBRAIC[txs1]  = 326.9+0.4/( 2.326E-4*exp( ( STATES[V]+wt[3] )/17.8 ) +
                     1.292E-3*exp( -( STATES[V]+wt[4] )/230 ) );
  ALGEBRAIC[txs2]  = 5/( 0.01*exp( ( STATES[V]-50 )/100 )+
                     0.0193*exp( -( STATES[V]+66.54 )/155 ) );
  xs1temp = ( ( STATES[xs1]-ALGEBRAIC[xs1ss] )*exp( -0.01/ALGEBRAIC[txs1] ) ) + ALGEBRAIC[xs1ss];
  xs2temp = ( ( STATES[xs2]-ALGEBRAIC[xs2ss] )*exp( -0.01/ALGEBRAIC[txs2] ) ) + ALGEBRAIC[xs2ss];
  ALGEBRAIC[i_Ks]  = CONSTANTS[GKs] * xs1temp * xs2temp * ( STATES[V]-ALGEBRAIC[EKs] ) *
                     ( 1+0.6/( 1+pow( ( 3.8E-5/STATES[Ca_i] ), 1.4) ) );
}
