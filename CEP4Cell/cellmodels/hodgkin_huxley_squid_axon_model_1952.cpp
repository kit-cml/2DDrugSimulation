/*
   There are a total of 10 entries in the algebraic variable array.
   There are a total of 4 entries in each of the rate and state variable arrays.
   There are a total of 8+3 entries in the constant variable array.
 */

#include "hodgkin_huxley_squid_axon_model_1952.hpp"
#include <cmath>

/*
 * VOI is time in component environment (millisecond).
 * STATES[0] is V in component membrane (millivolt).
 * CONSTANTS[0] is E_R in component membrane (millivolt).
 * CONSTANTS[1] is Cm in component membrane (microF_per_cm2).
 * ALGEBRAIC[4] is i_Na in component sodium_channel (microA_per_cm2).
 * ALGEBRAIC[8] is i_K in component potassium_channel (microA_per_cm2).
 * ALGEBRAIC[9] is i_L in component leakage_current (microA_per_cm2).
 * ALGEBRAIC[0] is i_Stim in component membrane (microA_per_cm2).
 * CONSTANTS[2] is g_Na in component sodium_channel (milliS_per_cm2).
 * CONSTANTS[5] is E_Na in component sodium_channel (millivolt).
 * STATES[1] is m in component sodium_channel_m_gate (dimensionless).
 * STATES[2] is h in component sodium_channel_h_gate (dimensionless).
 * ALGEBRAIC[1] is alpha_m in component sodium_channel_m_gate (per_millisecond).
 * ALGEBRAIC[5] is beta_m in component sodium_channel_m_gate (per_millisecond).
 * ALGEBRAIC[2] is alpha_h in component sodium_channel_h_gate (per_millisecond).
 * ALGEBRAIC[6] is beta_h in component sodium_channel_h_gate (per_millisecond).
 * CONSTANTS[3] is g_K in component potassium_channel (milliS_per_cm2).
 * CONSTANTS[6] is E_K in component potassium_channel (millivolt).
 * STATES[3] is n in component potassium_channel_n_gate (dimensionless).
 * ALGEBRAIC[3] is alpha_n in component potassium_channel_n_gate (per_millisecond).
 * ALGEBRAIC[7] is beta_n in component potassium_channel_n_gate (per_millisecond).
 * CONSTANTS[4] is g_L in component leakage_current (milliS_per_cm2).
 * CONSTANTS[7] is E_L in component leakage_current (millivolt).
 * RATES[0] is d/dt V in component membrane (millivolt).
 * RATES[1] is d/dt m in component sodium_channel_m_gate (dimensionless).
 * RATES[2] is d/dt h in component sodium_channel_h_gate (dimensionless).
 * RATES[3] is d/dt n in component potassium_channel_n_gate (dimensionless).
 */


hodgkin_huxley_squid_axon_model_1952::hodgkin_huxley_squid_axon_model_1952()
{
algebraic_size = 10;
constants_size = 8;
states_size = 4;
ALGEBRAIC = new double[algebraic_size];
CONSTANTS = new double[constants_size];
RATES = new double[states_size];
STATES = new double[states_size];
}

hodgkin_huxley_squid_axon_model_1952::~hodgkin_huxley_squid_axon_model_1952()
{
delete mutation;
delete []ALGEBRAIC;
delete []CONSTANTS;
delete []RATES;
delete []STATES;
}

void hodgkin_huxley_squid_axon_model_1952::initConsts()
{
STATES[0] = 0;
CONSTANTS[0] = 0;
CONSTANTS[1] = 1;
CONSTANTS[2] = 120;
STATES[1] = 0.05;
STATES[2] = 0.6;
CONSTANTS[3] = 36;
STATES[3] = 0.325;
CONSTANTS[4] = 0.3;
CONSTANTS[5] = CONSTANTS[0] - 115.000;
CONSTANTS[6] = CONSTANTS[0]+12.0000;
CONSTANTS[7] = CONSTANTS[0] - 10.6130;

CONSTANTS[8] = 0.;
CONSTANTS[9] = 1000.;
CONSTANTS[10] = 0.;
}

void hodgkin_huxley_squid_axon_model_1952::computeRates( double TIME, double *CONSTANTS, double *RATES, double *STATES, double *ALGEBRAIC )
{
ALGEBRAIC[1] = ( 0.100000*(STATES[0]+25.0000))/(exp((STATES[0]+25.0000)/10.0000) - 1.00000);
ALGEBRAIC[5] =  4.00000*exp(STATES[0]/18.0000);
ALGEBRAIC[2] =  0.0700000*exp(STATES[0]/20.0000);
ALGEBRAIC[6] = 1.00000/(exp((STATES[0]+30.0000)/10.0000)+1.00000);
ALGEBRAIC[3] = ( 0.0100000*(STATES[0]+10.0000))/(exp((STATES[0]+10.0000)/10.0000) - 1.00000);
ALGEBRAIC[7] =  0.125000*exp(STATES[0]/80.0000);
ALGEBRAIC[4] =  CONSTANTS[2]*pow(STATES[1], 3.00000)*STATES[2]*(STATES[0] - CONSTANTS[5]);
ALGEBRAIC[8] =  CONSTANTS[3]*pow(STATES[3], 4.00000)*(STATES[0] - CONSTANTS[6]);
ALGEBRAIC[9] =  CONSTANTS[4]*(STATES[0] - CONSTANTS[7]);
ALGEBRAIC[0] = (TIME>=10.0000&&TIME<=10.5000 ? - 20.0000 : 0.000000);

if( mutation != 0){
  mutation->mutate( ALGEBRAIC, STATES, CONSTANTS );
}

RATES[1] =  ALGEBRAIC[1]*(1.00000 - STATES[1]) -  ALGEBRAIC[5]*STATES[1];
RATES[2] =  ALGEBRAIC[2]*(1.00000 - STATES[2]) -  ALGEBRAIC[6]*STATES[2];
RATES[3] =  ALGEBRAIC[3]*(1.00000 - STATES[3]) -  ALGEBRAIC[7]*STATES[3];
RATES[0] = - (- ALGEBRAIC[0]+ALGEBRAIC[4]+ALGEBRAIC[8]+ALGEBRAIC[9])/CONSTANTS[1];
}

void hodgkin_huxley_squid_axon_model_1952::solveAnalytical(double dt)
{
}
