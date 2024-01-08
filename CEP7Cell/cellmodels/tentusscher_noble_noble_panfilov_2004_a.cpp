/*
   There are a total of 67 entries in the algebraic variable array.
   There are a total of 17 entries in each of the rate and state variable arrays.
   There are a total of 46 entries in the constant variable array.
 */

#include "tentusscher_noble_noble_panfilov_2004_a.hpp"
#include <cmath>
#include <cstdlib>
#include <cstdio>

/*
 * TIME is time in component environment (millisecond).
 * STATES[V] is V in component membrane (millivolt).
 * CONSTANTS[R] is R in component membrane (joule_per_mole_kelvin).
 * CONSTANTS[T] is T in component membrane (kelvin).
 * CONSTANTS[F] is F in component membrane (coulomb_per_millimole).
 * CONSTANTS[Cm] is Cm in component membrane (microF).
 * CONSTANTS[V_c] is V_c in component membrane (micrometre3).
 * ALGEBRAIC[i_K1] is i_K1 in component inward_rectifier_potassium_current (picoA_per_picoF).
 * ALGEBRAIC[i_to] is i_to in component transient_outward_current (picoA_per_picoF).
 * ALGEBRAIC[i_Kr] is i_Kr in component rapid_time_dependent_potassium_current (picoA_per_picoF).
 * ALGEBRAIC[i_Ks] is i_Ks in component slow_time_dependent_potassium_current (picoA_per_picoF).
 * ALGEBRAIC[i_CaL] is i_CaL in component L_type_Ca_current (picoA_per_picoF).
 * ALGEBRAIC[i_NaK] is i_NaK in component sodium_potassium_pump_current (picoA_per_picoF).
 * ALGEBRAIC[i_Na] is i_Na in component fast_sodium_current (picoA_per_picoF).
 * ALGEBRAIC[i_b_Na] is i_b_Na in component sodium_background_current (picoA_per_picoF).
 * ALGEBRAIC[i_NaCa] is i_NaCa in component sodium_calcium_exchanger_current (picoA_per_picoF).
 * ALGEBRAIC[i_b_Ca] is i_b_Ca in component calcium_background_current (picoA_per_picoF).
 * ALGEBRAIC[i_p_K] is i_p_K in component potassium_pump_current (picoA_per_picoF).
 * ALGEBRAIC[i_p_Ca] is i_p_Ca in component calcium_pump_current (picoA_per_picoF).
 * ALGEBRAIC[i_Stim] is i_Stim in component membrane (picoA_per_picoF).
 * CONSTANTS[stim_start] is stim_start in component membrane (millisecond).
 * CONSTANTS[stim_period] is stim_period in component membrane (millisecond).
 * CONSTANTS[stim_duration] is stim_duration in component membrane (millisecond).
 * CONSTANTS[stim_amplitude] is stim_amplitude in component membrane (picoA_per_picoF).
 * ALGEBRAIC[E_Na] is E_Na in component reversal_potentials (millivolt).
 * ALGEBRAIC[E_K] is E_K in component reversal_potentials (millivolt).
 * ALGEBRAIC[E_Ks] is E_Ks in component reversal_potentials (millivolt).
 * ALGEBRAIC[E_Ca] is E_Ca in component reversal_potentials (millivolt).
 * CONSTANTS[P_kna] is P_kna in component reversal_potentials (dimensionless).
 * CONSTANTS[K_o] is K_o in component potassium_dynamics (millimolar).
 * CONSTANTS[Na_o] is Na_o in component sodium_dynamics (millimolar).
 * STATES[K_i] is K_i in component potassium_dynamics (millimolar).
 * STATES[Na_i] is Na_i in component sodium_dynamics (millimolar).
 * CONSTANTS[Ca_o] is Ca_o in component calcium_dynamics (millimolar).
 * STATES[Ca_i] is Ca_i in component calcium_dynamics (millimolar).
 * CONSTANTS[g_K1] is g_K1 in component inward_rectifier_potassium_current (nanoS_per_picoF).
 * ALGEBRAIC[xK1_inf] is xK1_inf in component inward_rectifier_potassium_current (dimensionless).
 * ALGEBRAIC[alpha_K1] is alpha_K1 in component inward_rectifier_potassium_current (dimensionless).
 * ALGEBRAIC[beta_K1] is beta_K1 in component inward_rectifier_potassium_current (dimensionless).
 * CONSTANTS[g_Kr] is g_Kr in component rapid_time_dependent_potassium_current (nanoS_per_picoF).
 * STATES[Xr1] is Xr1 in component rapid_time_dependent_potassium_current_Xr1_gate (dimensionless).
 * STATES[Xr2] is Xr2 in component rapid_time_dependent_potassium_current_Xr2_gate (dimensionless).
 * ALGEBRAIC[xr1_inf] is xr1_inf in component rapid_time_dependent_potassium_current_Xr1_gate (dimensionless).
 * ALGEBRAIC[alpha_xr1] is alpha_xr1 in component rapid_time_dependent_potassium_current_Xr1_gate (dimensionless).
 * ALGEBRAIC[beta_xr1] is beta_xr1 in component rapid_time_dependent_potassium_current_Xr1_gate (dimensionless).
 * ALGEBRAIC[tau_xr1] is tau_xr1 in component rapid_time_dependent_potassium_current_Xr1_gate (millisecond).
 * ALGEBRAIC[xr2_inf] is xr2_inf in component rapid_time_dependent_potassium_current_Xr2_gate (dimensionless).
 * ALGEBRAIC[alpha_xr2] is alpha_xr2 in component rapid_time_dependent_potassium_current_Xr2_gate (dimensionless).
 * ALGEBRAIC[beta_xr2] is beta_xr2 in component rapid_time_dependent_potassium_current_Xr2_gate (dimensionless).
 * ALGEBRAIC[tau_xr2] is tau_xr2 in component rapid_time_dependent_potassium_current_Xr2_gate (millisecond).
 * CONSTANTS[g_Ks] is g_Ks in component slow_time_dependent_potassium_current (nanoS_per_picoF).
 * STATES[Xs] is Xs in component slow_time_dependent_potassium_current_Xs_gate (dimensionless).
 * ALGEBRAIC[xs_inf] is xs_inf in component slow_time_dependent_potassium_current_Xs_gate (dimensionless).
 * ALGEBRAIC[alpha_xs] is alpha_xs in component slow_time_dependent_potassium_current_Xs_gate (dimensionless).
 * ALGEBRAIC[beta_xs] is beta_xs in component slow_time_dependent_potassium_current_Xs_gate (dimensionless).
 * ALGEBRAIC[tau_xs] is tau_xs in component slow_time_dependent_potassium_current_Xs_gate (millisecond).
 * CONSTANTS[g_Na] is g_Na in component fast_sodium_current (nanoS_per_picoF).
 * STATES[m] is m in component fast_sodium_current_m_gate (dimensionless).
 * STATES[h] is h in component fast_sodium_current_h_gate (dimensionless).
 * STATES[j] is j in component fast_sodium_current_j_gate (dimensionless).
 * ALGEBRAIC[m_inf] is m_inf in component fast_sodium_current_m_gate (dimensionless).
 * ALGEBRAIC[alpha_m] is alpha_m in component fast_sodium_current_m_gate (dimensionless).
 * ALGEBRAIC[beta_m] is beta_m in component fast_sodium_current_m_gate (dimensionless).
 * ALGEBRAIC[tau_m] is tau_m in component fast_sodium_current_m_gate (millisecond).
 * ALGEBRAIC[h_inf] is h_inf in component fast_sodium_current_h_gate (dimensionless).
 * ALGEBRAIC[alpha_h] is alpha_h in component fast_sodium_current_h_gate (per_millisecond).
 * ALGEBRAIC[beta_h] is beta_h in component fast_sodium_current_h_gate (per_millisecond).
 * ALGEBRAIC[tau_h] is tau_h in component fast_sodium_current_h_gate (millisecond).
 * ALGEBRAIC[j_inf] is j_inf in component fast_sodium_current_j_gate (dimensionless).
 * ALGEBRAIC[alpha_j] is alpha_j in component fast_sodium_current_j_gate (per_millisecond).
 * ALGEBRAIC[beta_j] is beta_j in component fast_sodium_current_j_gate (per_millisecond).
 * ALGEBRAIC[tau_j] is tau_j in component fast_sodium_current_j_gate (millisecond).
 * CONSTANTS[g_bna] is g_bna in component sodium_background_current (nanoS_per_picoF).
 * CONSTANTS[g_CaL] is g_CaL in component L_type_Ca_current (litre_per_farad_second).
 * STATES[d] is d in component L_type_Ca_current_d_gate (dimensionless).
 * STATES[f] is f in component L_type_Ca_current_f_gate (dimensionless).
 * STATES[fCa] is fCa in component L_type_Ca_current_fCa_gate (dimensionless).
 * ALGEBRAIC[d_inf] is d_inf in component L_type_Ca_current_d_gate (dimensionless).
 * ALGEBRAIC[alpha_d] is alpha_d in component L_type_Ca_current_d_gate (dimensionless).
 * ALGEBRAIC[beta_d] is beta_d in component L_type_Ca_current_d_gate (dimensionless).
 * ALGEBRAIC[gamma_d] is gamma_d in component L_type_Ca_current_d_gate (millisecond).
 * ALGEBRAIC[tau_d] is tau_d in component L_type_Ca_current_d_gate (millisecond).
 * ALGEBRAIC[f_inf] is f_inf in component L_type_Ca_current_f_gate (dimensionless).
 * ALGEBRAIC[tau_f] is tau_f in component L_type_Ca_current_f_gate (millisecond).
 * ALGEBRAIC[alpha_fCa] is alpha_fCa in component L_type_Ca_current_fCa_gate (dimensionless).
 * ALGEBRAIC[beta_fCa] is beta_fCa in component L_type_Ca_current_fCa_gate (dimensionless).
 * ALGEBRAIC[gama_fCa] is gama_fCa in component L_type_Ca_current_fCa_gate (dimensionless).
 * ALGEBRAIC[fCa_inf] is fCa_inf in component L_type_Ca_current_fCa_gate (dimensionless).
 * CONSTANTS[tau_fCa] is tau_fCa in component L_type_Ca_current_fCa_gate (millisecond).
 * ALGEBRAIC[d_fCa] is d_fCa in component L_type_Ca_current_fCa_gate (per_millisecond).
 * CONSTANTS[g_bca] is g_bca in component calcium_background_current (nanoS_per_picoF).
 * CONSTANTS[g_to] is g_to in component transient_outward_current (nanoS_per_picoF).
 * STATES[s] is s in component transient_outward_current_s_gate (dimensionless).
 * STATES[r] is r in component transient_outward_current_r_gate (dimensionless).
 * ALGEBRAIC[s_inf] is s_inf in component transient_outward_current_s_gate (dimensionless).
 * ALGEBRAIC[tau_s] is tau_s in component transient_outward_current_s_gate (millisecond).
 * ALGEBRAIC[r_inf] is r_inf in component transient_outward_current_r_gate (dimensionless).
 * ALGEBRAIC[tau_r] is tau_r in component transient_outward_current_r_gate (millisecond).
 * CONSTANTS[P_NaK] is P_NaK in component sodium_potassium_pump_current (picoA_per_picoF).
 * CONSTANTS[K_mk] is K_mk in component sodium_potassium_pump_current (millimolar).
 * CONSTANTS[K_mNa] is K_mNa in component sodium_potassium_pump_current (millimolar).
 * CONSTANTS[K_NaCa] is K_NaCa in component sodium_calcium_exchanger_current (picoA_per_picoF).
 * CONSTANTS[K_sat] is K_sat in component sodium_calcium_exchanger_current (dimensionless).
 * CONSTANTS[alpha] is alpha in component sodium_calcium_exchanger_current (dimensionless).
 * CONSTANTS[gama] is gama in component sodium_calcium_exchanger_current (dimensionless).
 * CONSTANTS[Km_Ca] is Km_Ca in component sodium_calcium_exchanger_current (millimolar).
 * CONSTANTS[Km_Nai] is Km_Nai in component sodium_calcium_exchanger_current (millimolar).
 * CONSTANTS[g_pCa] is g_pCa in component calcium_pump_current (picoA_per_picoF).
 * CONSTANTS[K_pCa] is K_pCa in component calcium_pump_current (millimolar).
 * CONSTANTS[g_pK] is g_pK in component potassium_pump_current (nanoS_per_picoF).
 * STATES[Ca_SR] is Ca_SR in component calcium_dynamics (millimolar).
 * ALGEBRAIC[i_rel] is i_rel in component calcium_dynamics (millimolar_per_millisecond).
 * ALGEBRAIC[i_up] is i_up in component calcium_dynamics (millimolar_per_millisecond).
 * ALGEBRAIC[i_leak] is i_leak in component calcium_dynamics (millimolar_per_millisecond).
 * STATES[g] is g in component calcium_dynamics (dimensionless).
 * CONSTANTS[tau_g] is tau_g in component calcium_dynamics (millisecond).
 * ALGEBRAIC[g_inf] is g_inf in component calcium_dynamics (dimensionless).
 * CONSTANTS[a_rel] is a_rel in component calcium_dynamics (millimolar_per_millisecond).
 * CONSTANTS[b_rel] is b_rel in component calcium_dynamics (millimolar).
 * CONSTANTS[c_rel] is c_rel in component calcium_dynamics (millimolar_per_millisecond).
 * CONSTANTS[K_up] is K_up in component calcium_dynamics (millimolar).
 * CONSTANTS[V_leak] is V_leak in component calcium_dynamics (per_millisecond).
 * CONSTANTS[Vmax_up] is Vmax_up in component calcium_dynamics (millimolar_per_millisecond).
 * ALGEBRAIC[Ca_i_bufc] is Ca_i_bufc in component calcium_dynamics (dimensionless).
 * ALGEBRAIC[Ca_sr_bufsr] is Ca_sr_bufsr in component calcium_dynamics (dimensionless).
 * CONSTANTS[Buf_c] is Buf_c in component calcium_dynamics (millimolar).
 * CONSTANTS[K_buf_c] is K_buf_c in component calcium_dynamics (millimolar).
 * CONSTANTS[Buf_sr] is Buf_sr in component calcium_dynamics (millimolar).
 * CONSTANTS[K_buf_sr] is K_buf_sr in component calcium_dynamics (millimolar).
 * CONSTANTS[V_sr] is V_sr in component calcium_dynamics (micrometre3).
 * ALGEBRAIC[d_g] is d_g in component calcium_dynamics (per_millisecond).
 * RATES[V] is d/dt V in component membrane (millivolt).
 * RATES[Xr1] is d/dt Xr1 in component rapid_time_dependent_potassium_current_Xr1_gate (dimensionless).
 * RATES[Xr2] is d/dt Xr2 in component rapid_time_dependent_potassium_current_Xr2_gate (dimensionless).
 * RATES[Xs] is d/dt Xs in component slow_time_dependent_potassium_current_Xs_gate (dimensionless).
 * RATES[m] is d/dt m in component fast_sodium_current_m_gate (dimensionless).
 * RATES[h] is d/dt h in component fast_sodium_current_h_gate (dimensionless).
 * RATES[j] is d/dt j in component fast_sodium_current_j_gate (dimensionless).
 * RATES[d] is d/dt d in component L_type_Ca_current_d_gate (dimensionless).
 * RATES[f] is d/dt f in component L_type_Ca_current_f_gate (dimensionless).
 * RATES[fCa] is d/dt fCa in component L_type_Ca_current_fCa_gate (dimensionless).
 * RATES[s] is d/dt s in component transient_outward_current_s_gate (dimensionless).
 * RATES[r] is d/dt r in component transient_outward_current_r_gate (dimensionless).
 * RATES[g] is d/dt g in component calcium_dynamics (dimensionless).
 * RATES[Ca_i] is d/dt Ca_i in component calcium_dynamics (millimolar).
 * RATES[Ca_SR] is d/dt Ca_SR in component calcium_dynamics (millimolar).
 * RATES[Na_i] is d/dt Na_i in component sodium_dynamics (millimolar).
 * RATES[K_i] is d/dt K_i in component potassium_dynamics (millimolar).
 */


tentusscher_noble_noble_panfilov_2004_a::tentusscher_noble_noble_panfilov_2004_a()
{
algebraic_size = 67;
constants_size = 46;
states_size = 17;
ALGEBRAIC = new double[algebraic_size];
CONSTANTS = new double[constants_size];
RATES = new double[states_size];
STATES = new double[states_size];
}

tentusscher_noble_noble_panfilov_2004_a::~tentusscher_noble_noble_panfilov_2004_a()
{
if( mutation != 0 )
{
  delete mutation;
}
delete []ALGEBRAIC;
delete []CONSTANTS;
delete []RATES;
delete []STATES;
}

void tentusscher_noble_noble_panfilov_2004_a::initConsts()
{
STATES[V] = -86.2;
CONSTANTS[R] = 8314.472;
CONSTANTS[T] = 310;
CONSTANTS[F] = 96485.3415;
CONSTANTS[Cm] = 0.185;
CONSTANTS[V_c] = 0.016404;
CONSTANTS[stim_start] = 10;
CONSTANTS[stim_period] = 1000;
CONSTANTS[stim_duration] = 1;
CONSTANTS[stim_amplitude] = 52;
CONSTANTS[P_kna] = 0.03;
CONSTANTS[K_o] = 5.4;
CONSTANTS[Na_o] = 140;
STATES[K_i] = 138.3;
STATES[Na_i] = 11.6;
CONSTANTS[Ca_o] = 2;
STATES[Ca_i] = 0.0002;
CONSTANTS[g_K1] = 5.405;
CONSTANTS[g_Kr] = 0.096;
STATES[Xr1] = 0;
STATES[Xr2] = 1;
CONSTANTS[g_Ks] = 0.062;
STATES[Xs] = 0;
CONSTANTS[g_Na] = 14.838;
STATES[m] = 0;
STATES[h] = 0.75;
STATES[j] = 0.75;
CONSTANTS[g_bna] = 0.00029;
CONSTANTS[g_CaL] = 0.000175;
STATES[d] = 0;
STATES[f] = 1;
STATES[fCa] = 1;
CONSTANTS[g_bca] = 0.000592;
CONSTANTS[g_to] = 0.294;
STATES[s] = 1;
STATES[r] = 0;
CONSTANTS[P_NaK] = 1.362;
CONSTANTS[K_mk] = 1;
CONSTANTS[K_mNa] = 40;
CONSTANTS[K_NaCa] = 1000;
CONSTANTS[K_sat] = 0.1;
CONSTANTS[alpha] = 2.5;
CONSTANTS[gama] = 0.35;
CONSTANTS[Km_Ca] = 1.38;
CONSTANTS[Km_Nai] = 87.5;
CONSTANTS[g_pCa] = 0.825;
CONSTANTS[K_pCa] = 0.0005;
CONSTANTS[g_pK] = 0.0146;
STATES[Ca_SR] = 0.2;
STATES[g] = 1;
CONSTANTS[tau_g] = 2;
CONSTANTS[a_rel] = 0.016464;
CONSTANTS[b_rel] = 0.25;
CONSTANTS[c_rel] = 0.008232;
CONSTANTS[K_up] = 0.00025;
CONSTANTS[V_leak] = 8e-5;
CONSTANTS[Vmax_up] = 0.000425;
CONSTANTS[Buf_c] = 0.15;
CONSTANTS[K_buf_c] = 0.001;
CONSTANTS[Buf_sr] = 10;
CONSTANTS[K_buf_sr] = 0.3;
CONSTANTS[V_sr] = 0.001094;
CONSTANTS[tau_fCa] = 2.00000;
}

void tentusscher_noble_noble_panfilov_2004_a::computeRates( double TIME, double *CONSTANTS, double *RATES, double *STATES, double *ALGEBRAIC )
{
ALGEBRAIC[f_inf] = 1.00000/(1.00000+exp((STATES[V]+20.0000)/7.00000));
ALGEBRAIC[tau_f] =  1125.00*exp(- pow(STATES[V]+27.0000, 2.00000)/240.000)+80.0000+165.000/(1.00000+exp((25.0000 - STATES[V])/10.0000));
ALGEBRAIC[s_inf] = 1.00000/(1.00000+exp((STATES[V]+20.0000)/5.00000));
ALGEBRAIC[tau_s] =  85.0000*exp(- pow(STATES[V]+45.0000, 2.00000)/320.000)+5.00000/(1.00000+exp((STATES[V] - 20.0000)/5.00000))+3.00000;
ALGEBRAIC[r_inf] = 1.00000/(1.00000+exp((20.0000 - STATES[V])/6.00000));
ALGEBRAIC[tau_r] =  9.50000*exp(- pow(STATES[V]+40.0000, 2.00000)/1800.00)+0.800000;
ALGEBRAIC[g_inf] = (STATES[Ca_i]<0.000350000 ? 1.00000/(1.00000+pow(STATES[Ca_i]/0.000350000, 6.00000)) : 1.00000/(1.00000+pow(STATES[Ca_i]/0.000350000, 16.0000)));
ALGEBRAIC[d_g] = (ALGEBRAIC[g_inf] - STATES[g])/CONSTANTS[tau_g];
ALGEBRAIC[xr1_inf] = 1.00000/(1.00000+exp((- 26.0000 - STATES[V])/7.00000));
ALGEBRAIC[alpha_xr1] = 450.000/(1.00000+exp((- 45.0000 - STATES[V])/10.0000));
ALGEBRAIC[beta_xr1] = 6.00000/(1.00000+exp((STATES[V]+30.0000)/11.5000));
ALGEBRAIC[tau_xr1] =  1.00000*ALGEBRAIC[alpha_xr1]*ALGEBRAIC[beta_xr1];
ALGEBRAIC[xr2_inf] = 1.00000/(1.00000+exp((STATES[V]+88.0000)/24.0000));
ALGEBRAIC[alpha_xr2] = 3.00000/(1.00000+exp((- 60.0000 - STATES[V])/20.0000));
ALGEBRAIC[beta_xr2] = 1.12000/(1.00000+exp((STATES[V] - 60.0000)/20.0000));
ALGEBRAIC[tau_xr2] =  1.00000*ALGEBRAIC[alpha_xr2]*ALGEBRAIC[beta_xr2];
ALGEBRAIC[xs_inf] = 1.00000/(1.00000+exp((- 5.00000 - STATES[V])/14.0000));
ALGEBRAIC[alpha_xs] = 1100.00/ pow((1.00000+exp((- 10.0000 - STATES[V])/6.00000)), 1.0 / 2);
ALGEBRAIC[beta_xs] = 1.00000/(1.00000+exp((STATES[V] - 60.0000)/20.0000));
ALGEBRAIC[tau_xs] =  1.00000*ALGEBRAIC[alpha_xs]*ALGEBRAIC[beta_xs];
ALGEBRAIC[m_inf] = 1.00000/pow(1.00000+exp((- 56.8600 - STATES[V])/9.03000), 2.00000);
ALGEBRAIC[alpha_m] = 1.00000/(1.00000+exp((- 60.0000 - STATES[V])/5.00000));
ALGEBRAIC[beta_m] = 0.100000/(1.00000+exp((STATES[V]+35.0000)/5.00000))+0.100000/(1.00000+exp((STATES[V] - 50.0000)/200.000));
ALGEBRAIC[tau_m] =  1.00000*ALGEBRAIC[alpha_m]*ALGEBRAIC[beta_m];
ALGEBRAIC[h_inf] = 1.00000/pow(1.00000+exp((STATES[V]+71.5500)/7.43000), 2.00000);
ALGEBRAIC[alpha_h] = (STATES[V]<- 40.0000 ?  0.0570000*exp(- (STATES[V]+80.0000)/6.80000) : 0.000000);
ALGEBRAIC[beta_h] = (STATES[V]<- 40.0000 ?  2.70000*exp( 0.0790000*STATES[V])+ 310000.*exp( 0.348500*STATES[V]) : 0.770000/( 0.130000*(1.00000+exp((STATES[V]+10.6600)/- 11.1000))));
ALGEBRAIC[tau_h] = 1.00000/(ALGEBRAIC[alpha_h]+ALGEBRAIC[beta_h]);
ALGEBRAIC[j_inf] = 1.00000/pow(1.00000+exp((STATES[V]+71.5500)/7.43000), 2.00000);
ALGEBRAIC[alpha_j] = (STATES[V]<- 40.0000 ? (( ( - 25428.0*exp( 0.244400*STATES[V]) -  6.94800e-06*exp( - 0.0439100*STATES[V]))*(STATES[V]+37.7800))/1.00000)/(1.00000+exp( 0.311000*(STATES[V]+79.2300))) : 0.000000);
ALGEBRAIC[beta_j] = (STATES[V]<- 40.0000 ? ( 0.0242400*exp( - 0.0105200*STATES[V]))/(1.00000+exp( - 0.137800*(STATES[V]+40.1400))) : ( 0.600000*exp( 0.0570000*STATES[V]))/(1.00000+exp( - 0.100000*(STATES[V]+32.0000))));
ALGEBRAIC[tau_j] = 1.00000/(ALGEBRAIC[alpha_j]+ALGEBRAIC[beta_j]);
ALGEBRAIC[d_inf] = 1.00000/(1.00000+exp((- 5.00000 - STATES[V])/7.50000));
ALGEBRAIC[alpha_d] = 1.40000/(1.00000+exp((- 35.0000 - STATES[V])/13.0000))+0.250000;
ALGEBRAIC[beta_d] = 1.40000/(1.00000+exp((STATES[V]+5.00000)/5.00000));
ALGEBRAIC[gamma_d] = 1.00000/(1.00000+exp((50.0000 - STATES[V])/20.0000));
ALGEBRAIC[tau_d] =  1.00000*ALGEBRAIC[alpha_d]*ALGEBRAIC[beta_d]+ALGEBRAIC[gamma_d];
ALGEBRAIC[alpha_fCa] = 1.00000/(1.00000+pow(STATES[Ca_i]/0.000325000, 8.00000));
ALGEBRAIC[beta_fCa] = 0.100000/(1.00000+exp((STATES[Ca_i] - 0.000500000)/0.000100000));
ALGEBRAIC[gama_fCa] = 0.200000/(1.00000+exp((STATES[Ca_i] - 0.000750000)/0.000800000));
ALGEBRAIC[fCa_inf] = (ALGEBRAIC[alpha_fCa]+ALGEBRAIC[beta_fCa]+ALGEBRAIC[gama_fCa]+0.230000)/1.46000;
ALGEBRAIC[d_fCa] = (ALGEBRAIC[fCa_inf] - STATES[fCa])/CONSTANTS[tau_fCa];
ALGEBRAIC[i_NaK] = (( (( CONSTANTS[P_NaK]*CONSTANTS[K_o])/(CONSTANTS[K_o]+CONSTANTS[K_mk]))*STATES[Na_i])/(STATES[Na_i]+CONSTANTS[K_mNa]))/(1.00000+ 0.124500*exp(( - 0.100000*STATES[V]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[T]))+ 0.0353000*exp(( - STATES[V]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[T])));
ALGEBRAIC[E_Na] =  (( CONSTANTS[R]*CONSTANTS[T])/CONSTANTS[F])*log(CONSTANTS[Na_o]/STATES[Na_i]);
ALGEBRAIC[i_Na] =  CONSTANTS[g_Na]*pow(STATES[m], 3.00000)*STATES[h]*STATES[j]*(STATES[V] - ALGEBRAIC[E_Na]);
ALGEBRAIC[i_b_Na] =  CONSTANTS[g_bna]*(STATES[V] - ALGEBRAIC[E_Na]);
ALGEBRAIC[i_NaCa] = ( CONSTANTS[K_NaCa]*( exp(( CONSTANTS[gama]*STATES[V]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[T]))*pow(STATES[Na_i], 3.00000)*CONSTANTS[Ca_o] -  exp(( (CONSTANTS[gama] - 1.00000)*STATES[V]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[T]))*pow(CONSTANTS[Na_o], 3.00000)*STATES[Ca_i]*CONSTANTS[alpha]))/( (pow(CONSTANTS[Km_Nai], 3.00000)+pow(CONSTANTS[Na_o], 3.00000))*(CONSTANTS[Km_Ca]+CONSTANTS[Ca_o])*(1.00000+ CONSTANTS[K_sat]*exp(( (CONSTANTS[gama] - 1.00000)*STATES[V]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[T]))));
ALGEBRAIC[E_K] =  (( CONSTANTS[R]*CONSTANTS[T])/CONSTANTS[F])*log(CONSTANTS[K_o]/STATES[K_i]);
ALGEBRAIC[alpha_K1] = 0.100000/(1.00000+exp( 0.0600000*((STATES[V] - ALGEBRAIC[E_K]) - 200.000)));
ALGEBRAIC[beta_K1] = ( 3.00000*exp( 0.000200000*((STATES[V] - ALGEBRAIC[E_K])+100.000))+exp( 0.100000*((STATES[V] - ALGEBRAIC[E_K]) - 10.0000)))/(1.00000+exp( - 0.500000*(STATES[V] - ALGEBRAIC[E_K])));
ALGEBRAIC[xK1_inf] = ALGEBRAIC[alpha_K1]/(ALGEBRAIC[alpha_K1]+ALGEBRAIC[beta_K1]);
ALGEBRAIC[i_K1] =  CONSTANTS[g_K1]*ALGEBRAIC[xK1_inf]* pow((CONSTANTS[K_o]/5.40000), 1.0 / 2)*(STATES[V] - ALGEBRAIC[E_K]);
ALGEBRAIC[i_to] =  CONSTANTS[g_to]*STATES[r]*STATES[s]*(STATES[V] - ALGEBRAIC[E_K]);
ALGEBRAIC[i_Kr] =  CONSTANTS[g_Kr]* pow((CONSTANTS[K_o]/5.40000), 1.0 / 2)*STATES[Xr1]*STATES[Xr2]*(STATES[V] - ALGEBRAIC[E_K]);
ALGEBRAIC[E_Ks] =  (( CONSTANTS[R]*CONSTANTS[T])/CONSTANTS[F])*log((CONSTANTS[K_o]+ CONSTANTS[P_kna]*CONSTANTS[Na_o])/(STATES[K_i]+ CONSTANTS[P_kna]*STATES[Na_i]));
ALGEBRAIC[i_Ks] =  CONSTANTS[g_Ks]*pow(STATES[Xs], 2.00000)*(STATES[V] - ALGEBRAIC[E_Ks]);
ALGEBRAIC[i_CaL] = ( (( CONSTANTS[g_CaL]*STATES[d]*STATES[f]*STATES[fCa]*4.00000*STATES[V]*pow(CONSTANTS[F], 2.00000))/( CONSTANTS[R]*CONSTANTS[T]))*( STATES[Ca_i]*exp(( 2.00000*STATES[V]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[T])) -  0.341000*CONSTANTS[Ca_o]))/(exp(( 2.00000*STATES[V]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[T])) - 1.00000);
ALGEBRAIC[E_Ca] =  (( 0.500000*CONSTANTS[R]*CONSTANTS[T])/CONSTANTS[F])*log(CONSTANTS[Ca_o]/STATES[Ca_i]);
ALGEBRAIC[i_b_Ca] =  CONSTANTS[g_bca]*(STATES[V] - ALGEBRAIC[E_Ca]);
ALGEBRAIC[i_p_K] = ( CONSTANTS[g_pK]*(STATES[V] - ALGEBRAIC[E_K]))/(1.00000+exp((25.0000 - STATES[V])/5.98000));
ALGEBRAIC[i_p_Ca] = ( CONSTANTS[g_pCa]*STATES[Ca_i])/(STATES[Ca_i]+CONSTANTS[K_pCa]);
#if defined(SINGLE_CELL)
ALGEBRAIC[i_Stim] = (TIME -  floor(TIME/CONSTANTS[stim_period])*CONSTANTS[stim_period]>=CONSTANTS[stim_start]&&TIME -  floor(TIME/CONSTANTS[stim_period])*CONSTANTS[stim_period]<=CONSTANTS[stim_start]+CONSTANTS[stim_duration] ? - CONSTANTS[stim_amplitude] : 0.000000);
#else
if(isS1) ALGEBRAIC[0] = - CONSTANTS[8];
else ALGEBRAIC[0] = 0.0;
if(isEctopic) ALGEBRAIC[0] = - CONSTANTS[8];
else {
  if(!isS1) ALGEBRAIC[0] = 0.0;
}
#endif
ALGEBRAIC[i_rel] =  (( CONSTANTS[a_rel]*pow(STATES[Ca_SR], 2.00000))/(pow(CONSTANTS[b_rel], 2.00000)+pow(STATES[Ca_SR], 2.00000))+CONSTANTS[c_rel])*STATES[d]*STATES[g];
ALGEBRAIC[i_up] = CONSTANTS[Vmax_up]/(1.00000+pow(CONSTANTS[K_up], 2.00000)/pow(STATES[Ca_i], 2.00000));
ALGEBRAIC[i_leak] =  CONSTANTS[V_leak]*(STATES[Ca_SR] - STATES[Ca_i]);
ALGEBRAIC[Ca_i_bufc] = 1.00000/(1.00000+( CONSTANTS[Buf_c]*CONSTANTS[K_buf_c])/pow(STATES[Ca_i]+CONSTANTS[K_buf_c], 2.00000));
ALGEBRAIC[Ca_sr_bufsr] = 1.00000/(1.00000+( CONSTANTS[Buf_sr]*CONSTANTS[K_buf_sr])/pow(STATES[Ca_SR]+CONSTANTS[K_buf_sr], 2.00000));

if( mutation != 0){
  mutation->mutate( ALGEBRAIC, STATES, CONSTANTS );
}

RATES[f] = (ALGEBRAIC[f_inf] - STATES[f])/ALGEBRAIC[tau_f];
RATES[s] = (ALGEBRAIC[s_inf] - STATES[s])/ALGEBRAIC[tau_s];
RATES[r] = (ALGEBRAIC[r_inf] - STATES[r])/ALGEBRAIC[tau_r];
RATES[g] = (ALGEBRAIC[g_inf]>STATES[g]&&STATES[V]>- 60.0000 ? 0.000000 : ALGEBRAIC[d_g]);
RATES[Xr1] = (ALGEBRAIC[xr1_inf] - STATES[Xr1])/ALGEBRAIC[tau_xr1];
RATES[Xr2] = (ALGEBRAIC[xr2_inf] - STATES[Xr2])/ALGEBRAIC[tau_xr2];
RATES[Xs] = (ALGEBRAIC[xs_inf] - STATES[Xs])/ALGEBRAIC[tau_xs];
RATES[m] = (ALGEBRAIC[m_inf] - STATES[m])/ALGEBRAIC[tau_m];
RATES[h] = (ALGEBRAIC[h_inf] - STATES[h])/ALGEBRAIC[tau_h];
RATES[j] = (ALGEBRAIC[j_inf] - STATES[j])/ALGEBRAIC[tau_j];
RATES[d] = (ALGEBRAIC[d_inf] - STATES[d])/ALGEBRAIC[tau_d];
RATES[fCa] = (ALGEBRAIC[fCa_inf]>STATES[fCa]&&STATES[V]>- 60.0000 ? 0.000000 : ALGEBRAIC[d_fCa]);
RATES[Na_i] = ( - 1.00000*(ALGEBRAIC[i_Na]+ALGEBRAIC[i_b_Na]+ 3.00000*ALGEBRAIC[i_NaK]+ 3.00000*ALGEBRAIC[i_NaCa])*CONSTANTS[Cm])/( 1.00000*CONSTANTS[V_c]*CONSTANTS[F]);
RATES[V] =  (- 1.00000/1.00000)*(ALGEBRAIC[i_K1]+ALGEBRAIC[i_to]+ALGEBRAIC[i_Kr]+ALGEBRAIC[i_Ks]+ALGEBRAIC[i_CaL]+ALGEBRAIC[i_NaK]+ALGEBRAIC[i_Na]+ALGEBRAIC[i_b_Na]+ALGEBRAIC[i_NaCa]+ALGEBRAIC[i_b_Ca]+ALGEBRAIC[i_p_K]+ALGEBRAIC[i_p_Ca]+ALGEBRAIC[i_Stim]);
RATES[K_i] = ( - 1.00000*((ALGEBRAIC[i_K1]+ALGEBRAIC[i_to]+ALGEBRAIC[i_Kr]+ALGEBRAIC[i_Ks]+ALGEBRAIC[i_p_K]+ALGEBRAIC[i_Stim]) -  2.00000*ALGEBRAIC[i_NaK])*CONSTANTS[Cm])/( 1.00000*CONSTANTS[V_c]*CONSTANTS[F]);
RATES[Ca_i] =  ALGEBRAIC[Ca_i_bufc]*(((ALGEBRAIC[i_leak] - ALGEBRAIC[i_up])+ALGEBRAIC[i_rel]) -  (( 1.00000*((ALGEBRAIC[i_CaL]+ALGEBRAIC[i_b_Ca]+ALGEBRAIC[i_p_Ca]) -  2.00000*ALGEBRAIC[i_NaCa]))/( 2.00000*1.00000*CONSTANTS[V_c]*CONSTANTS[F]))*CONSTANTS[Cm]);
RATES[Ca_SR] =  (( ALGEBRAIC[Ca_sr_bufsr]*CONSTANTS[V_c])/CONSTANTS[V_sr])*(ALGEBRAIC[i_up] - (ALGEBRAIC[i_rel]+ALGEBRAIC[i_leak]));
}

void tentusscher_noble_noble_panfilov_2004_a::solveAnalytical(double dt)
{
STATES[11] = ((STATES[11] - ALGEBRAIC[8])*exp(-dt / ALGEBRAIC[21])) + ALGEBRAIC[8];
STATES[13] = ((STATES[13] - ALGEBRAIC[10])*exp(-dt / ALGEBRAIC[23])) + ALGEBRAIC[10];
STATES[14] = ((STATES[14] - ALGEBRAIC[11])*exp(-dt / ALGEBRAIC[24])) + ALGEBRAIC[11];
STATES[4] = ((STATES[4] - ALGEBRAIC[1])*exp(-dt / ALGEBRAIC[36])) + ALGEBRAIC[1];
STATES[5] = ((STATES[5] - ALGEBRAIC[2])*exp(-dt / ALGEBRAIC[37])) + ALGEBRAIC[2];
STATES[6] = ((STATES[6] - ALGEBRAIC[3])*exp(-dt / ALGEBRAIC[38])) + ALGEBRAIC[3];
STATES[7] = ((STATES[7] - ALGEBRAIC[4])*exp(-dt / ALGEBRAIC[39])) + ALGEBRAIC[4];
STATES[8] = ((STATES[8] - ALGEBRAIC[5])*exp(-dt / ALGEBRAIC[40])) + ALGEBRAIC[5];
STATES[9] = ((STATES[9] - ALGEBRAIC[6])*exp(-dt / ALGEBRAIC[41])) + ALGEBRAIC[6];
STATES[10] = ((STATES[10] - ALGEBRAIC[7])*exp(-dt / ALGEBRAIC[45])) + ALGEBRAIC[7];


STATES[16] += RATES[16] * dt;
STATES[12] += RATES[12] * dt;
STATES[2] += RATES[2] * dt;
STATES[0] += RATES[0] * dt;
STATES[1] += RATES[1] * dt;
STATES[3] += RATES[3] * dt;
STATES[15] += RATES[15] * dt;
}