/*
   There are a total of 75 entries in the algebraic variable array.
   There are a total of 21 entries in each of the rate and state variable arrays.
   There are a total of 49 entries in the constant variable array.
 */

#include "courtemanche_ramirez_nattel_1998.hpp"
#include <cmath>

/*
 * TIME is time in component environment (millisecond).
 * STATES[0] is V in component membrane (millivolt).
 * CONSTANTS[0] is R in component membrane (joule_per_mole_kelvin).
 * CONSTANTS[1] is T in component membrane (kelvin).
 * CONSTANTS[2] is F in component membrane (coulomb_per_millimole).
 * CONSTANTS[3] is Cm in component membrane (picoF).
 * ALGEBRAIC[0] is i_st in component membrane (picoA).
 * ALGEBRAIC[30] is i_Na in component fast_sodium_current (picoA).
 * ALGEBRAIC[50] is i_K1 in component time_independent_potassium_current (picoA).
 * ALGEBRAIC[51] is i_to in component transient_outward_K_current (picoA).
 * ALGEBRAIC[53] is i_Kur in component ultrarapid_delayed_rectifier_K_current (picoA).
 * ALGEBRAIC[54] is i_Kr in component rapid_delayed_rectifier_K_current (picoA).
 * ALGEBRAIC[55] is i_Ks in component slow_delayed_rectifier_K_current (picoA).
 * ALGEBRAIC[56] is i_CaL in component L_type_Ca_channel (picoA).
 * ALGEBRAIC[64] is i_CaP in component sarcolemmal_calcium_pump_current (picoA).
 * ALGEBRAIC[58] is i_NaK in component sodium_potassium_pump (picoA).
 * ALGEBRAIC[63] is i_NaCa in component Na_Ca_exchanger_current (picoA).
 * ALGEBRAIC[61] is i_B_Na in component background_currents (picoA).
 * ALGEBRAIC[62] is i_B_Ca in component background_currents (picoA).
 * CONSTANTS[4] is stim_start in component membrane (millisecond).
 * CONSTANTS[5] is stim_end in component membrane (millisecond).
 * CONSTANTS[6] is stim_period in component membrane (millisecond).
 * CONSTANTS[7] is stim_duration in component membrane (millisecond).
 * CONSTANTS[8] is stim_amplitude in component membrane (picoA).
 * ALGEBRAIC[17] is E_Na in component fast_sodium_current (millivolt).
 * CONSTANTS[9] is g_Na in component fast_sodium_current (nanoS_per_picoF).
 * STATES[1] is Na_i in component intracellular_ion_concentrations (millimolar).
 * CONSTANTS[10] is Na_o in component standard_ionic_concentrations (millimolar).
 * STATES[2] is m in component fast_sodium_current_m_gate (dimensionless).
 * STATES[3] is h in component fast_sodium_current_h_gate (dimensionless).
 * STATES[4] is j in component fast_sodium_current_j_gate (dimensionless).
 * ALGEBRAIC[1] is alpha_m in component fast_sodium_current_m_gate (per_millisecond).
 * ALGEBRAIC[18] is beta_m in component fast_sodium_current_m_gate (per_millisecond).
 * ALGEBRAIC[31] is m_inf in component fast_sodium_current_m_gate (dimensionless).
 * ALGEBRAIC[41] is tau_m in component fast_sodium_current_m_gate (millisecond).
 * ALGEBRAIC[2] is alpha_h in component fast_sodium_current_h_gate (per_millisecond).
 * ALGEBRAIC[19] is beta_h in component fast_sodium_current_h_gate (per_millisecond).
 * ALGEBRAIC[32] is h_inf in component fast_sodium_current_h_gate (dimensionless).
 * ALGEBRAIC[42] is tau_h in component fast_sodium_current_h_gate (millisecond).
 * ALGEBRAIC[3] is alpha_j in component fast_sodium_current_j_gate (per_millisecond).
 * ALGEBRAIC[20] is beta_j in component fast_sodium_current_j_gate (per_millisecond).
 * ALGEBRAIC[33] is j_inf in component fast_sodium_current_j_gate (dimensionless).
 * ALGEBRAIC[43] is tau_j in component fast_sodium_current_j_gate (millisecond).
 * ALGEBRAIC[40] is E_K in component time_independent_potassium_current (millivolt).
 * CONSTANTS[11] is g_K1 in component time_independent_potassium_current (nanoS_per_picoF).
 * CONSTANTS[12] is K_o in component standard_ionic_concentrations (millimolar).
 * STATES[5] is K_i in component intracellular_ion_concentrations (millimolar).
 * CONSTANTS[13] is K_Q10 in component transient_outward_K_current (dimensionless).
 * CONSTANTS[14] is g_to in component transient_outward_K_current (nanoS_per_picoF).
 * STATES[6] is oa in component transient_outward_K_current_oa_gate (dimensionless).
 * STATES[7] is oi in component transient_outward_K_current_oi_gate (dimensionless).
 * ALGEBRAIC[4] is alpha_oa in component transient_outward_K_current_oa_gate (per_millisecond).
 * ALGEBRAIC[21] is beta_oa in component transient_outward_K_current_oa_gate (per_millisecond).
 * ALGEBRAIC[34] is tau_oa in component transient_outward_K_current_oa_gate (millisecond).
 * ALGEBRAIC[44] is oa_infinity in component transient_outward_K_current_oa_gate (dimensionless).
 * ALGEBRAIC[5] is alpha_oi in component transient_outward_K_current_oi_gate (per_millisecond).
 * ALGEBRAIC[22] is beta_oi in component transient_outward_K_current_oi_gate (per_millisecond).
 * ALGEBRAIC[35] is tau_oi in component transient_outward_K_current_oi_gate (millisecond).
 * ALGEBRAIC[45] is oi_infinity in component transient_outward_K_current_oi_gate (dimensionless).
 * ALGEBRAIC[52] is g_Kur in component ultrarapid_delayed_rectifier_K_current (nanoS_per_picoF).
 * STATES[8] is ua in component ultrarapid_delayed_rectifier_K_current_ua_gate (dimensionless).
 * STATES[9] is ui in component ultrarapid_delayed_rectifier_K_current_ui_gate (dimensionless).
 * ALGEBRAIC[6] is alpha_ua in component ultrarapid_delayed_rectifier_K_current_ua_gate (per_millisecond).
 * ALGEBRAIC[23] is beta_ua in component ultrarapid_delayed_rectifier_K_current_ua_gate (per_millisecond).
 * ALGEBRAIC[36] is tau_ua in component ultrarapid_delayed_rectifier_K_current_ua_gate (millisecond).
 * ALGEBRAIC[46] is ua_infinity in component ultrarapid_delayed_rectifier_K_current_ua_gate (dimensionless).
 * ALGEBRAIC[7] is alpha_ui in component ultrarapid_delayed_rectifier_K_current_ui_gate (per_millisecond).
 * ALGEBRAIC[24] is beta_ui in component ultrarapid_delayed_rectifier_K_current_ui_gate (per_millisecond).
 * ALGEBRAIC[37] is tau_ui in component ultrarapid_delayed_rectifier_K_current_ui_gate (millisecond).
 * ALGEBRAIC[47] is ui_infinity in component ultrarapid_delayed_rectifier_K_current_ui_gate (dimensionless).
 * CONSTANTS[15] is g_Kr in component rapid_delayed_rectifier_K_current (nanoS_per_picoF).
 * STATES[10] is xr in component rapid_delayed_rectifier_K_current_xr_gate (dimensionless).
 * ALGEBRAIC[8] is alpha_xr in component rapid_delayed_rectifier_K_current_xr_gate (per_millisecond).
 * ALGEBRAIC[25] is beta_xr in component rapid_delayed_rectifier_K_current_xr_gate (per_millisecond).
 * ALGEBRAIC[38] is tau_xr in component rapid_delayed_rectifier_K_current_xr_gate (millisecond).
 * ALGEBRAIC[48] is xr_infinity in component rapid_delayed_rectifier_K_current_xr_gate (dimensionless).
 * CONSTANTS[16] is g_Ks in component slow_delayed_rectifier_K_current (nanoS_per_picoF).
 * STATES[11] is xs in component slow_delayed_rectifier_K_current_xs_gate (dimensionless).
 * ALGEBRAIC[9] is alpha_xs in component slow_delayed_rectifier_K_current_xs_gate (per_millisecond).
 * ALGEBRAIC[26] is beta_xs in component slow_delayed_rectifier_K_current_xs_gate (per_millisecond).
 * ALGEBRAIC[39] is tau_xs in component slow_delayed_rectifier_K_current_xs_gate (millisecond).
 * ALGEBRAIC[49] is xs_infinity in component slow_delayed_rectifier_K_current_xs_gate (dimensionless).
 * CONSTANTS[17] is g_CaL in component L_type_Ca_channel (nanoS_per_picoF).
 * STATES[12] is Ca_i in component intracellular_ion_concentrations (millimolar).
 * STATES[13] is d in component L_type_Ca_channel_d_gate (dimensionless).
 * STATES[14] is f in component L_type_Ca_channel_f_gate (dimensionless).
 * STATES[15] is f_Ca in component L_type_Ca_channel_f_Ca_gate (dimensionless).
 * ALGEBRAIC[10] is d_infinity in component L_type_Ca_channel_d_gate (dimensionless).
 * ALGEBRAIC[27] is tau_d in component L_type_Ca_channel_d_gate (millisecond).
 * ALGEBRAIC[11] is f_infinity in component L_type_Ca_channel_f_gate (dimensionless).
 * ALGEBRAIC[28] is tau_f in component L_type_Ca_channel_f_gate (millisecond).
 * ALGEBRAIC[12] is f_Ca_infinity in component L_type_Ca_channel_f_Ca_gate (dimensionless).
 * CONSTANTS[44] is tau_f_Ca in component L_type_Ca_channel_f_Ca_gate (millisecond).
 * CONSTANTS[18] is Km_Na_i in component sodium_potassium_pump (millimolar).
 * CONSTANTS[19] is Km_K_o in component sodium_potassium_pump (millimolar).
 * CONSTANTS[20] is i_NaK_max in component sodium_potassium_pump (picoA_per_picoF).
 * ALGEBRAIC[57] is f_NaK in component sodium_potassium_pump (dimensionless).
 * CONSTANTS[45] is sigma in component sodium_potassium_pump (dimensionless).
 * ALGEBRAIC[60] is i_B_K in component background_currents (picoA).
 * CONSTANTS[21] is g_B_Na in component background_currents (nanoS_per_picoF).
 * CONSTANTS[22] is g_B_Ca in component background_currents (nanoS_per_picoF).
 * CONSTANTS[23] is g_B_K in component background_currents (nanoS_per_picoF).
 * ALGEBRAIC[59] is E_Ca in component background_currents (millivolt).
 * CONSTANTS[24] is Ca_o in component standard_ionic_concentrations (millimolar).
 * CONSTANTS[25] is I_NaCa_max in component Na_Ca_exchanger_current (picoA_per_picoF).
 * CONSTANTS[26] is K_mNa in component Na_Ca_exchanger_current (millimolar).
 * CONSTANTS[27] is K_mCa in component Na_Ca_exchanger_current (millimolar).
 * CONSTANTS[28] is K_sat in component Na_Ca_exchanger_current (dimensionless).
 * CONSTANTS[29] is gmma in component Na_Ca_exchanger_current (dimensionless).
 * CONSTANTS[30] is i_CaP_max in component sarcolemmal_calcium_pump_current (picoA_per_picoF).
 * ALGEBRAIC[65] is i_rel in component Ca_release_current_from_JSR (millimolar_per_millisecond).
 * ALGEBRAIC[66] is Fn in component Ca_release_current_from_JSR (dimensionless).
 * CONSTANTS[31] is K_rel in component Ca_release_current_from_JSR (per_millisecond).
 * CONSTANTS[47] is V_rel in component intracellular_ion_concentrations (micrometre_3).
 * STATES[16] is Ca_rel in component intracellular_ion_concentrations (millimolar).
 * STATES[17] is u in component Ca_release_current_from_JSR_u_gate (dimensionless).
 * STATES[18] is v in component Ca_release_current_from_JSR_v_gate (dimensionless).
 * STATES[19] is w in component Ca_release_current_from_JSR_w_gate (dimensionless).
 * CONSTANTS[46] is tau_u in component Ca_release_current_from_JSR_u_gate (millisecond).
 * ALGEBRAIC[68] is u_infinity in component Ca_release_current_from_JSR_u_gate (dimensionless).
 * ALGEBRAIC[69] is tau_v in component Ca_release_current_from_JSR_v_gate (millisecond).
 * ALGEBRAIC[71] is v_infinity in component Ca_release_current_from_JSR_v_gate (dimensionless).
 * ALGEBRAIC[13] is tau_w in component Ca_release_current_from_JSR_w_gate (millisecond).
 * ALGEBRAIC[29] is w_infinity in component Ca_release_current_from_JSR_w_gate (dimensionless).
 * ALGEBRAIC[67] is i_tr in component transfer_current_from_NSR_to_JSR (millimolar_per_millisecond).
 * CONSTANTS[32] is tau_tr in component transfer_current_from_NSR_to_JSR (millisecond).
 * STATES[20] is Ca_up in component intracellular_ion_concentrations (millimolar).
 * CONSTANTS[33] is I_up_max in component Ca_uptake_current_by_the_NSR (millimolar_per_millisecond).
 * ALGEBRAIC[70] is i_up in component Ca_uptake_current_by_the_NSR (millimolar_per_millisecond).
 * CONSTANTS[34] is K_up in component Ca_uptake_current_by_the_NSR (millimolar).
 * ALGEBRAIC[72] is i_up_leak in component Ca_leak_current_by_the_NSR (millimolar_per_millisecond).
 * CONSTANTS[35] is Ca_up_max in component Ca_leak_current_by_the_NSR (millimolar).
 * CONSTANTS[36] is CMDN_max in component Ca_buffers (millimolar).
 * CONSTANTS[37] is TRPN_max in component Ca_buffers (millimolar).
 * CONSTANTS[38] is CSQN_max in component Ca_buffers (millimolar).
 * CONSTANTS[39] is Km_CMDN in component Ca_buffers (millimolar).
 * CONSTANTS[40] is Km_TRPN in component Ca_buffers (millimolar).
 * CONSTANTS[41] is Km_CSQN in component Ca_buffers (millimolar).
 * ALGEBRAIC[14] is Ca_CMDN in component Ca_buffers (millimolar).
 * ALGEBRAIC[15] is Ca_TRPN in component Ca_buffers (millimolar).
 * ALGEBRAIC[16] is Ca_CSQN in component Ca_buffers (millimolar).
 * CONSTANTS[42] is V_cell in component intracellular_ion_concentrations (micrometre_3).
 * CONSTANTS[43] is V_i in component intracellular_ion_concentrations (micrometre_3).
 * CONSTANTS[48] is V_up in component intracellular_ion_concentrations (micrometre_3).
 * ALGEBRAIC[73] is B1 in component intracellular_ion_concentrations (millimolar_per_millisecond).
 * ALGEBRAIC[74] is B2 in component intracellular_ion_concentrations (dimensionless).
 * RATES[0] is d/dt V in component membrane (millivolt).
 * RATES[2] is d/dt m in component fast_sodium_current_m_gate (dimensionless).
 * RATES[3] is d/dt h in component fast_sodium_current_h_gate (dimensionless).
 * RATES[4] is d/dt j in component fast_sodium_current_j_gate (dimensionless).
 * RATES[6] is d/dt oa in component transient_outward_K_current_oa_gate (dimensionless).
 * RATES[7] is d/dt oi in component transient_outward_K_current_oi_gate (dimensionless).
 * RATES[8] is d/dt ua in component ultrarapid_delayed_rectifier_K_current_ua_gate (dimensionless).
 * RATES[9] is d/dt ui in component ultrarapid_delayed_rectifier_K_current_ui_gate (dimensionless).
 * RATES[10] is d/dt xr in component rapid_delayed_rectifier_K_current_xr_gate (dimensionless).
 * RATES[11] is d/dt xs in component slow_delayed_rectifier_K_current_xs_gate (dimensionless).
 * RATES[13] is d/dt d in component L_type_Ca_channel_d_gate (dimensionless).
 * RATES[14] is d/dt f in component L_type_Ca_channel_f_gate (dimensionless).
 * RATES[15] is d/dt f_Ca in component L_type_Ca_channel_f_Ca_gate (dimensionless).
 * RATES[17] is d/dt u in component Ca_release_current_from_JSR_u_gate (dimensionless).
 * RATES[18] is d/dt v in component Ca_release_current_from_JSR_v_gate (dimensionless).
 * RATES[19] is d/dt w in component Ca_release_current_from_JSR_w_gate (dimensionless).
 * RATES[1] is d/dt Na_i in component intracellular_ion_concentrations (millimolar).
 * RATES[5] is d/dt K_i in component intracellular_ion_concentrations (millimolar).
 * RATES[12] is d/dt Ca_i in component intracellular_ion_concentrations (millimolar).
 * RATES[20] is d/dt Ca_up in component intracellular_ion_concentrations (millimolar).
 * RATES[16] is d/dt Ca_rel in component intracellular_ion_concentrations (millimolar).
 */


courtemanche_ramirez_nattel_1998::courtemanche_ramirez_nattel_1998()
{
algebraic_size = 75;
constants_size = 49;
states_size = 21;
ALGEBRAIC = new double[algebraic_size];
CONSTANTS = new double[constants_size];
RATES = new double[states_size];
STATES = new double[states_size];
}

courtemanche_ramirez_nattel_1998::~courtemanche_ramirez_nattel_1998()
{
delete mutation;
delete []ALGEBRAIC;
delete []CONSTANTS;
delete []RATES;
delete []STATES;
}

void courtemanche_ramirez_nattel_1998::initConsts( )
{
STATES[0] = -81.18;
CONSTANTS[0] = 8.3143;
CONSTANTS[1] = 310;
CONSTANTS[2] = 96.4867;
CONSTANTS[3] = 100;
CONSTANTS[4] = 5;
CONSTANTS[5] = 50000;
CONSTANTS[6] = 1000;
CONSTANTS[7] = 2;
CONSTANTS[8] = -2000;
CONSTANTS[9] = 7.8;
STATES[1] = 1.117e+01;
CONSTANTS[10] = 140;
STATES[2] = 2.908e-3;
STATES[3] = 9.649e-1;
STATES[4] = 9.775e-1;
CONSTANTS[11] = 0.09;
CONSTANTS[12] = 5.4;
STATES[5] = 1.39e+02;
CONSTANTS[13] = 3;
CONSTANTS[14] = 0.1652;
STATES[6] = 3.043e-2;
STATES[7] = 9.992e-1;
STATES[8] = 4.966e-3;
STATES[9] = 9.986e-1;
CONSTANTS[15] = 0.029411765;
STATES[10] = 3.296e-5;
CONSTANTS[16] = 0.12941176;
STATES[11] = 1.869e-2;
CONSTANTS[17] = 0.12375;
STATES[12] = 1.013e-4;
STATES[13] = 1.367e-4;
STATES[14] = 9.996e-1;
STATES[15] = 7.755e-1;
CONSTANTS[18] = 10;
CONSTANTS[19] = 1.5;
CONSTANTS[20] = 0.59933874;
CONSTANTS[21] = 0.0006744375;
CONSTANTS[22] = 0.001131;
CONSTANTS[23] = 0;
CONSTANTS[24] = 1.8;
CONSTANTS[25] = 1600;
CONSTANTS[26] = 87.5;
CONSTANTS[27] = 1.38;
CONSTANTS[28] = 0.1;
CONSTANTS[29] = 0.35;
CONSTANTS[30] = 0.275;
CONSTANTS[31] = 30;
STATES[16] = 1.488;
STATES[17] = 2.35e-112;
STATES[18] = 1;
STATES[19] = 0.9992;
CONSTANTS[32] = 180;
STATES[20] = 1.488;
CONSTANTS[33] = 0.005;
CONSTANTS[34] = 0.00092;
CONSTANTS[35] = 15;
CONSTANTS[36] = 0.05;
CONSTANTS[37] = 0.07;
CONSTANTS[38] = 10;
CONSTANTS[39] = 0.00238;
CONSTANTS[40] = 0.0005;
CONSTANTS[41] = 0.8;
CONSTANTS[42] = 20100;
CONSTANTS[43] =  CONSTANTS[42]*0.680000;
CONSTANTS[44] = 2.00000;
CONSTANTS[45] =  (1.00000/7.00000)*(exp(CONSTANTS[10]/67.3000) - 1.00000);
CONSTANTS[46] = 8.00000;
CONSTANTS[47] =  0.00480000*CONSTANTS[42];
CONSTANTS[48] =  0.0552000*CONSTANTS[42];
}

void courtemanche_ramirez_nattel_1998::computeRates( double TIME, double *CONSTANTS, double *RATES, double *STATES, double *ALGEBRAIC )
{
ALGEBRAIC[12] = pow(1.00000+STATES[12]/0.000350000, - 1.00000);
ALGEBRAIC[10] = pow(1.00000+exp((STATES[0]+10.0000)/- 8.00000), - 1.00000);
ALGEBRAIC[27] = (fabs(STATES[0]+10.0000)<1.00000e-10 ? 4.57900/(1.00000+exp((STATES[0]+10.0000)/- 6.24000)) : (1.00000 - exp((STATES[0]+10.0000)/- 6.24000))/( 0.0350000*(STATES[0]+10.0000)*(1.00000+exp((STATES[0]+10.0000)/- 6.24000))));
ALGEBRAIC[11] = exp(- (STATES[0]+28.0000)/6.90000)/(1.00000+exp(- (STATES[0]+28.0000)/6.90000));
ALGEBRAIC[28] =  9.00000*pow( 0.0197000*exp( - pow(0.0337000, 2.00000)*pow(STATES[0]+10.0000, 2.00000))+0.0200000, - 1.00000);
ALGEBRAIC[13] = (fabs(STATES[0] - 7.90000)<1.00000e-10 ? ( 6.00000*0.200000)/1.30000 : ( 6.00000*(1.00000 - exp(- (STATES[0] - 7.90000)/5.00000)))/( (1.00000+ 0.300000*exp(- (STATES[0] - 7.90000)/5.00000))*1.00000*(STATES[0] - 7.90000)));
ALGEBRAIC[29] = 1.00000 - pow(1.00000+exp(- (STATES[0] - 40.0000)/17.0000), - 1.00000);
ALGEBRAIC[1] = (STATES[0]==- 47.1300 ? 3.20000 : ( 0.320000*(STATES[0]+47.1300))/(1.00000 - exp( - 0.100000*(STATES[0]+47.1300))));
ALGEBRAIC[18] =  0.0800000*exp(- STATES[0]/11.0000);
ALGEBRAIC[31] = ALGEBRAIC[1]/(ALGEBRAIC[1]+ALGEBRAIC[18]);
ALGEBRAIC[41] = 1.00000/(ALGEBRAIC[1]+ALGEBRAIC[18]);
ALGEBRAIC[2] = (STATES[0]<- 40.0000 ?  0.135000*exp((STATES[0]+80.0000)/- 6.80000) : 0.000000);
ALGEBRAIC[19] = (STATES[0]<- 40.0000 ?  3.56000*exp( 0.0790000*STATES[0])+ 310000.*exp( 0.350000*STATES[0]) : 1.00000/( 0.130000*(1.00000+exp((STATES[0]+10.6600)/- 11.1000))));
ALGEBRAIC[32] = ALGEBRAIC[2]/(ALGEBRAIC[2]+ALGEBRAIC[19]);
ALGEBRAIC[42] = 1.00000/(ALGEBRAIC[2]+ALGEBRAIC[19]);
ALGEBRAIC[3] = (STATES[0]<- 40.0000 ? ( ( - 127140.*exp( 0.244400*STATES[0]) -  3.47400e-05*exp( - 0.0439100*STATES[0]))*(STATES[0]+37.7800))/(1.00000+exp( 0.311000*(STATES[0]+79.2300))) : 0.000000);
ALGEBRAIC[20] = (STATES[0]<- 40.0000 ? ( 0.121200*exp( - 0.0105200*STATES[0]))/(1.00000+exp( - 0.137800*(STATES[0]+40.1400))) : ( 0.300000*exp( - 2.53500e-07*STATES[0]))/(1.00000+exp( - 0.100000*(STATES[0]+32.0000))));
ALGEBRAIC[33] = ALGEBRAIC[3]/(ALGEBRAIC[3]+ALGEBRAIC[20]);
ALGEBRAIC[43] = 1.00000/(ALGEBRAIC[3]+ALGEBRAIC[20]);
ALGEBRAIC[4] =  0.650000*pow(exp((STATES[0] - - 10.0000)/- 8.50000)+exp(((STATES[0] - - 10.0000) - 40.0000)/- 59.0000), - 1.00000);
ALGEBRAIC[21] =  0.650000*pow(2.50000+exp(((STATES[0] - - 10.0000)+72.0000)/17.0000), - 1.00000);
ALGEBRAIC[34] = pow(ALGEBRAIC[4]+ALGEBRAIC[21], - 1.00000)/CONSTANTS[13];
ALGEBRAIC[44] = pow(1.00000+exp(((STATES[0] - - 10.0000)+10.4700)/- 17.5400), - 1.00000);
ALGEBRAIC[5] = pow(18.5300+ 1.00000*exp(((STATES[0] - - 10.0000)+103.700)/10.9500), - 1.00000);
ALGEBRAIC[22] = pow(35.5600+ 1.00000*exp(((STATES[0] - - 10.0000) - 8.74000)/- 7.44000), - 1.00000);
ALGEBRAIC[35] = pow(ALGEBRAIC[5]+ALGEBRAIC[22], - 1.00000)/CONSTANTS[13];
ALGEBRAIC[45] = pow(1.00000+exp(((STATES[0] - - 10.0000)+33.1000)/5.30000), - 1.00000);
ALGEBRAIC[6] =  0.650000*pow(exp((STATES[0] - - 10.0000)/- 8.50000)+exp(((STATES[0] - - 10.0000) - 40.0000)/- 59.0000), - 1.00000);
ALGEBRAIC[23] =  0.650000*pow(2.50000+exp(((STATES[0] - - 10.0000)+72.0000)/17.0000), - 1.00000);
ALGEBRAIC[36] = pow(ALGEBRAIC[6]+ALGEBRAIC[23], - 1.00000)/CONSTANTS[13];
ALGEBRAIC[46] = pow(1.00000+exp(((STATES[0] - - 10.0000)+20.3000)/- 9.60000), - 1.00000);
ALGEBRAIC[7] = pow(21.0000+ 1.00000*exp(((STATES[0] - - 10.0000) - 195.000)/- 28.0000), - 1.00000);
ALGEBRAIC[24] = 1.00000/exp(((STATES[0] - - 10.0000) - 168.000)/- 16.0000);
ALGEBRAIC[37] = pow(ALGEBRAIC[7]+ALGEBRAIC[24], - 1.00000)/CONSTANTS[13];
ALGEBRAIC[47] = pow(1.00000+exp(((STATES[0] - - 10.0000) - 109.450)/27.4800), - 1.00000);
ALGEBRAIC[8] = (fabs(STATES[0]+14.1000)<1.00000e-10 ? 0.00150000 : ( 0.000300000*(STATES[0]+14.1000))/(1.00000 - exp((STATES[0]+14.1000)/- 5.00000)));
ALGEBRAIC[25] = (fabs(STATES[0] - 3.33280)<1.00000e-10 ? 0.000378361 : ( 7.38980e-05*(STATES[0] - 3.33280))/(exp((STATES[0] - 3.33280)/5.12370) - 1.00000));
ALGEBRAIC[38] = pow(ALGEBRAIC[8]+ALGEBRAIC[25], - 1.00000);
ALGEBRAIC[48] = pow(1.00000+exp((STATES[0]+14.1000)/- 6.50000), - 1.00000);
ALGEBRAIC[9] = (fabs(STATES[0] - 19.9000)<1.00000e-10 ? 0.000680000 : ( 4.00000e-05*(STATES[0] - 19.9000))/(1.00000 - exp((STATES[0] - 19.9000)/- 17.0000)));
ALGEBRAIC[26] = (fabs(STATES[0] - 19.9000)<1.00000e-10 ? 0.000315000 : ( 3.50000e-05*(STATES[0] - 19.9000))/(exp((STATES[0] - 19.9000)/9.00000) - 1.00000));
ALGEBRAIC[39] =  0.500000*pow(ALGEBRAIC[9]+ALGEBRAIC[26], - 1.00000);
ALGEBRAIC[49] = pow(1.00000+exp((STATES[0] - 19.9000)/- 12.7000), - 0.500000);
ALGEBRAIC[40] =  (( CONSTANTS[0]*CONSTANTS[1])/CONSTANTS[2])*log(CONSTANTS[12]/STATES[5]);
ALGEBRAIC[50] = ( CONSTANTS[3]*CONSTANTS[11]*(STATES[0] - ALGEBRAIC[40]))/(1.00000+exp( 0.0700000*(STATES[0]+80.0000)));
ALGEBRAIC[51] =  CONSTANTS[3]*CONSTANTS[14]*pow(STATES[6], 3.00000)*STATES[7]*(STATES[0] - ALGEBRAIC[40]);
ALGEBRAIC[52] = 0.00500000+0.0500000/(1.00000+exp((STATES[0] - 15.0000)/- 13.0000));
ALGEBRAIC[53] =  CONSTANTS[3]*ALGEBRAIC[52]*pow(STATES[8], 3.00000)*STATES[9]*(STATES[0] - ALGEBRAIC[40]);
ALGEBRAIC[54] = ( CONSTANTS[3]*CONSTANTS[15]*STATES[10]*(STATES[0] - ALGEBRAIC[40]))/(1.00000+exp((STATES[0]+15.0000)/22.4000));
ALGEBRAIC[55] =  CONSTANTS[3]*CONSTANTS[16]*pow(STATES[11], 2.00000)*(STATES[0] - ALGEBRAIC[40]);
ALGEBRAIC[57] = pow(1.00000+ 0.124500*exp(( - 0.100000*CONSTANTS[2]*STATES[0])/( CONSTANTS[0]*CONSTANTS[1]))+ 0.0365000*CONSTANTS[45]*exp(( - CONSTANTS[2]*STATES[0])/( CONSTANTS[0]*CONSTANTS[1])), - 1.00000);
ALGEBRAIC[58] = ( (( CONSTANTS[3]*CONSTANTS[20]*ALGEBRAIC[57]*1.00000)/(1.00000+pow(CONSTANTS[18]/STATES[1], 1.50000)))*CONSTANTS[12])/(CONSTANTS[12]+CONSTANTS[19]);
ALGEBRAIC[60] =  CONSTANTS[3]*CONSTANTS[23]*(STATES[0] - ALGEBRAIC[40]);
ALGEBRAIC[17] =  (( CONSTANTS[0]*CONSTANTS[1])/CONSTANTS[2])*log(CONSTANTS[10]/STATES[1]);
ALGEBRAIC[30] =  CONSTANTS[3]*CONSTANTS[9]*pow(STATES[2], 3.00000)*STATES[3]*STATES[4]*(STATES[0] - ALGEBRAIC[17]);
ALGEBRAIC[63] = ( CONSTANTS[3]*CONSTANTS[25]*( exp(( CONSTANTS[29]*CONSTANTS[2]*STATES[0])/( CONSTANTS[0]*CONSTANTS[1]))*pow(STATES[1], 3.00000)*CONSTANTS[24] -  exp(( (CONSTANTS[29] - 1.00000)*CONSTANTS[2]*STATES[0])/( CONSTANTS[0]*CONSTANTS[1]))*pow(CONSTANTS[10], 3.00000)*STATES[12]))/( (pow(CONSTANTS[26], 3.00000)+pow(CONSTANTS[10], 3.00000))*(CONSTANTS[27]+CONSTANTS[24])*(1.00000+ CONSTANTS[28]*exp(( (CONSTANTS[29] - 1.00000)*STATES[0]*CONSTANTS[2])/( CONSTANTS[0]*CONSTANTS[1]))));
ALGEBRAIC[61] =  CONSTANTS[3]*CONSTANTS[21]*(STATES[0] - ALGEBRAIC[17]);

#if defined(SINGLE_CELL)
ALGEBRAIC[0] = (TIME>=CONSTANTS[4]&&TIME<=CONSTANTS[5]&&(TIME - CONSTANTS[4]) -  floor((TIME - CONSTANTS[4])/CONSTANTS[6])*CONSTANTS[6]<=CONSTANTS[7] ? CONSTANTS[8] : 0.000000);
#else
if(isS1) ALGEBRAIC[0] = CONSTANTS[8];
else ALGEBRAIC[0] = 0.0;
if(isEctopic) ALGEBRAIC[0] = CONSTANTS[8];
else {
  if(!isS1) ALGEBRAIC[0] = 0.0;
}
#endif

ALGEBRAIC[56] =  CONSTANTS[3]*CONSTANTS[17]*STATES[13]*STATES[14]*STATES[15]*(STATES[0] - 65.0000);
ALGEBRAIC[64] = ( CONSTANTS[3]*CONSTANTS[30]*STATES[12])/(0.000500000+STATES[12]);
ALGEBRAIC[59] =  (( CONSTANTS[0]*CONSTANTS[1])/( 2.00000*CONSTANTS[2]))*log(CONSTANTS[24]/STATES[12]);
ALGEBRAIC[62] =  CONSTANTS[3]*CONSTANTS[22]*(STATES[0] - ALGEBRAIC[59]);
ALGEBRAIC[65] =  CONSTANTS[31]*pow(STATES[17], 2.00000)*STATES[18]*STATES[19]*(STATES[16] - STATES[12]);
ALGEBRAIC[67] = (STATES[20] - STATES[16])/CONSTANTS[32];
ALGEBRAIC[66] =  1000.00*( 1.00000e-15*CONSTANTS[47]*ALGEBRAIC[65] -  (1.00000e-15/( 2.00000*CONSTANTS[2]))*( 0.500000*ALGEBRAIC[56] -  0.200000*ALGEBRAIC[63]));
ALGEBRAIC[68] = pow(1.00000+exp(- (ALGEBRAIC[66] - 3.41750e-13)/1.36700e-15), - 1.00000);
ALGEBRAIC[69] = 1.91000+ 2.09000*pow(1.00000+exp(- (ALGEBRAIC[66] - 3.41750e-13)/1.36700e-15), - 1.00000);
ALGEBRAIC[71] = 1.00000 - pow(1.00000+exp(- (ALGEBRAIC[66] - 6.83500e-14)/1.36700e-15), - 1.00000);
ALGEBRAIC[70] = CONSTANTS[33]/(1.00000+CONSTANTS[34]/STATES[12]);
ALGEBRAIC[72] = ( CONSTANTS[33]*STATES[20])/CONSTANTS[35];
ALGEBRAIC[73] = ( 2.00000*ALGEBRAIC[63] - (ALGEBRAIC[64]+ALGEBRAIC[56]+ALGEBRAIC[62]))/( 2.00000*CONSTANTS[43]*CONSTANTS[2])+( CONSTANTS[48]*(ALGEBRAIC[72] - ALGEBRAIC[70])+ ALGEBRAIC[65]*CONSTANTS[47])/CONSTANTS[43];
ALGEBRAIC[74] = 1.00000+( CONSTANTS[37]*CONSTANTS[40])/pow(STATES[12]+CONSTANTS[40], 2.00000)+( CONSTANTS[36]*CONSTANTS[39])/pow(STATES[12]+CONSTANTS[39], 2.00000);
ALGEBRAIC[14] = ( CONSTANTS[36]*STATES[12])/(STATES[12]+CONSTANTS[39]);
ALGEBRAIC[15] = ( CONSTANTS[37]*STATES[12])/(STATES[12]+CONSTANTS[40]);
ALGEBRAIC[16] = ( CONSTANTS[38]*STATES[16])/(STATES[16]+CONSTANTS[41]);

if(mutation != 0 && isMutated == true){
  mutation->mutate( ALGEBRAIC, STATES, CONSTANTS );
}

RATES[15] = (ALGEBRAIC[12] - STATES[15])/CONSTANTS[44];
RATES[13] = (ALGEBRAIC[10] - STATES[13])/ALGEBRAIC[27];
RATES[14] = (ALGEBRAIC[11] - STATES[14])/ALGEBRAIC[28];
RATES[19] = (ALGEBRAIC[29] - STATES[19])/ALGEBRAIC[13];
RATES[2] = (ALGEBRAIC[31] - STATES[2])/ALGEBRAIC[41];
RATES[3] = (ALGEBRAIC[32] - STATES[3])/ALGEBRAIC[42];
RATES[4] = (ALGEBRAIC[33] - STATES[4])/ALGEBRAIC[43];
RATES[6] = (ALGEBRAIC[44] - STATES[6])/ALGEBRAIC[34];
RATES[7] = (ALGEBRAIC[45] - STATES[7])/ALGEBRAIC[35];
RATES[8] = (ALGEBRAIC[46] - STATES[8])/ALGEBRAIC[36];
RATES[9] = (ALGEBRAIC[47] - STATES[9])/ALGEBRAIC[37];
RATES[10] = (ALGEBRAIC[48] - STATES[10])/ALGEBRAIC[38];
RATES[11] = (ALGEBRAIC[49] - STATES[11])/ALGEBRAIC[39];
RATES[5] = ( 2.00000*ALGEBRAIC[58] - (ALGEBRAIC[50]+ALGEBRAIC[51]+ALGEBRAIC[53]+ALGEBRAIC[54]+ALGEBRAIC[55]+ALGEBRAIC[60]))/( CONSTANTS[43]*CONSTANTS[2]);
RATES[1] = ( - 3.00000*ALGEBRAIC[58] - ( 3.00000*ALGEBRAIC[63]+ALGEBRAIC[61]+ALGEBRAIC[30]))/( CONSTANTS[43]*CONSTANTS[2]);
RATES[0] = - (ALGEBRAIC[30]+ALGEBRAIC[50]+ALGEBRAIC[51]+ALGEBRAIC[53]+ALGEBRAIC[54]+ALGEBRAIC[55]+ALGEBRAIC[61]+ALGEBRAIC[62]+ALGEBRAIC[58]+ALGEBRAIC[64]+ALGEBRAIC[63]+ALGEBRAIC[56]+ALGEBRAIC[0])/CONSTANTS[3];
RATES[16] =  (ALGEBRAIC[67] - ALGEBRAIC[65])*pow(1.00000+( CONSTANTS[38]*CONSTANTS[41])/pow(STATES[16]+CONSTANTS[41], 2.00000), - 1.00000);
RATES[17] = (ALGEBRAIC[68] - STATES[17])/CONSTANTS[46];
RATES[18] = (ALGEBRAIC[71] - STATES[18])/ALGEBRAIC[69];
RATES[20] = ALGEBRAIC[70] - (ALGEBRAIC[72]+( ALGEBRAIC[67]*CONSTANTS[47])/CONSTANTS[48]);
RATES[12] = ALGEBRAIC[73]/ALGEBRAIC[74];
}

void courtemanche_ramirez_nattel_1998::solveAnalytical(double dt)
{
STATES[15] = ((STATES[15] - ALGEBRAIC[12])*exp(-dt / ALGEBRAIC[44])) + ALGEBRAIC[12];
STATES[13] = ((STATES[13] - ALGEBRAIC[10])*exp(-dt / ALGEBRAIC[27])) + ALGEBRAIC[10];
STATES[14] = ((STATES[14] - ALGEBRAIC[11])*exp(-dt / ALGEBRAIC[28])) + ALGEBRAIC[11];
STATES[19] = ((STATES[19] - ALGEBRAIC[29])*exp(-dt / ALGEBRAIC[13])) + ALGEBRAIC[29];
STATES[2] = ((STATES[2] - ALGEBRAIC[31])*exp(-dt / ALGEBRAIC[41])) + ALGEBRAIC[31];
STATES[3] = ((STATES[3] - ALGEBRAIC[32])*exp(-dt / ALGEBRAIC[42])) + ALGEBRAIC[32];
STATES[4] = ((STATES[4] - ALGEBRAIC[33])*exp(-dt / ALGEBRAIC[43])) + ALGEBRAIC[33];
STATES[6] = ((STATES[6] - ALGEBRAIC[44])*exp(-dt / ALGEBRAIC[34])) + ALGEBRAIC[44];
STATES[7] = ((STATES[7] - ALGEBRAIC[45])*exp(-dt / ALGEBRAIC[35])) + ALGEBRAIC[45];
STATES[8] = ((STATES[8] - ALGEBRAIC[46])*exp(-dt / ALGEBRAIC[36])) + ALGEBRAIC[46];
STATES[9] = ((STATES[9] - ALGEBRAIC[47])*exp(-dt / ALGEBRAIC[37])) + ALGEBRAIC[47];
STATES[10] = ((STATES[10] - ALGEBRAIC[48])*exp(-dt / ALGEBRAIC[38])) + ALGEBRAIC[48];
STATES[11] = ((STATES[11] - ALGEBRAIC[49])*exp(-dt / ALGEBRAIC[39])) + ALGEBRAIC[49];
STATES[17] = ((STATES[17] - ALGEBRAIC[68])*exp(-dt / ALGEBRAIC[46])) + ALGEBRAIC[68];
STATES[18] = ((STATES[18] - ALGEBRAIC[71])*exp(-dt / ALGEBRAIC[69])) + ALGEBRAIC[71];
STATES[0] = STATES[0] + ( (-dt/CONSTANTS[3]) * (ALGEBRAIC[30]+ALGEBRAIC[50]+ALGEBRAIC[51]+ALGEBRAIC[53]+ALGEBRAIC[54]+ALGEBRAIC[55]+ALGEBRAIC[61]+ALGEBRAIC[62]+ALGEBRAIC[58]+ALGEBRAIC[64]+ALGEBRAIC[63]+ALGEBRAIC[56]+ALGEBRAIC[0]) );
STATES[1] = STATES[1] + RATES[1] * dt;
STATES[5] = STATES[5] + RATES[5] * dt;
STATES[16] = STATES[16] + RATES[16] * dt;
STATES[20] = STATES[20] + RATES[20] * dt;
STATES[12] = STATES[12] + RATES[12] * dt;
}
