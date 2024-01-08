/*
   There are a total of 198 entries in the algebraic variable array.
   There are a total of 41 entries in each of the rate and state variable arrays.
   There are a total of 139+2 entries in the constant variable array.
 */

#include "Ohara_Rudy_2011.hpp"
#include "../modules/globals.hpp"
#include <cmath>
#include <cstdlib>
#include <cstdio>

/*
 * TIME is time in component environment (millisecond).
 * CONSTANTS[0] is celltype in component environment (dimensionless).
 * CONSTANTS[1] is nao in component extracellular (millimolar).
 * CONSTANTS[2] is cao in component extracellular (millimolar).
 * CONSTANTS[3] is ko in component extracellular (millimolar).
 * CONSTANTS[4] is R in component physical_constants (joule_per_kilomole_kelvin).
 * CONSTANTS[5] is T in component physical_constants (kelvin).
 * CONSTANTS[6] is F in component physical_constants (coulomb_per_mole).
 * CONSTANTS[7] is zna in component physical_constants (dimensionless).
 * CONSTANTS[8] is zca in component physical_constants (dimensionless).
 * CONSTANTS[9] is zk in component physical_constants (dimensionless).
 * CONSTANTS[10] is L in component cell_geometry (centimeter).
 * CONSTANTS[11] is rad in component cell_geometry (centimeter).
 * CONSTANTS[104] is vcell in component cell_geometry (microliter).
 * CONSTANTS[113] is Ageo in component cell_geometry (centimeter_squared).
 * CONSTANTS[117] is Acap in component cell_geometry (centimeter_squared).
 * CONSTANTS[118] is vmyo in component cell_geometry (microliter).
 * CONSTANTS[119] is vnsr in component cell_geometry (microliter).
 * CONSTANTS[120] is vjsr in component cell_geometry (microliter).
 * CONSTANTS[121] is vss in component cell_geometry (microliter).
 * STATES[0] is v in component membrane (millivolt).
 * ALGEBRAIC[29] is vffrt in component membrane (coulomb_per_mole).
 * ALGEBRAIC[39] is vfrt in component membrane (dimensionless).
 * ALGEBRAIC[62] is INa in component INa (microA_per_microF).
 * ALGEBRAIC[64] is INaL in component INaL (microA_per_microF).
 * ALGEBRAIC[70] is Ito in component Ito (microA_per_microF).
 * ALGEBRAIC[81] is ICaL in component ICaL (microA_per_microF).
 * ALGEBRAIC[82] is ICaNa in component ICaL (microA_per_microF).
 * ALGEBRAIC[85] is ICaK in component ICaL (microA_per_microF).
 * ALGEBRAIC[96] is IKr in component IKr (microA_per_microF).
 * ALGEBRAIC[98] is IKs in component IKs (microA_per_microF).
 * ALGEBRAIC[100] is IK1 in component IK1 (microA_per_microF).
 * ALGEBRAIC[132] is INaCa_i in component INaCa_i (microA_per_microF).
 * ALGEBRAIC[162] is INaCa_ss in component INaCa_i (microA_per_microF).
 * ALGEBRAIC[181] is INaK in component INaK (microA_per_microF).
 * ALGEBRAIC[184] is INab in component INab (microA_per_microF).
 * ALGEBRAIC[183] is IKb in component IKb (microA_per_microF).
 * ALGEBRAIC[188] is IpCa in component IpCa (microA_per_microF).
 * ALGEBRAIC[186] is ICab in component ICab (microA_per_microF).
 * ALGEBRAIC[12] is Istim in component membrane (microA_per_microF).
 * CONSTANTS[12] is amp in component membrane (microA_per_microF).
 * CONSTANTS[13] is duration in component membrane (millisecond).
 * CONSTANTS[14] is KmCaMK in component CaMK (millimolar).
 * CONSTANTS[15] is aCaMK in component CaMK (per_millimolar_per_millisecond).
 * CONSTANTS[16] is bCaMK in component CaMK (per_millisecond).
 * CONSTANTS[17] is CaMKo in component CaMK (dimensionless).
 * CONSTANTS[18] is KmCaM in component CaMK (millimolar).
 * ALGEBRAIC[45] is CaMKb in component CaMK (millimolar).
 * ALGEBRAIC[47] is CaMKa in component CaMK (millimolar).
 * STATES[1] is CaMKt in component CaMK (millimolar).
 * STATES[2] is cass in component intracellular_ions (millimolar).
 * CONSTANTS[19] is cmdnmax_b in component intracellular_ions (millimolar).
 * CONSTANTS[93] is cmdnmax in component intracellular_ions (millimolar).
 * CONSTANTS[20] is kmcmdn in component intracellular_ions (millimolar).
 * CONSTANTS[21] is trpnmax in component intracellular_ions (millimolar).
 * CONSTANTS[22] is kmtrpn in component intracellular_ions (millimolar).
 * CONSTANTS[23] is BSRmax in component intracellular_ions (millimolar).
 * CONSTANTS[24] is KmBSR in component intracellular_ions (millimolar).
 * CONSTANTS[25] is BSLmax in component intracellular_ions (millimolar).
 * CONSTANTS[26] is KmBSL in component intracellular_ions (millimolar).
 * CONSTANTS[27] is csqnmax in component intracellular_ions (millimolar).
 * CONSTANTS[28] is kmcsqn in component intracellular_ions (millimolar).
 * STATES[3] is nai in component intracellular_ions (millimolar).
 * STATES[4] is nass in component intracellular_ions (millimolar).
 * STATES[5] is ki in component intracellular_ions (millimolar).
 * STATES[6] is kss in component intracellular_ions (millimolar).
 * STATES[7] is cansr in component intracellular_ions (millimolar).
 * STATES[8] is cajsr in component intracellular_ions (millimolar).
 * STATES[9] is cai in component intracellular_ions (millimolar).
 * ALGEBRAIC[187] is JdiffNa in component diff (millimolar_per_millisecond).
 * ALGEBRAIC[189] is Jdiff in component diff (millimolar_per_millisecond).
 * ALGEBRAIC[196] is Jup in component SERCA (millimolar_per_millisecond).
 * ALGEBRAIC[185] is JdiffK in component diff (millimolar_per_millisecond).
 * ALGEBRAIC[191] is Jrel in component ryr (millimolar_per_millisecond).
 * ALGEBRAIC[197] is Jtr in component trans_flux (millimolar_per_millisecond).
 * ALGEBRAIC[49] is Bcai in component intracellular_ions (dimensionless).
 * ALGEBRAIC[53] is Bcajsr in component intracellular_ions (dimensionless).
 * ALGEBRAIC[51] is Bcass in component intracellular_ions (dimensionless).
 * CONSTANTS[29] is cm in component intracellular_ions (microF_per_centimeter_squared).
 * CONSTANTS[30] is PKNa in component reversal_potentials (dimensionless).
 * ALGEBRAIC[56] is ENa in component reversal_potentials (millivolt).
 * ALGEBRAIC[57] is EK in component reversal_potentials (millivolt).
 * ALGEBRAIC[58] is EKs in component reversal_potentials (millivolt).
 * ALGEBRAIC[0] is mss in component INa (dimensionless).
 * ALGEBRAIC[13] is tm in component INa (millisecond).
 * CONSTANTS[31] is mssV1 in component INa (millivolt).
 * CONSTANTS[32] is mssV2 in component INa (millivolt).
 * CONSTANTS[33] is mtV1 in component INa (millivolt).
 * CONSTANTS[34] is mtV2 in component INa (millivolt).
 * CONSTANTS[35] is mtD1 in component INa (dimensionless).
 * CONSTANTS[36] is mtD2 in component INa (dimensionless).
 * CONSTANTS[37] is mtV3 in component INa (millivolt).
 * CONSTANTS[38] is mtV4 in component INa (millivolt).
 * STATES[10] is m in component INa (dimensionless).
 * ALGEBRAIC[1] is hss in component INa (dimensionless).
 * ALGEBRAIC[14] is thf in component INa (millisecond).
 * ALGEBRAIC[15] is ths in component INa (millisecond).
 * CONSTANTS[39] is hssV1 in component INa (millivolt).
 * CONSTANTS[40] is hssV2 in component INa (millivolt).
 * CONSTANTS[94] is Ahs in component INa (dimensionless).
 * CONSTANTS[41] is Ahf in component INa (dimensionless).
 * STATES[11] is hf in component INa (dimensionless).
 * STATES[12] is hs in component INa (dimensionless).
 * ALGEBRAIC[59] is h in component INa (dimensionless).
 * CONSTANTS[42] is GNa in component INa (milliS_per_microF).
 * ALGEBRAIC[16] is jss in component INa (dimensionless).
 * ALGEBRAIC[30] is tj in component INa (millisecond).
 * STATES[13] is j in component INa (dimensionless).
 * ALGEBRAIC[31] is hssp in component INa (dimensionless).
 * ALGEBRAIC[40] is thsp in component INa (millisecond).
 * STATES[14] is hsp in component INa (dimensionless).
 * ALGEBRAIC[60] is hp in component INa (dimensionless).
 * ALGEBRAIC[41] is tjp in component INa (millisecond).
 * STATES[15] is jp in component INa (dimensionless).
 * ALGEBRAIC[61] is fINap in component INa (dimensionless).
 * ALGEBRAIC[32] is mLss in component INaL (dimensionless).
 * ALGEBRAIC[42] is tmL in component INaL (millisecond).
 * STATES[16] is mL in component INaL (dimensionless).
 * CONSTANTS[43] is thL in component INaL (millisecond).
 * ALGEBRAIC[2] is hLss in component INaL (dimensionless).
 * STATES[17] is hL in component INaL (dimensionless).
 * ALGEBRAIC[3] is hLssp in component INaL (dimensionless).
 * CONSTANTS[95] is thLp in component INaL (millisecond).
 * STATES[18] is hLp in component INaL (dimensionless).
 * CONSTANTS[44] is GNaL_b in component INaL (milliS_per_microF).
 * CONSTANTS[96] is GNaL in component INaL (milliS_per_microF).
 * ALGEBRAIC[63] is fINaLp in component INaL (dimensionless).
 * CONSTANTS[45] is Gto_b in component Ito (milliS_per_microF).
 * ALGEBRAIC[4] is ass in component Ito (dimensionless).
 * ALGEBRAIC[17] is ta in component Ito (millisecond).
 * STATES[19] is a in component Ito (dimensionless).
 * ALGEBRAIC[5] is iss in component Ito (dimensionless).
 * ALGEBRAIC[18] is delta_epi in component Ito (dimensionless).
 * ALGEBRAIC[33] is tiF_b in component Ito (millisecond).
 * ALGEBRAIC[43] is tiS_b in component Ito (millisecond).
 * ALGEBRAIC[46] is tiF in component Ito (millisecond).
 * ALGEBRAIC[48] is tiS in component Ito (millisecond).
 * ALGEBRAIC[65] is AiF in component Ito (dimensionless).
 * ALGEBRAIC[66] is AiS in component Ito (dimensionless).
 * STATES[20] is iF in component Ito (dimensionless).
 * STATES[21] is iS in component Ito (dimensionless).
 * ALGEBRAIC[67] is i in component Ito (dimensionless).
 * ALGEBRAIC[34] is assp in component Ito (dimensionless).
 * STATES[22] is ap in component Ito (dimensionless).
 * ALGEBRAIC[50] is dti_develop in component Ito (dimensionless).
 * ALGEBRAIC[52] is dti_recover in component Ito (dimensionless).
 * ALGEBRAIC[54] is tiFp in component Ito (millisecond).
 * ALGEBRAIC[55] is tiSp in component Ito (millisecond).
 * STATES[23] is iFp in component Ito (dimensionless).
 * STATES[24] is iSp in component Ito (dimensionless).
 * ALGEBRAIC[68] is ip in component Ito (dimensionless).
 * CONSTANTS[97] is Gto in component Ito (milliS_per_microF).
 * ALGEBRAIC[69] is fItop in component Ito (dimensionless).
 * CONSTANTS[46] is Kmn in component ICaL (millimolar).
 * CONSTANTS[47] is k2n in component ICaL (per_millisecond).
 * CONSTANTS[48] is PCa_b in component ICaL (dimensionless).
 * ALGEBRAIC[6] is dss in component ICaL (dimensionless).
 * STATES[25] is d in component ICaL (dimensionless).
 * ALGEBRAIC[7] is fss in component ICaL (dimensionless).
 * CONSTANTS[98] is Aff in component ICaL (dimensionless).
 * CONSTANTS[109] is Afs in component ICaL (dimensionless).
 * STATES[26] is ff in component ICaL (dimensionless).
 * STATES[27] is fs in component ICaL (dimensionless).
 * ALGEBRAIC[71] is f in component ICaL (dimensionless).
 * ALGEBRAIC[19] is fcass in component ICaL (dimensionless).
 * ALGEBRAIC[72] is Afcaf in component ICaL (dimensionless).
 * ALGEBRAIC[73] is Afcas in component ICaL (dimensionless).
 * STATES[28] is fcaf in component ICaL (dimensionless).
 * STATES[29] is fcas in component ICaL (dimensionless).
 * ALGEBRAIC[74] is fca in component ICaL (dimensionless).
 * STATES[30] is jca in component ICaL (dimensionless).
 * STATES[31] is ffp in component ICaL (dimensionless).
 * ALGEBRAIC[75] is fp in component ICaL (dimensionless).
 * STATES[32] is fcafp in component ICaL (dimensionless).
 * ALGEBRAIC[76] is fcap in component ICaL (dimensionless).
 * ALGEBRAIC[8] is km2n in component ICaL (per_millisecond).
 * ALGEBRAIC[20] is anca in component ICaL (dimensionless).
 * STATES[33] is nca in component ICaL (dimensionless).
 * ALGEBRAIC[77] is PhiCaL in component ICaL (dimensionless).
 * ALGEBRAIC[78] is PhiCaNa in component ICaL (dimensionless).
 * ALGEBRAIC[79] is PhiCaK in component ICaL (dimensionless).
 * CONSTANTS[99] is PCa in component ICaL (dimensionless).
 * CONSTANTS[110] is PCap in component ICaL (dimensionless).
 * CONSTANTS[111] is PCaNa in component ICaL (dimensionless).
 * CONSTANTS[112] is PCaK in component ICaL (dimensionless).
 * CONSTANTS[115] is PCaNap in component ICaL (dimensionless).
 * CONSTANTS[116] is PCaKp in component ICaL (dimensionless).
 * ALGEBRAIC[80] is fICaLp in component ICaL (dimensionless).
 * ALGEBRAIC[21] is td in component ICaL (millisecond).
 * ALGEBRAIC[22] is tff in component ICaL (millisecond).
 * ALGEBRAIC[23] is tfs in component ICaL (millisecond).
 * ALGEBRAIC[35] is tfcaf in component ICaL (millisecond).
 * ALGEBRAIC[36] is tfcas in component ICaL (millisecond).
 * CONSTANTS[100] is tjca in component ICaL (millisecond).
 * ALGEBRAIC[37] is tffp in component ICaL (millisecond).
 * ALGEBRAIC[44] is tfcafp in component ICaL (millisecond).
 * CONSTANTS[49] is GKr_b in component IKr (milliS_per_microF).
 * CONSTANTS[101] is GKr in component IKr (milliS_per_microF).
 * ALGEBRAIC[9] is xrss in component IKr (dimensionless).
 * ALGEBRAIC[24] is txrf in component IKr (millisecond).
 * ALGEBRAIC[25] is txrs in component IKr (millisecond).
 * ALGEBRAIC[88] is Axrf in component IKr (dimensionless).
 * ALGEBRAIC[91] is Axrs in component IKr (dimensionless).
 * STATES[34] is xrf in component IKr (dimensionless).
 * STATES[35] is xrs in component IKr (dimensionless).
 * ALGEBRAIC[94] is xr in component IKr (dimensionless).
 * ALGEBRAIC[95] is rkr in component IKr (dimensionless).
 * CONSTANTS[50] is GKs_b in component IKs (milliS_per_microF).
 * CONSTANTS[102] is GKs in component IKs (milliS_per_microF).
 * ALGEBRAIC[10] is xs1ss in component IKs (dimensionless).
 * ALGEBRAIC[26] is xs2ss in component IKs (dimensionless).
 * ALGEBRAIC[27] is txs1 in component IKs (millisecond).
 * STATES[36] is xs1 in component IKs (dimensionless).
 * STATES[37] is xs2 in component IKs (dimensionless).
 * ALGEBRAIC[97] is KsCa in component IKs (dimensionless).
 * ALGEBRAIC[38] is txs2 in component IKs (millisecond).
 * CONSTANTS[103] is GK1 in component IK1 (milliS_per_microF).
 * CONSTANTS[51] is GK1_b in component IK1 (milliS_per_microF).
 * ALGEBRAIC[11] is xk1ss in component IK1 (dimensionless).
 * ALGEBRAIC[28] is txk1 in component IK1 (millisecond).
 * STATES[38] is xk1 in component IK1 (dimensionless).
 * ALGEBRAIC[99] is rk1 in component IK1 (millisecond).
 * CONSTANTS[52] is kna1 in component INaCa_i (per_millisecond).
 * CONSTANTS[53] is kna2 in component INaCa_i (per_millisecond).
 * CONSTANTS[54] is kna3 in component INaCa_i (per_millisecond).
 * CONSTANTS[55] is kasymm in component INaCa_i (dimensionless).
 * CONSTANTS[56] is wna in component INaCa_i (dimensionless).
 * CONSTANTS[57] is wca in component INaCa_i (dimensionless).
 * CONSTANTS[58] is wnaca in component INaCa_i (dimensionless).
 * CONSTANTS[59] is kcaon in component INaCa_i (per_millisecond).
 * CONSTANTS[60] is kcaoff in component INaCa_i (per_millisecond).
 * CONSTANTS[61] is qna in component INaCa_i (dimensionless).
 * CONSTANTS[62] is qca in component INaCa_i (dimensionless).
 * ALGEBRAIC[102] is hna in component INaCa_i (dimensionless).
 * ALGEBRAIC[101] is hca in component INaCa_i (dimensionless).
 * CONSTANTS[63] is KmCaAct in component INaCa_i (millimolar).
 * CONSTANTS[64] is Gncx_b in component INaCa_i (milliS_per_microF).
 * CONSTANTS[128] is Gncx in component INaCa_i (milliS_per_microF).
 * ALGEBRAIC[103] is h1_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[104] is h2_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[105] is h3_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[106] is h4_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[107] is h5_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[108] is h6_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[109] is h7_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[110] is h8_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[111] is h9_i in component INaCa_i (dimensionless).
 * CONSTANTS[122] is h10_i in component INaCa_i (dimensionless).
 * CONSTANTS[123] is h11_i in component INaCa_i (dimensionless).
 * CONSTANTS[124] is h12_i in component INaCa_i (dimensionless).
 * CONSTANTS[125] is k1_i in component INaCa_i (dimensionless).
 * CONSTANTS[126] is k2_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[112] is k3p_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[113] is k3pp_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[114] is k3_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[117] is k4_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[115] is k4p_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[116] is k4pp_i in component INaCa_i (dimensionless).
 * CONSTANTS[127] is k5_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[118] is k6_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[119] is k7_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[120] is k8_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[121] is x1_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[122] is x2_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[123] is x3_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[124] is x4_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[125] is E1_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[126] is E2_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[127] is E3_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[128] is E4_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[129] is allo_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[130] is JncxNa_i in component INaCa_i (millimolar_per_millisecond).
 * ALGEBRAIC[131] is JncxCa_i in component INaCa_i (millimolar_per_millisecond).
 * ALGEBRAIC[133] is h1_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[134] is h2_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[135] is h3_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[136] is h4_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[137] is h5_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[138] is h6_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[139] is h7_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[140] is h8_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[141] is h9_ss in component INaCa_i (dimensionless).
 * CONSTANTS[129] is h10_ss in component INaCa_i (dimensionless).
 * CONSTANTS[130] is h11_ss in component INaCa_i (dimensionless).
 * CONSTANTS[131] is h12_ss in component INaCa_i (dimensionless).
 * CONSTANTS[132] is k1_ss in component INaCa_i (dimensionless).
 * CONSTANTS[133] is k2_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[142] is k3p_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[143] is k3pp_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[144] is k3_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[147] is k4_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[145] is k4p_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[146] is k4pp_ss in component INaCa_i (dimensionless).
 * CONSTANTS[134] is k5_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[148] is k6_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[149] is k7_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[150] is k8_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[151] is x1_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[152] is x2_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[153] is x3_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[154] is x4_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[155] is E1_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[156] is E2_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[157] is E3_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[158] is E4_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[159] is allo_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[160] is JncxNa_ss in component INaCa_i (millimolar_per_millisecond).
 * ALGEBRAIC[161] is JncxCa_ss in component INaCa_i (millimolar_per_millisecond).
 * CONSTANTS[65] is k1p in component INaK (per_millisecond).
 * CONSTANTS[66] is k1m in component INaK (per_millisecond).
 * CONSTANTS[67] is k2p in component INaK (per_millisecond).
 * CONSTANTS[68] is k2m in component INaK (per_millisecond).
 * CONSTANTS[69] is k3p in component INaK (per_millisecond).
 * CONSTANTS[70] is k3m in component INaK (per_millisecond).
 * CONSTANTS[71] is k4p in component INaK (per_millisecond).
 * CONSTANTS[72] is k4m in component INaK (per_millisecond).
 * CONSTANTS[73] is Knai0 in component INaK (millimolar).
 * CONSTANTS[74] is Knao0 in component INaK (millimolar).
 * CONSTANTS[75] is delta in component INaK (millivolt).
 * CONSTANTS[76] is Kki in component INaK (per_millisecond).
 * CONSTANTS[77] is Kko in component INaK (per_millisecond).
 * CONSTANTS[78] is MgADP in component INaK (millimolar).
 * CONSTANTS[79] is MgATP in component INaK (millimolar).
 * CONSTANTS[80] is Kmgatp in component INaK (millimolar).
 * CONSTANTS[81] is H in component INaK (millimolar).
 * CONSTANTS[82] is eP in component INaK (dimensionless).
 * CONSTANTS[83] is Khp in component INaK (millimolar).
 * CONSTANTS[84] is Knap in component INaK (millimolar).
 * CONSTANTS[85] is Kxkur in component INaK (millimolar).
 * CONSTANTS[86] is Pnak_b in component INaK (milliS_per_microF).
 * CONSTANTS[138] is Pnak in component INaK (milliS_per_microF).
 * ALGEBRAIC[163] is Knai in component INaK (millimolar).
 * ALGEBRAIC[164] is Knao in component INaK (millimolar).
 * ALGEBRAIC[165] is P in component INaK (dimensionless).
 * ALGEBRAIC[166] is a1 in component INaK (dimensionless).
 * CONSTANTS[135] is b1 in component INaK (dimensionless).
 * CONSTANTS[136] is a2 in component INaK (dimensionless).
 * ALGEBRAIC[167] is b2 in component INaK (dimensionless).
 * ALGEBRAIC[168] is a3 in component INaK (dimensionless).
 * ALGEBRAIC[169] is b3 in component INaK (dimensionless).
 * CONSTANTS[137] is a4 in component INaK (dimensionless).
 * ALGEBRAIC[170] is b4 in component INaK (dimensionless).
 * ALGEBRAIC[171] is x1 in component INaK (dimensionless).
 * ALGEBRAIC[172] is x2 in component INaK (dimensionless).
 * ALGEBRAIC[173] is x3 in component INaK (dimensionless).
 * ALGEBRAIC[174] is x4 in component INaK (dimensionless).
 * ALGEBRAIC[175] is E1 in component INaK (dimensionless).
 * ALGEBRAIC[176] is E2 in component INaK (dimensionless).
 * ALGEBRAIC[177] is E3 in component INaK (dimensionless).
 * ALGEBRAIC[178] is E4 in component INaK (dimensionless).
 * ALGEBRAIC[179] is JnakNa in component INaK (millimolar_per_millisecond).
 * ALGEBRAIC[180] is JnakK in component INaK (millimolar_per_millisecond).
 * ALGEBRAIC[182] is xkb in component IKb (dimensionless).
 * CONSTANTS[87] is GKb_b in component IKb (milliS_per_microF).
 * CONSTANTS[105] is GKb in component IKb (milliS_per_microF).
 * CONSTANTS[88] is PNab in component INab (milliS_per_microF).
 * CONSTANTS[89] is PCab in component ICab (milliS_per_microF).
 * CONSTANTS[90] is GpCa in component IpCa (milliS_per_microF).
 * CONSTANTS[91] is KmCap in component IpCa (millimolar).
 * CONSTANTS[92] is bt in component ryr (millisecond).
 * CONSTANTS[106] is a_rel in component ryr (millisecond).
 * ALGEBRAIC[86] is Jrel_inf in component ryr (dimensionless).
 * ALGEBRAIC[92] is tau_rel in component ryr (millisecond).
 * ALGEBRAIC[87] is Jrel_infp in component ryr (dimensionless).
 * ALGEBRAIC[84] is Jrel_temp in component ryr (dimensionless).
 * ALGEBRAIC[93] is tau_relp in component ryr (millisecond).
 * STATES[39] is Jrelnp in component ryr (dimensionless).
 * STATES[40] is Jrelp in component ryr (dimensionless).
 * CONSTANTS[107] is btp in component ryr (millisecond).
 * CONSTANTS[114] is a_relp in component ryr (millisecond).
 * ALGEBRAIC[83] is Jrel_inf_temp in component ryr (dimensionless).
 * ALGEBRAIC[190] is fJrelp in component ryr (dimensionless).
 * ALGEBRAIC[89] is tau_rel_temp in component ryr (millisecond).
 * ALGEBRAIC[90] is tau_relp_temp in component ryr (millisecond).
 * CONSTANTS[108] is upScale in component SERCA (dimensionless).
 * ALGEBRAIC[192] is Jupnp in component SERCA (millimolar_per_millisecond).
 * ALGEBRAIC[193] is Jupp in component SERCA (millimolar_per_millisecond).
 * ALGEBRAIC[194] is fJupp in component SERCA (dimensionless).
 * ALGEBRAIC[195] is Jleak in component SERCA (millimolar_per_millisecond).
 * RATES[0] is d/dt v in component membrane (millivolt).
 * RATES[1] is d/dt CaMKt in component CaMK (millimolar).
 * RATES[3] is d/dt nai in component intracellular_ions (millimolar).
 * RATES[4] is d/dt nass in component intracellular_ions (millimolar).
 * RATES[5] is d/dt ki in component intracellular_ions (millimolar).
 * RATES[6] is d/dt kss in component intracellular_ions (millimolar).
 * RATES[9] is d/dt cai in component intracellular_ions (millimolar).
 * RATES[2] is d/dt cass in component intracellular_ions (millimolar).
 * RATES[7] is d/dt cansr in component intracellular_ions (millimolar).
 * RATES[8] is d/dt cajsr in component intracellular_ions (millimolar).
 * RATES[10] is d/dt m in component INa (dimensionless).
 * RATES[11] is d/dt hf in component INa (dimensionless).
 * RATES[12] is d/dt hs in component INa (dimensionless).
 * RATES[13] is d/dt j in component INa (dimensionless).
 * RATES[14] is d/dt hsp in component INa (dimensionless).
 * RATES[15] is d/dt jp in component INa (dimensionless).
 * RATES[16] is d/dt mL in component INaL (dimensionless).
 * RATES[17] is d/dt hL in component INaL (dimensionless).
 * RATES[18] is d/dt hLp in component INaL (dimensionless).
 * RATES[19] is d/dt a in component Ito (dimensionless).
 * RATES[20] is d/dt iF in component Ito (dimensionless).
 * RATES[21] is d/dt iS in component Ito (dimensionless).
 * RATES[22] is d/dt ap in component Ito (dimensionless).
 * RATES[23] is d/dt iFp in component Ito (dimensionless).
 * RATES[24] is d/dt iSp in component Ito (dimensionless).
 * RATES[25] is d/dt d in component ICaL (dimensionless).
 * RATES[26] is d/dt ff in component ICaL (dimensionless).
 * RATES[27] is d/dt fs in component ICaL (dimensionless).
 * RATES[28] is d/dt fcaf in component ICaL (dimensionless).
 * RATES[29] is d/dt fcas in component ICaL (dimensionless).
 * RATES[30] is d/dt jca in component ICaL (dimensionless).
 * RATES[31] is d/dt ffp in component ICaL (dimensionless).
 * RATES[32] is d/dt fcafp in component ICaL (dimensionless).
 * RATES[33] is d/dt nca in component ICaL (dimensionless).
 * RATES[34] is d/dt xrf in component IKr (dimensionless).
 * RATES[35] is d/dt xrs in component IKr (dimensionless).
 * RATES[36] is d/dt xs1 in component IKs (dimensionless).
 * RATES[37] is d/dt xs2 in component IKs (dimensionless).
 * RATES[38] is d/dt xk1 in component IK1 (dimensionless).
 * RATES[39] is d/dt Jrelnp in component ryr (dimensionless).
 * RATES[40] is d/dt Jrelp in component ryr (dimensionless).
 */

Ohara_Rudy_2011::Ohara_Rudy_2011()
{
algebraic_size = 198;
constants_size = 139+3;
states_size = 41;
ALGEBRAIC = new double[algebraic_size];
CONSTANTS = new double[constants_size];
RATES = new double[states_size];
STATES = new double[states_size];
}

Ohara_Rudy_2011::~Ohara_Rudy_2011()
{
delete []ALGEBRAIC;
delete []CONSTANTS;
delete []RATES;
delete []STATES;
}

void Ohara_Rudy_2011::___applyDutta()
{
CONSTANTS[102] *= 1.870;
CONSTANTS[101] *= 1.013;
CONSTANTS[103] *= 1.698;
CONSTANTS[99] *= 1.007;
CONSTANTS[96] *= 2.661;
}

void Ohara_Rudy_2011::___initConsts()
{
CONSTANTS[1] = 140;
CONSTANTS[2] = 1.8;
CONSTANTS[3] = 5.4;
CONSTANTS[4] = 8314;
CONSTANTS[5] = 310;
CONSTANTS[6] = 96485;
CONSTANTS[7] = 1;
CONSTANTS[8] = 2;
CONSTANTS[9] = 1;
CONSTANTS[10] = 0.01;
CONSTANTS[11] = 0.0011;
CONSTANTS[stim_start] = 10.0;
CONSTANTS[BCL] = 1000.0;
CONSTANTS[cumm_BCL] = 0.0;
STATES[0] = -87;
CONSTANTS[12] = -80;
CONSTANTS[13] = 0.5;
CONSTANTS[14] = 0.15;
CONSTANTS[15] = 0.05;
CONSTANTS[16] = 0.00068;
CONSTANTS[17] = 0.05;
CONSTANTS[18] = 0.0015;
STATES[1] = 0;
STATES[2] = 1e-4;
CONSTANTS[19] = 0.05;
CONSTANTS[20] = 0.00238;
CONSTANTS[21] = 0.07;
CONSTANTS[22] = 0.0005;
CONSTANTS[23] = 0.047;
CONSTANTS[24] = 0.00087;
CONSTANTS[25] = 1.124;
CONSTANTS[26] = 0.0087;
CONSTANTS[27] = 10;
CONSTANTS[28] = 0.8;
STATES[3] = 7;
STATES[4] = 7;
STATES[5] = 145;
STATES[6] = 145;
STATES[7] = 1.2;
STATES[8] = 1.2;
STATES[9] = 1e-4;
CONSTANTS[29] = 1;
CONSTANTS[30] = 0.01833;
CONSTANTS[31] = 39.57;
CONSTANTS[32] = 9.871;
CONSTANTS[33] = 11.64;
CONSTANTS[34] = 34.77;
CONSTANTS[35] = 6.765;
CONSTANTS[36] = 8.552;
CONSTANTS[37] = 77.42;
CONSTANTS[38] = 5.955;
STATES[10] = 0;
CONSTANTS[39] = 82.9;
CONSTANTS[40] = 6.086;
CONSTANTS[41] = 0.99;
STATES[11] = 1;
STATES[12] = 1;
CONSTANTS[42] = 75;
STATES[13] = 1;
STATES[14] = 1;
STATES[15] = 1;
STATES[16] = 0;
CONSTANTS[43] = 200;
STATES[17] = 1;
STATES[18] = 1;
CONSTANTS[44] = 0.0075;
CONSTANTS[45] = 0.02;
STATES[19] = 0;
STATES[20] = 1;
STATES[21] = 1;
STATES[22] = 0;
STATES[23] = 1;
STATES[24] = 1;
CONSTANTS[46] = 0.002;
CONSTANTS[47] = 1000;
CONSTANTS[48] = 0.0001;
STATES[25] = 0;
STATES[26] = 1;
STATES[27] = 1;
STATES[28] = 1;
STATES[29] = 1;
STATES[30] = 1;
STATES[31] = 1;
STATES[32] = 1;
STATES[33] = 0;
CONSTANTS[49] = 0.046;
STATES[34] = 0;
STATES[35] = 0;
CONSTANTS[50] = 0.0034;
STATES[36] = 0;
STATES[37] = 0;
CONSTANTS[51] = 0.1908;
STATES[38] = 1;
CONSTANTS[52] = 15;
CONSTANTS[53] = 5;
CONSTANTS[54] = 88.12;
CONSTANTS[55] = 12.5;
CONSTANTS[56] = 6e4;
CONSTANTS[57] = 6e4;
CONSTANTS[58] = 5e3;
CONSTANTS[59] = 1.5e6;
CONSTANTS[60] = 5e3;
CONSTANTS[61] = 0.5224;
CONSTANTS[62] = 0.167;
CONSTANTS[63] = 150e-6;
CONSTANTS[64] = 0.0008;
CONSTANTS[65] = 949.5;
CONSTANTS[66] = 182.4;
CONSTANTS[67] = 687.2;
CONSTANTS[68] = 39.4;
CONSTANTS[69] = 1899;
CONSTANTS[70] = 79300;
CONSTANTS[71] = 639;
CONSTANTS[72] = 40;
CONSTANTS[73] = 9.073;
CONSTANTS[74] = 27.78;
CONSTANTS[75] = -0.155;
CONSTANTS[76] = 0.5;
CONSTANTS[77] = 0.3582;
CONSTANTS[78] = 0.05;
CONSTANTS[79] = 9.8;
CONSTANTS[80] = 1.698e-7;
CONSTANTS[81] = 1e-7;
CONSTANTS[82] = 4.2;
CONSTANTS[83] = 1.698e-7;
CONSTANTS[84] = 224;
CONSTANTS[85] = 292;
CONSTANTS[86] = 30;
CONSTANTS[87] = 0.003;
CONSTANTS[88] = 3.75e-10;
CONSTANTS[89] = 2.5e-8;
CONSTANTS[90] = 0.0005;
CONSTANTS[91] = 0.0005;
CONSTANTS[92] = 4.75;
STATES[39] = 0;
STATES[40] = 0;
CONSTANTS[93] = (CONSTANTS[0]==1.00000 ?  CONSTANTS[19]*1.30000 : CONSTANTS[19]);
CONSTANTS[94] = 1.00000 - CONSTANTS[41];
CONSTANTS[95] =  3.00000*CONSTANTS[43];
CONSTANTS[96] = (CONSTANTS[0]==1.00000 ?  CONSTANTS[44]*0.600000 : CONSTANTS[44]);
CONSTANTS[97] = (CONSTANTS[0]==1.00000 ?  CONSTANTS[45]*4.00000 : CONSTANTS[0]==2.00000 ?  CONSTANTS[45]*4.00000 : CONSTANTS[45]);
CONSTANTS[98] = 0.600000;
CONSTANTS[99] = (CONSTANTS[0]==1.00000 ?  CONSTANTS[48]*1.20000 : CONSTANTS[0]==2.00000 ?  CONSTANTS[48]*2.50000 : CONSTANTS[48]);
CONSTANTS[100] = 75.0000;
CONSTANTS[101] = (CONSTANTS[0]==1.00000 ?  CONSTANTS[49]*1.30000 : CONSTANTS[0]==2.00000 ?  CONSTANTS[49]*0.800000 : CONSTANTS[49]);
CONSTANTS[102] = (CONSTANTS[0]==1.00000 ?  CONSTANTS[50]*1.40000 : CONSTANTS[50]);
CONSTANTS[103] = (CONSTANTS[0]==1.00000 ?  CONSTANTS[51]*1.20000 : CONSTANTS[0]==2.00000 ?  CONSTANTS[51]*1.30000 : CONSTANTS[51]);
CONSTANTS[104] =  1000.00*3.14000*CONSTANTS[11]*CONSTANTS[11]*CONSTANTS[10];
CONSTANTS[105] = (CONSTANTS[0]==1.00000 ?  CONSTANTS[87]*0.600000 : CONSTANTS[87]);
CONSTANTS[106] =  0.500000*CONSTANTS[92];
CONSTANTS[107] =  1.25000*CONSTANTS[92];
CONSTANTS[108] = (CONSTANTS[0]==1.00000 ? 1.30000 : 1.00000);
CONSTANTS[109] = 1.00000 - CONSTANTS[98];
CONSTANTS[110] =  1.10000*CONSTANTS[99];
CONSTANTS[111] =  0.00125000*CONSTANTS[99];
CONSTANTS[112] =  0.000357400*CONSTANTS[99];
CONSTANTS[113] =  2.00000*3.14000*CONSTANTS[11]*CONSTANTS[11]+ 2.00000*3.14000*CONSTANTS[11]*CONSTANTS[10];
CONSTANTS[114] =  0.500000*CONSTANTS[107];
CONSTANTS[115] =  0.00125000*CONSTANTS[110];
CONSTANTS[116] =  0.000357400*CONSTANTS[110];
CONSTANTS[117] =  2.00000*CONSTANTS[113];
CONSTANTS[118] =  0.680000*CONSTANTS[104];
CONSTANTS[119] =  0.0552000*CONSTANTS[104];
CONSTANTS[120] =  0.00480000*CONSTANTS[104];
CONSTANTS[121] =  0.0200000*CONSTANTS[104];
CONSTANTS[122] = CONSTANTS[55]+1.00000+ (CONSTANTS[1]/CONSTANTS[52])*(1.00000+CONSTANTS[1]/CONSTANTS[53]);
CONSTANTS[123] = ( CONSTANTS[1]*CONSTANTS[1])/( CONSTANTS[122]*CONSTANTS[52]*CONSTANTS[53]);
CONSTANTS[124] = 1.00000/CONSTANTS[122];
CONSTANTS[125] =  CONSTANTS[124]*CONSTANTS[2]*CONSTANTS[59];
CONSTANTS[126] = CONSTANTS[60];
CONSTANTS[127] = CONSTANTS[60];
CONSTANTS[128] = (CONSTANTS[0]==1.00000 ?  CONSTANTS[64]*1.10000 : CONSTANTS[0]==2.00000 ?  CONSTANTS[64]*1.40000 : CONSTANTS[64]);
CONSTANTS[129] = CONSTANTS[55]+1.00000+ (CONSTANTS[1]/CONSTANTS[52])*(1.00000+CONSTANTS[1]/CONSTANTS[53]);
CONSTANTS[130] = ( CONSTANTS[1]*CONSTANTS[1])/( CONSTANTS[129]*CONSTANTS[52]*CONSTANTS[53]);
CONSTANTS[131] = 1.00000/CONSTANTS[129];
CONSTANTS[132] =  CONSTANTS[131]*CONSTANTS[2]*CONSTANTS[59];
CONSTANTS[133] = CONSTANTS[60];
CONSTANTS[134] = CONSTANTS[60];
CONSTANTS[135] =  CONSTANTS[66]*CONSTANTS[78];
CONSTANTS[136] = CONSTANTS[67];
CONSTANTS[137] = (( CONSTANTS[71]*CONSTANTS[79])/CONSTANTS[80])/(1.00000+CONSTANTS[79]/CONSTANTS[80]);
CONSTANTS[138] = (CONSTANTS[0]==1.00000 ?  CONSTANTS[86]*0.900000 : CONSTANTS[0]==2.00000 ?  CONSTANTS[86]*0.700000 : CONSTANTS[86]);
}

void Ohara_Rudy_2011::initConsts()
{
CONSTANTS[0] = 0;
___initConsts();
}

void Ohara_Rudy_2011::initConsts(double type)
{
CONSTANTS[0] = type;
___initConsts();
}

void Ohara_Rudy_2011::initConsts(double type, bool is_dutta)
{
CONSTANTS[0] = type;
___initConsts();
if(is_dutta == true) ___applyDutta();
}


void Ohara_Rudy_2011::initConsts(double type, double D, double *hill, bool is_dutta)
{
initConsts(type, is_dutta);

CONSTANTS[PCa]  *= (hill[0] > 10E-14 && hill[1] > 10E-14) ?
                  pow(1.+pow(D/hill[0], hill[1]), -1) : 1.;
CONSTANTS[GK1] *= (hill[2] > 10E-14 && hill[3] > 10E-14) ?
                  pow(1.+pow(D/hill[2], hill[3]), -1) : 1.;
CONSTANTS[GKs] *= (hill[4] > 10E-14 && hill[5] > 10E-14) ?
                  pow(1.+pow(D/hill[4], hill[5]), -1) : 1.;
CONSTANTS[GNa]  *= (hill[6] > 10E-14 && hill[7] > 10E-14) ?
                  pow(1.+pow(D/hill[6], hill[7]), -1) : 1.;
CONSTANTS[GNaL]  *= (hill[8] > 10E-14 && hill[9] > 10E-14) ?
                  pow(1.+pow(D/hill[8], hill[9]), -1) : 1.;
CONSTANTS[Gto]  *= (hill[10] > 10E-14 && hill[11] > 10E-14) ?
                  pow(1.+pow(D/hill[10], hill[11]), -1) : 1.;
CONSTANTS[GKr] *= (hill[12] > 10E-14 && hill[13] > 10E-14) ?
                  pow(1.+pow(D/hill[12], hill[13]), -1) : 1.;
}


void Ohara_Rudy_2011::computeRates( double TIME, double *CONSTANTS, double *RATES, double *STATES, double *ALGEBRAIC )
{
if( glob_var::is_hrv == true ){
ALGEBRAIC[Istim] = ( CONSTANTS[stim_start] + CONSTANTS[cumm_BCL] <= TIME && TIME < CONSTANTS[stim_start]+CONSTANTS[cumm_BCL]+CONSTANTS[duration] ? CONSTANTS[amp] : 0.000000);
}
else{
#if defined(SINGLE_CELL) && !defined(TISSUE)
ALGEBRAIC[Istim] = (TIME>=CONSTANTS[stim_start] && (TIME - CONSTANTS[stim_start]) - floor((TIME - CONSTANTS[stim_start])/CONSTANTS[BCL])*CONSTANTS[BCL]<=CONSTANTS[duration] ? CONSTANTS[amp] : 0.000000);
#elif defined(TISSUE)
if(isS1) ALGEBRAIC[Istim] = CONSTANTS[amp];
else ALGEBRAIC[Istim] = 0.0;
#endif
}
ALGEBRAIC[2] = 1.00000/(1.00000+exp((STATES[0]+87.6100)/7.48800));
RATES[17] = (ALGEBRAIC[2] - STATES[17])/CONSTANTS[43];
ALGEBRAIC[3] = 1.00000/(1.00000+exp((STATES[0]+93.8100)/7.48800));
RATES[18] = (ALGEBRAIC[3] - STATES[18])/CONSTANTS[95];
ALGEBRAIC[0] = 1.00000/(1.00000+exp(- (STATES[0]+CONSTANTS[31])/CONSTANTS[32]));
ALGEBRAIC[13] = 1.00000/( CONSTANTS[35]*exp((STATES[0]+CONSTANTS[33])/CONSTANTS[34])+ CONSTANTS[36]*exp(- (STATES[0]+CONSTANTS[37])/CONSTANTS[38]));
RATES[10] = (ALGEBRAIC[0] - STATES[10])/ALGEBRAIC[13];
ALGEBRAIC[1] = 1.00000/(1.00000+exp((STATES[0]+CONSTANTS[39])/CONSTANTS[40]));
ALGEBRAIC[14] = 1.00000/( 1.43200e-05*exp(- (STATES[0]+1.19600)/6.28500)+ 6.14900*exp((STATES[0]+0.509600)/20.2700));
RATES[11] = (ALGEBRAIC[1] - STATES[11])/ALGEBRAIC[14];
ALGEBRAIC[15] = 1.00000/( 0.00979400*exp(- (STATES[0]+17.9500)/28.0500)+ 0.334300*exp((STATES[0]+5.73000)/56.6600));
RATES[12] = (ALGEBRAIC[1] - STATES[12])/ALGEBRAIC[15];
ALGEBRAIC[4] = 1.00000/(1.00000+exp(- (STATES[0] - 14.3400)/14.8200));
ALGEBRAIC[17] = 1.05150/(1.00000/( 1.20890*(1.00000+exp(- (STATES[0] - 18.4099)/29.3814)))+3.50000/(1.00000+exp((STATES[0]+100.000)/29.3814)));
RATES[19] = (ALGEBRAIC[4] - STATES[19])/ALGEBRAIC[17];
ALGEBRAIC[6] = 1.00000/(1.00000+exp(- (STATES[0]+3.94000)/4.23000));
ALGEBRAIC[21] = 0.600000+1.00000/(exp( - 0.0500000*(STATES[0]+6.00000))+exp( 0.0900000*(STATES[0]+14.0000)));
RATES[25] = (ALGEBRAIC[6] - STATES[25])/ALGEBRAIC[21];
ALGEBRAIC[7] = 1.00000/(1.00000+exp((STATES[0]+19.5800)/3.69600));
ALGEBRAIC[22] = 7.00000+1.00000/( 0.00450000*exp(- (STATES[0]+20.0000)/10.0000)+ 0.00450000*exp((STATES[0]+20.0000)/10.0000));
RATES[26] = (ALGEBRAIC[7] - STATES[26])/ALGEBRAIC[22];
ALGEBRAIC[23] = 1000.00+1.00000/( 3.50000e-05*exp(- (STATES[0]+5.00000)/4.00000)+ 3.50000e-05*exp((STATES[0]+5.00000)/6.00000));
RATES[27] = (ALGEBRAIC[7] - STATES[27])/ALGEBRAIC[23];
ALGEBRAIC[19] = ALGEBRAIC[7];
RATES[30] = (ALGEBRAIC[19] - STATES[30])/CONSTANTS[100];
ALGEBRAIC[8] =  STATES[30]*1.00000;
ALGEBRAIC[20] = 1.00000/(CONSTANTS[47]/ALGEBRAIC[8]+pow(1.00000+CONSTANTS[46]/STATES[2], 4.00000));
RATES[33] =  ALGEBRAIC[20]*CONSTANTS[47] -  STATES[33]*ALGEBRAIC[8];
ALGEBRAIC[9] = 1.00000/(1.00000+exp(- (STATES[0]+8.33700)/6.78900));
ALGEBRAIC[24] = 12.9800+1.00000/( 0.365200*exp((STATES[0] - 31.6600)/3.86900)+ 4.12300e-05*exp(- (STATES[0] - 47.7800)/20.3800));
RATES[34] = (ALGEBRAIC[9] - STATES[34])/ALGEBRAIC[24];
ALGEBRAIC[25] = 1.86500+1.00000/( 0.0662900*exp((STATES[0] - 34.7000)/7.35500)+ 1.12800e-05*exp(- (STATES[0] - 29.7400)/25.9400));
RATES[35] = (ALGEBRAIC[9] - STATES[35])/ALGEBRAIC[25];
ALGEBRAIC[10] = 1.00000/(1.00000+exp(- (STATES[0]+11.6000)/8.93200));
ALGEBRAIC[27] = 817.300+1.00000/( 0.000232600*exp((STATES[0]+48.2800)/17.8000)+ 0.00129200*exp(- (STATES[0]+210.000)/230.000));
RATES[36] = (ALGEBRAIC[10] - STATES[36])/ALGEBRAIC[27];
ALGEBRAIC[11] = 1.00000/(1.00000+exp(- (STATES[0]+ 2.55380*CONSTANTS[3]+144.590)/( 1.56920*CONSTANTS[3]+3.81150)));
ALGEBRAIC[28] = 122.200/(exp(- (STATES[0]+127.200)/20.3600)+exp((STATES[0]+236.800)/69.3300));
RATES[38] = (ALGEBRAIC[11] - STATES[38])/ALGEBRAIC[28];
ALGEBRAIC[16] = ALGEBRAIC[1];
ALGEBRAIC[30] = 2.03800+1.00000/( 0.0213600*exp(- (STATES[0]+100.600)/8.28100)+ 0.305200*exp((STATES[0]+0.994100)/38.4500));
RATES[13] = (ALGEBRAIC[16] - STATES[13])/ALGEBRAIC[30];
ALGEBRAIC[34] = 1.00000/(1.00000+exp(- (STATES[0] - 24.3400)/14.8200));
RATES[22] = (ALGEBRAIC[34] - STATES[22])/ALGEBRAIC[17];
ALGEBRAIC[35] = 7.00000+1.00000/( 0.0400000*exp(- (STATES[0] - 4.00000)/7.00000)+ 0.0400000*exp((STATES[0] - 4.00000)/7.00000));
RATES[28] = (ALGEBRAIC[19] - STATES[28])/ALGEBRAIC[35];
ALGEBRAIC[36] = 100.000+1.00000/( 0.000120000*exp(- STATES[0]/3.00000)+ 0.000120000*exp(STATES[0]/7.00000));
RATES[29] = (ALGEBRAIC[19] - STATES[29])/ALGEBRAIC[36];
ALGEBRAIC[37] =  2.50000*ALGEBRAIC[22];
RATES[31] = (ALGEBRAIC[7] - STATES[31])/ALGEBRAIC[37];
ALGEBRAIC[26] = ALGEBRAIC[10];
ALGEBRAIC[38] = 1.00000/( 0.0100000*exp((STATES[0] - 50.0000)/20.0000)+ 0.0193000*exp(- (STATES[0]+66.5400)/31.0000));
RATES[37] = (ALGEBRAIC[26] - STATES[37])/ALGEBRAIC[38];
ALGEBRAIC[45] = ( CONSTANTS[17]*(1.00000 - STATES[1]))/(1.00000+CONSTANTS[18]/STATES[2]);
RATES[1] =  CONSTANTS[15]*ALGEBRAIC[45]*(ALGEBRAIC[45]+STATES[1]) -  CONSTANTS[16]*STATES[1];
ALGEBRAIC[31] = 1.00000/(1.00000+exp((STATES[0]+89.1000)/6.08600));
ALGEBRAIC[40] =  3.00000*ALGEBRAIC[15];
RATES[14] = (ALGEBRAIC[31] - STATES[14])/ALGEBRAIC[40];
ALGEBRAIC[41] =  1.46000*ALGEBRAIC[30];
RATES[15] = (ALGEBRAIC[16] - STATES[15])/ALGEBRAIC[41];
ALGEBRAIC[32] = 1.00000/(1.00000+exp(- (STATES[0]+42.8500)/5.26400));
ALGEBRAIC[42] = ALGEBRAIC[13];
RATES[16] = (ALGEBRAIC[32] - STATES[16])/ALGEBRAIC[42];
ALGEBRAIC[44] =  2.50000*ALGEBRAIC[35];
RATES[32] = (ALGEBRAIC[19] - STATES[32])/ALGEBRAIC[44];
ALGEBRAIC[5] = 1.00000/(1.00000+exp((STATES[0]+43.9400)/5.71100));
ALGEBRAIC[18] = (CONSTANTS[0]==1.00000 ? 1.00000 - 0.950000/(1.00000+exp((STATES[0]+70.0000)/5.00000)) : 1.00000);
ALGEBRAIC[33] = 4.56200+1.00000/( 0.393300*exp(- (STATES[0]+100.000)/100.000)+ 0.0800400*exp((STATES[0]+50.0000)/16.5900));
ALGEBRAIC[46] =  ALGEBRAIC[33]*ALGEBRAIC[18];
RATES[20] = (ALGEBRAIC[5] - STATES[20])/ALGEBRAIC[46];
ALGEBRAIC[43] = 23.6200+1.00000/( 0.00141600*exp(- (STATES[0]+96.5200)/59.0500)+ 1.78000e-08*exp((STATES[0]+114.100)/8.07900));
ALGEBRAIC[48] =  ALGEBRAIC[43]*ALGEBRAIC[18];
RATES[21] = (ALGEBRAIC[5] - STATES[21])/ALGEBRAIC[48];
ALGEBRAIC[50] = 1.35400+0.000100000/(exp((STATES[0] - 167.400)/15.8900)+exp(- (STATES[0] - 12.2300)/0.215400));
ALGEBRAIC[52] = 1.00000 - 0.500000/(1.00000+exp((STATES[0]+70.0000)/20.0000));
ALGEBRAIC[54] =  ALGEBRAIC[50]*ALGEBRAIC[52]*ALGEBRAIC[46];
RATES[23] = (ALGEBRAIC[5] - STATES[23])/ALGEBRAIC[54];
ALGEBRAIC[55] =  ALGEBRAIC[50]*ALGEBRAIC[52]*ALGEBRAIC[48];
RATES[24] = (ALGEBRAIC[5] - STATES[24])/ALGEBRAIC[55];
ALGEBRAIC[71] =  CONSTANTS[98]*STATES[26]+ CONSTANTS[109]*STATES[27];
ALGEBRAIC[72] = 0.300000+0.600000/(1.00000+exp((STATES[0] - 10.0000)/10.0000));
ALGEBRAIC[73] = 1.00000 - ALGEBRAIC[72];
ALGEBRAIC[74] =  ALGEBRAIC[72]*STATES[28]+ ALGEBRAIC[73]*STATES[29];
ALGEBRAIC[75] =  CONSTANTS[98]*STATES[31]+ CONSTANTS[109]*STATES[27];
ALGEBRAIC[76] =  ALGEBRAIC[72]*STATES[32]+ ALGEBRAIC[73]*STATES[29];
ALGEBRAIC[29] = ( STATES[0]*CONSTANTS[6]*CONSTANTS[6])/( CONSTANTS[4]*CONSTANTS[5]);
ALGEBRAIC[39] = ( STATES[0]*CONSTANTS[6])/( CONSTANTS[4]*CONSTANTS[5]);
ALGEBRAIC[77] = ( 4.00000*ALGEBRAIC[29]*( STATES[2]*exp( 2.00000*ALGEBRAIC[39]) -  0.341000*CONSTANTS[2]))/(exp( 2.00000*ALGEBRAIC[39]) - 1.00000);
ALGEBRAIC[47] = ALGEBRAIC[45]+STATES[1];
ALGEBRAIC[80] = 1.00000/(1.00000+CONSTANTS[14]/ALGEBRAIC[47]);
ALGEBRAIC[81] =  (1.00000 - ALGEBRAIC[80])*CONSTANTS[99]*ALGEBRAIC[77]*STATES[25]*( ALGEBRAIC[71]*(1.00000 - STATES[33])+ STATES[30]*ALGEBRAIC[74]*STATES[33])+ ALGEBRAIC[80]*CONSTANTS[110]*ALGEBRAIC[77]*STATES[25]*( ALGEBRAIC[75]*(1.00000 - STATES[33])+ STATES[30]*ALGEBRAIC[76]*STATES[33]);
ALGEBRAIC[83] = ( CONSTANTS[106]*- ALGEBRAIC[81])/(1.00000+ 1.00000*pow(1.50000/STATES[8], 8.00000));
ALGEBRAIC[86] = (CONSTANTS[0]==2.00000 ?  ALGEBRAIC[83]*1.70000 : ALGEBRAIC[83]);
ALGEBRAIC[89] = CONSTANTS[92]/(1.00000+0.0123000/STATES[8]);
ALGEBRAIC[92] = (ALGEBRAIC[89]<0.00100000 ? 0.00100000 : ALGEBRAIC[89]);
RATES[39] = (ALGEBRAIC[86] - STATES[39])/ALGEBRAIC[92];
ALGEBRAIC[84] = ( CONSTANTS[114]*- ALGEBRAIC[81])/(1.00000+pow(1.50000/STATES[8], 8.00000));
ALGEBRAIC[87] = (CONSTANTS[0]==2.00000 ?  ALGEBRAIC[84]*1.70000 : ALGEBRAIC[84]);
ALGEBRAIC[90] = CONSTANTS[107]/(1.00000+0.0123000/STATES[8]);
ALGEBRAIC[93] = (ALGEBRAIC[90]<0.00100000 ? 0.00100000 : ALGEBRAIC[90]);
RATES[40] = (ALGEBRAIC[87] - STATES[40])/ALGEBRAIC[93];
ALGEBRAIC[57] =  (( CONSTANTS[4]*CONSTANTS[5])/CONSTANTS[6])*log(CONSTANTS[3]/STATES[5]);
ALGEBRAIC[65] = 1.00000/(1.00000+exp((STATES[0] - 213.600)/151.200));
ALGEBRAIC[66] = 1.00000 - ALGEBRAIC[65];
ALGEBRAIC[67] =  ALGEBRAIC[65]*STATES[20]+ ALGEBRAIC[66]*STATES[21];
ALGEBRAIC[68] =  ALGEBRAIC[65]*STATES[23]+ ALGEBRAIC[66]*STATES[24];
ALGEBRAIC[69] = 1.00000/(1.00000+CONSTANTS[14]/ALGEBRAIC[47]);
ALGEBRAIC[70] =  CONSTANTS[97]*(STATES[0] - ALGEBRAIC[57])*( (1.00000 - ALGEBRAIC[69])*STATES[19]*ALGEBRAIC[67]+ ALGEBRAIC[69]*STATES[22]*ALGEBRAIC[68]);
ALGEBRAIC[88] = 1.00000/(1.00000+exp((STATES[0]+54.8100)/38.2100));
ALGEBRAIC[91] = 1.00000 - ALGEBRAIC[88];
ALGEBRAIC[94] =  ALGEBRAIC[88]*STATES[34]+ ALGEBRAIC[91]*STATES[35];
ALGEBRAIC[95] = ( (1.00000/(1.00000+exp((STATES[0]+55.0000)/75.0000)))*1.00000)/(1.00000+exp((STATES[0] - 10.0000)/30.0000));
ALGEBRAIC[96] =  CONSTANTS[101]* pow((CONSTANTS[3]/5.40000), 1.0 / 2)*ALGEBRAIC[94]*ALGEBRAIC[95]*(STATES[0] - ALGEBRAIC[57]);
ALGEBRAIC[58] =  (( CONSTANTS[4]*CONSTANTS[5])/CONSTANTS[6])*log((CONSTANTS[3]+ CONSTANTS[30]*CONSTANTS[1])/(STATES[5]+ CONSTANTS[30]*STATES[3]));
ALGEBRAIC[97] = 1.00000+0.600000/(1.00000+pow(3.80000e-05/STATES[9], 1.40000));
ALGEBRAIC[98] =  CONSTANTS[102]*ALGEBRAIC[97]*STATES[36]*STATES[37]*(STATES[0] - ALGEBRAIC[58]);
ALGEBRAIC[99] = 1.00000/(1.00000+exp(((STATES[0]+105.800) -  2.60000*CONSTANTS[3])/9.49300));
ALGEBRAIC[100] =  CONSTANTS[103]* pow(CONSTANTS[3], 1.0 / 2)*ALGEBRAIC[99]*STATES[38]*(STATES[0] - ALGEBRAIC[57]);
ALGEBRAIC[164] =  CONSTANTS[74]*exp(( (1.00000 - CONSTANTS[75])*STATES[0]*CONSTANTS[6])/( 3.00000*CONSTANTS[4]*CONSTANTS[5]));
ALGEBRAIC[168] = ( CONSTANTS[69]*pow(CONSTANTS[3]/CONSTANTS[77], 2.00000))/((pow(1.00000+CONSTANTS[1]/ALGEBRAIC[164], 3.00000)+pow(1.00000+CONSTANTS[3]/CONSTANTS[77], 2.00000)) - 1.00000);
ALGEBRAIC[165] = CONSTANTS[82]/(1.00000+CONSTANTS[81]/CONSTANTS[83]+STATES[3]/CONSTANTS[84]+STATES[5]/CONSTANTS[85]);
ALGEBRAIC[169] = ( CONSTANTS[70]*ALGEBRAIC[165]*CONSTANTS[81])/(1.00000+CONSTANTS[79]/CONSTANTS[80]);
ALGEBRAIC[163] =  CONSTANTS[73]*exp(( CONSTANTS[75]*STATES[0]*CONSTANTS[6])/( 3.00000*CONSTANTS[4]*CONSTANTS[5]));
ALGEBRAIC[166] = ( CONSTANTS[65]*pow(STATES[3]/ALGEBRAIC[163], 3.00000))/((pow(1.00000+STATES[3]/ALGEBRAIC[163], 3.00000)+pow(1.00000+STATES[5]/CONSTANTS[76], 2.00000)) - 1.00000);
ALGEBRAIC[167] = ( CONSTANTS[68]*pow(CONSTANTS[1]/ALGEBRAIC[164], 3.00000))/((pow(1.00000+CONSTANTS[1]/ALGEBRAIC[164], 3.00000)+pow(1.00000+CONSTANTS[3]/CONSTANTS[77], 2.00000)) - 1.00000);
ALGEBRAIC[170] = ( CONSTANTS[72]*pow(STATES[5]/CONSTANTS[76], 2.00000))/((pow(1.00000+STATES[3]/ALGEBRAIC[163], 3.00000)+pow(1.00000+STATES[5]/CONSTANTS[76], 2.00000)) - 1.00000);
ALGEBRAIC[171] =  CONSTANTS[137]*ALGEBRAIC[166]*CONSTANTS[136]+ ALGEBRAIC[167]*ALGEBRAIC[170]*ALGEBRAIC[169]+ CONSTANTS[136]*ALGEBRAIC[170]*ALGEBRAIC[169]+ ALGEBRAIC[169]*ALGEBRAIC[166]*CONSTANTS[136];
ALGEBRAIC[172] =  ALGEBRAIC[167]*CONSTANTS[135]*ALGEBRAIC[170]+ ALGEBRAIC[166]*CONSTANTS[136]*ALGEBRAIC[168]+ ALGEBRAIC[168]*CONSTANTS[135]*ALGEBRAIC[170]+ CONSTANTS[136]*ALGEBRAIC[168]*ALGEBRAIC[170];
ALGEBRAIC[173] =  CONSTANTS[136]*ALGEBRAIC[168]*CONSTANTS[137]+ ALGEBRAIC[169]*ALGEBRAIC[167]*CONSTANTS[135]+ ALGEBRAIC[167]*CONSTANTS[135]*CONSTANTS[137]+ ALGEBRAIC[168]*CONSTANTS[137]*CONSTANTS[135];
ALGEBRAIC[174] =  ALGEBRAIC[170]*ALGEBRAIC[169]*ALGEBRAIC[167]+ ALGEBRAIC[168]*CONSTANTS[137]*ALGEBRAIC[166]+ ALGEBRAIC[167]*CONSTANTS[137]*ALGEBRAIC[166]+ ALGEBRAIC[169]*ALGEBRAIC[167]*ALGEBRAIC[166];
ALGEBRAIC[175] = ALGEBRAIC[171]/(ALGEBRAIC[171]+ALGEBRAIC[172]+ALGEBRAIC[173]+ALGEBRAIC[174]);
ALGEBRAIC[176] = ALGEBRAIC[172]/(ALGEBRAIC[171]+ALGEBRAIC[172]+ALGEBRAIC[173]+ALGEBRAIC[174]);
ALGEBRAIC[179] =  3.00000*( ALGEBRAIC[175]*ALGEBRAIC[168] -  ALGEBRAIC[176]*ALGEBRAIC[169]);
ALGEBRAIC[177] = ALGEBRAIC[173]/(ALGEBRAIC[171]+ALGEBRAIC[172]+ALGEBRAIC[173]+ALGEBRAIC[174]);
ALGEBRAIC[178] = ALGEBRAIC[174]/(ALGEBRAIC[171]+ALGEBRAIC[172]+ALGEBRAIC[173]+ALGEBRAIC[174]);
ALGEBRAIC[180] =  2.00000*( ALGEBRAIC[178]*CONSTANTS[135] -  ALGEBRAIC[177]*ALGEBRAIC[166]);
ALGEBRAIC[181] =  CONSTANTS[138]*( CONSTANTS[7]*ALGEBRAIC[179]+ CONSTANTS[9]*ALGEBRAIC[180]);
ALGEBRAIC[182] = 1.00000/(1.00000+exp(- (STATES[0] - 14.4800)/18.3400));
ALGEBRAIC[183] =  CONSTANTS[105]*ALGEBRAIC[182]*(STATES[0] - ALGEBRAIC[57]);
ALGEBRAIC[185] = (STATES[6] - STATES[5])/2.00000;
RATES[5] = ( - ((ALGEBRAIC[70]+ALGEBRAIC[96]+ALGEBRAIC[98]+ALGEBRAIC[100]+ALGEBRAIC[183]+ALGEBRAIC[12]) -  2.00000*ALGEBRAIC[181])*CONSTANTS[29]*CONSTANTS[117])/( CONSTANTS[6]*CONSTANTS[118])+( ALGEBRAIC[185]*CONSTANTS[121])/CONSTANTS[118];
ALGEBRAIC[79] = ( 1.00000*ALGEBRAIC[29]*( 0.750000*STATES[6]*exp( 1.00000*ALGEBRAIC[39]) -  0.750000*CONSTANTS[3]))/(exp( 1.00000*ALGEBRAIC[39]) - 1.00000);
ALGEBRAIC[85] =  (1.00000 - ALGEBRAIC[80])*CONSTANTS[112]*ALGEBRAIC[79]*STATES[25]*( ALGEBRAIC[71]*(1.00000 - STATES[33])+ STATES[30]*ALGEBRAIC[74]*STATES[33])+ ALGEBRAIC[80]*CONSTANTS[116]*ALGEBRAIC[79]*STATES[25]*( ALGEBRAIC[75]*(1.00000 - STATES[33])+ STATES[30]*ALGEBRAIC[76]*STATES[33]);
RATES[6] = ( - ALGEBRAIC[85]*CONSTANTS[29]*CONSTANTS[117])/( CONSTANTS[6]*CONSTANTS[121]) - ALGEBRAIC[185];
ALGEBRAIC[56] =  (( CONSTANTS[4]*CONSTANTS[5])/CONSTANTS[6])*log(CONSTANTS[1]/STATES[3]);
ALGEBRAIC[59] =  CONSTANTS[41]*STATES[11]+ CONSTANTS[94]*STATES[12];
ALGEBRAIC[60] =  CONSTANTS[41]*STATES[11]+ CONSTANTS[94]*STATES[14];
ALGEBRAIC[61] = 1.00000/(1.00000+CONSTANTS[14]/ALGEBRAIC[47]);
ALGEBRAIC[62] =  CONSTANTS[42]*(STATES[0] - ALGEBRAIC[56])*pow(STATES[10], 3.00000)*( (1.00000 - ALGEBRAIC[61])*ALGEBRAIC[59]*STATES[13]+ ALGEBRAIC[61]*ALGEBRAIC[60]*STATES[15]);
ALGEBRAIC[63] = 1.00000/(1.00000+CONSTANTS[14]/ALGEBRAIC[47]);
ALGEBRAIC[64] =  CONSTANTS[96]*(STATES[0] - ALGEBRAIC[56])*STATES[16]*( (1.00000 - ALGEBRAIC[63])*STATES[17]+ ALGEBRAIC[63]*STATES[18]);
ALGEBRAIC[129] = 1.00000/(1.00000+pow(CONSTANTS[63]/STATES[9], 2.00000));
ALGEBRAIC[102] = exp(( CONSTANTS[61]*STATES[0]*CONSTANTS[6])/( CONSTANTS[4]*CONSTANTS[5]));
ALGEBRAIC[109] = 1.00000+ (CONSTANTS[1]/CONSTANTS[54])*(1.00000+1.00000/ALGEBRAIC[102]);
ALGEBRAIC[110] = CONSTANTS[1]/( CONSTANTS[54]*ALGEBRAIC[102]*ALGEBRAIC[109]);
ALGEBRAIC[113] =  ALGEBRAIC[110]*CONSTANTS[58];
ALGEBRAIC[103] = 1.00000+ (STATES[3]/CONSTANTS[54])*(1.00000+ALGEBRAIC[102]);
ALGEBRAIC[104] = ( STATES[3]*ALGEBRAIC[102])/( CONSTANTS[54]*ALGEBRAIC[103]);
ALGEBRAIC[116] =  ALGEBRAIC[104]*CONSTANTS[58];
ALGEBRAIC[106] = 1.00000+ (STATES[3]/CONSTANTS[52])*(1.00000+STATES[3]/CONSTANTS[53]);
ALGEBRAIC[107] = ( STATES[3]*STATES[3])/( ALGEBRAIC[106]*CONSTANTS[52]*CONSTANTS[53]);
ALGEBRAIC[119] =  ALGEBRAIC[107]*ALGEBRAIC[104]*CONSTANTS[56];
ALGEBRAIC[120] =  ALGEBRAIC[110]*CONSTANTS[123]*CONSTANTS[56];
ALGEBRAIC[111] = 1.00000/ALGEBRAIC[109];
ALGEBRAIC[112] =  ALGEBRAIC[111]*CONSTANTS[57];
ALGEBRAIC[114] = ALGEBRAIC[112]+ALGEBRAIC[113];
ALGEBRAIC[101] = exp(( CONSTANTS[62]*STATES[0]*CONSTANTS[6])/( CONSTANTS[4]*CONSTANTS[5]));
ALGEBRAIC[105] = 1.00000/ALGEBRAIC[103];
ALGEBRAIC[115] = ( ALGEBRAIC[105]*CONSTANTS[57])/ALGEBRAIC[101];
ALGEBRAIC[117] = ALGEBRAIC[115]+ALGEBRAIC[116];
ALGEBRAIC[108] = 1.00000/ALGEBRAIC[106];
ALGEBRAIC[118] =  ALGEBRAIC[108]*STATES[9]*CONSTANTS[59];
ALGEBRAIC[121] =  CONSTANTS[126]*ALGEBRAIC[117]*(ALGEBRAIC[119]+ALGEBRAIC[118])+ CONSTANTS[127]*ALGEBRAIC[119]*(CONSTANTS[126]+ALGEBRAIC[114]);
ALGEBRAIC[122] =  CONSTANTS[125]*ALGEBRAIC[119]*(ALGEBRAIC[117]+CONSTANTS[127])+ ALGEBRAIC[117]*ALGEBRAIC[118]*(CONSTANTS[125]+ALGEBRAIC[120]);
ALGEBRAIC[123] =  CONSTANTS[125]*ALGEBRAIC[114]*(ALGEBRAIC[119]+ALGEBRAIC[118])+ ALGEBRAIC[120]*ALGEBRAIC[118]*(CONSTANTS[126]+ALGEBRAIC[114]);
ALGEBRAIC[124] =  CONSTANTS[126]*ALGEBRAIC[120]*(ALGEBRAIC[117]+CONSTANTS[127])+ ALGEBRAIC[114]*CONSTANTS[127]*(CONSTANTS[125]+ALGEBRAIC[120]);
ALGEBRAIC[125] = ALGEBRAIC[121]/(ALGEBRAIC[121]+ALGEBRAIC[122]+ALGEBRAIC[123]+ALGEBRAIC[124]);
ALGEBRAIC[126] = ALGEBRAIC[122]/(ALGEBRAIC[121]+ALGEBRAIC[122]+ALGEBRAIC[123]+ALGEBRAIC[124]);
ALGEBRAIC[127] = ALGEBRAIC[123]/(ALGEBRAIC[121]+ALGEBRAIC[122]+ALGEBRAIC[123]+ALGEBRAIC[124]);
ALGEBRAIC[128] = ALGEBRAIC[124]/(ALGEBRAIC[121]+ALGEBRAIC[122]+ALGEBRAIC[123]+ALGEBRAIC[124]);
ALGEBRAIC[130] = ( 3.00000*( ALGEBRAIC[128]*ALGEBRAIC[119] -  ALGEBRAIC[125]*ALGEBRAIC[120])+ ALGEBRAIC[127]*ALGEBRAIC[116]) -  ALGEBRAIC[126]*ALGEBRAIC[113];
ALGEBRAIC[131] =  ALGEBRAIC[126]*CONSTANTS[126] -  ALGEBRAIC[125]*CONSTANTS[125];
ALGEBRAIC[132] =  0.800000*CONSTANTS[128]*ALGEBRAIC[129]*( CONSTANTS[7]*ALGEBRAIC[130]+ CONSTANTS[8]*ALGEBRAIC[131]);
ALGEBRAIC[184] = ( CONSTANTS[88]*ALGEBRAIC[29]*( STATES[3]*exp(ALGEBRAIC[39]) - CONSTANTS[1]))/(exp(ALGEBRAIC[39]) - 1.00000);
ALGEBRAIC[187] = (STATES[4] - STATES[3])/2.00000;
RATES[3] = ( - (ALGEBRAIC[62]+ALGEBRAIC[64]+ 3.00000*ALGEBRAIC[132]+ 3.00000*ALGEBRAIC[181]+ALGEBRAIC[184])*CONSTANTS[117]*CONSTANTS[29])/( CONSTANTS[6]*CONSTANTS[118])+( ALGEBRAIC[187]*CONSTANTS[121])/CONSTANTS[118];
ALGEBRAIC[78] = ( 1.00000*ALGEBRAIC[29]*( 0.750000*STATES[4]*exp( 1.00000*ALGEBRAIC[39]) -  0.750000*CONSTANTS[1]))/(exp( 1.00000*ALGEBRAIC[39]) - 1.00000);
ALGEBRAIC[82] =  (1.00000 - ALGEBRAIC[80])*CONSTANTS[111]*ALGEBRAIC[78]*STATES[25]*( ALGEBRAIC[71]*(1.00000 - STATES[33])+ STATES[30]*ALGEBRAIC[74]*STATES[33])+ ALGEBRAIC[80]*CONSTANTS[115]*ALGEBRAIC[78]*STATES[25]*( ALGEBRAIC[75]*(1.00000 - STATES[33])+ STATES[30]*ALGEBRAIC[76]*STATES[33]);
ALGEBRAIC[159] = 1.00000/(1.00000+pow(CONSTANTS[63]/STATES[2], 2.00000));
ALGEBRAIC[139] = 1.00000+ (CONSTANTS[1]/CONSTANTS[54])*(1.00000+1.00000/ALGEBRAIC[102]);
ALGEBRAIC[140] = CONSTANTS[1]/( CONSTANTS[54]*ALGEBRAIC[102]*ALGEBRAIC[139]);
ALGEBRAIC[143] =  ALGEBRAIC[140]*CONSTANTS[58];
ALGEBRAIC[133] = 1.00000+ (STATES[4]/CONSTANTS[54])*(1.00000+ALGEBRAIC[102]);
ALGEBRAIC[134] = ( STATES[4]*ALGEBRAIC[102])/( CONSTANTS[54]*ALGEBRAIC[133]);
ALGEBRAIC[146] =  ALGEBRAIC[134]*CONSTANTS[58];
ALGEBRAIC[136] = 1.00000+ (STATES[4]/CONSTANTS[52])*(1.00000+STATES[4]/CONSTANTS[53]);
ALGEBRAIC[137] = ( STATES[4]*STATES[4])/( ALGEBRAIC[136]*CONSTANTS[52]*CONSTANTS[53]);
ALGEBRAIC[149] =  ALGEBRAIC[137]*ALGEBRAIC[134]*CONSTANTS[56];
ALGEBRAIC[150] =  ALGEBRAIC[140]*CONSTANTS[130]*CONSTANTS[56];
ALGEBRAIC[141] = 1.00000/ALGEBRAIC[139];
ALGEBRAIC[142] =  ALGEBRAIC[141]*CONSTANTS[57];
ALGEBRAIC[144] = ALGEBRAIC[142]+ALGEBRAIC[143];
ALGEBRAIC[135] = 1.00000/ALGEBRAIC[133];
ALGEBRAIC[145] = ( ALGEBRAIC[135]*CONSTANTS[57])/ALGEBRAIC[101];
ALGEBRAIC[147] = ALGEBRAIC[145]+ALGEBRAIC[146];
ALGEBRAIC[138] = 1.00000/ALGEBRAIC[136];
ALGEBRAIC[148] =  ALGEBRAIC[138]*STATES[2]*CONSTANTS[59];
ALGEBRAIC[151] =  CONSTANTS[133]*ALGEBRAIC[147]*(ALGEBRAIC[149]+ALGEBRAIC[148])+ CONSTANTS[134]*ALGEBRAIC[149]*(CONSTANTS[133]+ALGEBRAIC[144]);
ALGEBRAIC[152] =  CONSTANTS[132]*ALGEBRAIC[149]*(ALGEBRAIC[147]+CONSTANTS[134])+ ALGEBRAIC[147]*ALGEBRAIC[148]*(CONSTANTS[132]+ALGEBRAIC[150]);
ALGEBRAIC[153] =  CONSTANTS[132]*ALGEBRAIC[144]*(ALGEBRAIC[149]+ALGEBRAIC[148])+ ALGEBRAIC[150]*ALGEBRAIC[148]*(CONSTANTS[133]+ALGEBRAIC[144]);
ALGEBRAIC[154] =  CONSTANTS[133]*ALGEBRAIC[150]*(ALGEBRAIC[147]+CONSTANTS[134])+ ALGEBRAIC[144]*CONSTANTS[134]*(CONSTANTS[132]+ALGEBRAIC[150]);
ALGEBRAIC[155] = ALGEBRAIC[151]/(ALGEBRAIC[151]+ALGEBRAIC[152]+ALGEBRAIC[153]+ALGEBRAIC[154]);
ALGEBRAIC[156] = ALGEBRAIC[152]/(ALGEBRAIC[151]+ALGEBRAIC[152]+ALGEBRAIC[153]+ALGEBRAIC[154]);
ALGEBRAIC[157] = ALGEBRAIC[153]/(ALGEBRAIC[151]+ALGEBRAIC[152]+ALGEBRAIC[153]+ALGEBRAIC[154]);
ALGEBRAIC[158] = ALGEBRAIC[154]/(ALGEBRAIC[151]+ALGEBRAIC[152]+ALGEBRAIC[153]+ALGEBRAIC[154]);
ALGEBRAIC[160] = ( 3.00000*( ALGEBRAIC[158]*ALGEBRAIC[149] -  ALGEBRAIC[155]*ALGEBRAIC[150])+ ALGEBRAIC[157]*ALGEBRAIC[146]) -  ALGEBRAIC[156]*ALGEBRAIC[143];
ALGEBRAIC[161] =  ALGEBRAIC[156]*CONSTANTS[133] -  ALGEBRAIC[155]*CONSTANTS[132];
ALGEBRAIC[162] =  0.200000*CONSTANTS[128]*ALGEBRAIC[159]*( CONSTANTS[7]*ALGEBRAIC[160]+ CONSTANTS[8]*ALGEBRAIC[161]);
RATES[4] = ( - (ALGEBRAIC[82]+ 3.00000*ALGEBRAIC[162])*CONSTANTS[29]*CONSTANTS[117])/( CONSTANTS[6]*CONSTANTS[121]) - ALGEBRAIC[187];
ALGEBRAIC[188] = ( CONSTANTS[90]*STATES[9])/(CONSTANTS[91]+STATES[9]);
ALGEBRAIC[186] = ( CONSTANTS[89]*4.00000*ALGEBRAIC[29]*( STATES[9]*exp( 2.00000*ALGEBRAIC[39]) -  0.341000*CONSTANTS[2]))/(exp( 2.00000*ALGEBRAIC[39]) - 1.00000);
RATES[0] = - (ALGEBRAIC[62]+ALGEBRAIC[64]+ALGEBRAIC[70]+ALGEBRAIC[81]+ALGEBRAIC[82]+ALGEBRAIC[85]+ALGEBRAIC[96]+ALGEBRAIC[98]+ALGEBRAIC[100]+ALGEBRAIC[132]+ALGEBRAIC[162]+ALGEBRAIC[181]+ALGEBRAIC[184]+ALGEBRAIC[183]+ALGEBRAIC[188]+ALGEBRAIC[186]+ALGEBRAIC[12]);
ALGEBRAIC[189] = (STATES[2] - STATES[9])/0.200000;
ALGEBRAIC[190] = 1.00000/(1.00000+CONSTANTS[14]/ALGEBRAIC[47]);
ALGEBRAIC[191] =  (1.00000 - ALGEBRAIC[190])*STATES[39]+ ALGEBRAIC[190]*STATES[40];
ALGEBRAIC[51] = 1.00000/(1.00000+( CONSTANTS[23]*CONSTANTS[24])/pow(CONSTANTS[24]+STATES[2], 2.00000)+( CONSTANTS[25]*CONSTANTS[26])/pow(CONSTANTS[26]+STATES[2], 2.00000));
RATES[2] =  ALGEBRAIC[51]*((( - (ALGEBRAIC[81] -  2.00000*ALGEBRAIC[162])*CONSTANTS[29]*CONSTANTS[117])/( 2.00000*CONSTANTS[6]*CONSTANTS[121])+( ALGEBRAIC[191]*CONSTANTS[120])/CONSTANTS[121]) - ALGEBRAIC[189]);
ALGEBRAIC[192] = ( CONSTANTS[108]*0.00437500*STATES[9])/(STATES[9]+0.000920000);
ALGEBRAIC[193] = ( CONSTANTS[108]*2.75000*0.00437500*STATES[9])/((STATES[9]+0.000920000) - 0.000170000);
ALGEBRAIC[194] = 1.00000/(1.00000+CONSTANTS[14]/ALGEBRAIC[47]);
ALGEBRAIC[195] = ( 0.00393750*STATES[7])/15.0000;
ALGEBRAIC[196] = ( (1.00000 - ALGEBRAIC[194])*ALGEBRAIC[192]+ ALGEBRAIC[194]*ALGEBRAIC[193]) - ALGEBRAIC[195];
ALGEBRAIC[49] = 1.00000/(1.00000+( CONSTANTS[93]*CONSTANTS[20])/pow(CONSTANTS[20]+STATES[9], 2.00000)+( CONSTANTS[21]*CONSTANTS[22])/pow(CONSTANTS[22]+STATES[9], 2.00000));
RATES[9] =  ALGEBRAIC[49]*((( - ((ALGEBRAIC[188]+ALGEBRAIC[186]) -  2.00000*ALGEBRAIC[132])*CONSTANTS[29]*CONSTANTS[117])/( 2.00000*CONSTANTS[6]*CONSTANTS[118]) - ( ALGEBRAIC[196]*CONSTANTS[119])/CONSTANTS[118])+( ALGEBRAIC[189]*CONSTANTS[121])/CONSTANTS[118]);
ALGEBRAIC[197] = (STATES[7] - STATES[8])/100.000;
RATES[7] = ALGEBRAIC[196] - ( ALGEBRAIC[197]*CONSTANTS[120])/CONSTANTS[119];
ALGEBRAIC[53] = 1.00000/(1.00000+( CONSTANTS[27]*CONSTANTS[28])/pow(CONSTANTS[28]+STATES[8], 2.00000));
RATES[8] =  ALGEBRAIC[53]*(ALGEBRAIC[197] - ALGEBRAIC[191]);
}

void Ohara_Rudy_2011::solveAnalytical(double dt)
{
STATES[0] += -dt * (ALGEBRAIC[62]+ALGEBRAIC[64]+ALGEBRAIC[70]+ALGEBRAIC[81]+ALGEBRAIC[82]+ALGEBRAIC[85]+ALGEBRAIC[96]+ALGEBRAIC[98]+ALGEBRAIC[100]+ALGEBRAIC[132]+ALGEBRAIC[162]+ALGEBRAIC[181]+ALGEBRAIC[184]+ALGEBRAIC[183]+ALGEBRAIC[188]+ALGEBRAIC[186]+ALGEBRAIC[12]);
STATES[1] +=  dt*(CONSTANTS[15]*ALGEBRAIC[45]*(ALGEBRAIC[45]+STATES[1]) -  CONSTANTS[16]*STATES[1]);
STATES[2] +=  dt*(ALGEBRAIC[51]*((( - (ALGEBRAIC[81] -  2.00000*ALGEBRAIC[162])*CONSTANTS[29]*CONSTANTS[117])/( 2.00000*CONSTANTS[6]*CONSTANTS[121])+( ALGEBRAIC[191]*CONSTANTS[120])/CONSTANTS[121]) - ALGEBRAIC[189]));
STATES[3] = dt*(( - (ALGEBRAIC[62]+ALGEBRAIC[64]+ 3.00000*ALGEBRAIC[132]+ 3.00000*ALGEBRAIC[181]+ALGEBRAIC[184])*CONSTANTS[117]*CONSTANTS[29])/( CONSTANTS[6]*CONSTANTS[118])+( ALGEBRAIC[187]*CONSTANTS[121])/CONSTANTS[118]);
STATES[4] += dt*(( - (ALGEBRAIC[82]+ 3.00000*ALGEBRAIC[162])*CONSTANTS[29]*CONSTANTS[117])/( CONSTANTS[6]*CONSTANTS[121]) - ALGEBRAIC[187]);
STATES[5] += dt*(( - ((ALGEBRAIC[70]+ALGEBRAIC[96]+ALGEBRAIC[98]+ALGEBRAIC[100]+ALGEBRAIC[183]+ALGEBRAIC[12]) -  2.00000*ALGEBRAIC[181])*CONSTANTS[29]*CONSTANTS[117])/( CONSTANTS[6]*CONSTANTS[118])+( ALGEBRAIC[185]*CONSTANTS[121])/CONSTANTS[118]);
STATES[6] += dt*(( - ALGEBRAIC[85]*CONSTANTS[29]*CONSTANTS[117])/( CONSTANTS[6]*CONSTANTS[121]) - ALGEBRAIC[185]);
STATES[7] += dt*(ALGEBRAIC[196] - ( ALGEBRAIC[197]*CONSTANTS[120])/CONSTANTS[119]);
STATES[8] += dt*(ALGEBRAIC[53]*(ALGEBRAIC[197] - ALGEBRAIC[191]));
STATES[9] +=  dt*(ALGEBRAIC[49]*((( - ((ALGEBRAIC[188]+ALGEBRAIC[186]) -  2.00000*ALGEBRAIC[132])*CONSTANTS[29]*CONSTANTS[117])/( 2.00000*CONSTANTS[6]*CONSTANTS[118]) - ( ALGEBRAIC[196]*CONSTANTS[119])/CONSTANTS[118])+( ALGEBRAIC[189]*CONSTANTS[121])/CONSTANTS[118]));
STATES[10] = ALGEBRAIC[0] - ((ALGEBRAIC[0] - STATES[10])/exp(-dt/ALGEBRAIC[13]));
STATES[11] = ALGEBRAIC[1] - ((ALGEBRAIC[1] - STATES[11])/exp(-dt/ALGEBRAIC[14]));
STATES[12] = ALGEBRAIC[1] - ((ALGEBRAIC[1] - STATES[12])/exp(-dt/ALGEBRAIC[15]));
STATES[13] = ALGEBRAIC[16] - ((ALGEBRAIC[16] - STATES[13])/exp(-dt/ALGEBRAIC[30]));
STATES[14] = ALGEBRAIC[31] - ((ALGEBRAIC[31] - STATES[14])/exp(-dt/ALGEBRAIC[40]));
STATES[15] = ALGEBRAIC[16] - ((ALGEBRAIC[16] - STATES[15])/exp(-dt/ALGEBRAIC[41]));
STATES[16] = ALGEBRAIC[32] - ((ALGEBRAIC[32] - STATES[16])/exp(-dt/ALGEBRAIC[42]));
STATES[17] = ALGEBRAIC[2] - ((ALGEBRAIC[2] - STATES[17])/exp(-dt/CONSTANTS[43]));
STATES[18] = ALGEBRAIC[3] - ((ALGEBRAIC[3] - STATES[18])/exp(-dt/CONSTANTS[95]));
STATES[19] = ALGEBRAIC[4] - ((ALGEBRAIC[4] - STATES[19])/exp(-dt/ALGEBRAIC[17]));
STATES[20] = ALGEBRAIC[5] - ((ALGEBRAIC[5] - STATES[20])/exp(-dt/ALGEBRAIC[46]));
STATES[21] = ALGEBRAIC[5] - ((ALGEBRAIC[5] - STATES[21])/exp(-dt/ALGEBRAIC[48]));
STATES[22] = ALGEBRAIC[34] - ((ALGEBRAIC[34] - STATES[22])/exp(-dt/ALGEBRAIC[17]));
STATES[23] = ALGEBRAIC[5] - ((ALGEBRAIC[5] - STATES[23])/exp(-dt/ALGEBRAIC[54]));
STATES[24] = ALGEBRAIC[5] - ((ALGEBRAIC[5] - STATES[24])/exp(-dt/ALGEBRAIC[55]));
STATES[25] = ALGEBRAIC[6] - ((ALGEBRAIC[6] - STATES[25])/exp(-dt/ALGEBRAIC[21]));
STATES[26] = ALGEBRAIC[7] - ((ALGEBRAIC[7] - STATES[26])/exp(-dt/ALGEBRAIC[22]));
STATES[27] = ALGEBRAIC[7] - ((ALGEBRAIC[7] - STATES[27])/exp(-dt/ALGEBRAIC[23]));
STATES[28] = ALGEBRAIC[19] - ((ALGEBRAIC[19] - STATES[28])/exp(-dt/ALGEBRAIC[35]));
STATES[29] = ALGEBRAIC[19] - ((ALGEBRAIC[19] - STATES[29])/exp(-dt/ALGEBRAIC[36]));
STATES[30] = ALGEBRAIC[19] - ((ALGEBRAIC[19] - STATES[30])/exp(-dt/CONSTANTS[100]));
STATES[31] = ALGEBRAIC[7] - ((ALGEBRAIC[7] - STATES[31])/exp(-dt/ALGEBRAIC[37]));
STATES[32] = ALGEBRAIC[19] - ((ALGEBRAIC[19] - STATES[32])/exp(-dt/ALGEBRAIC[44]));
STATES[33] = (ALGEBRAIC[20]*CONSTANTS[47]/ALGEBRAIC[8]) - ((ALGEBRAIC[20]*CONSTANTS[47])/(ALGEBRAIC[8]-STATES[33])*exp(-dt*ALGEBRAIC[8]));
STATES[34] = ALGEBRAIC[9] - ((ALGEBRAIC[9] - STATES[34])/exp(-dt/ALGEBRAIC[24]));
STATES[35] = ALGEBRAIC[9] - ((ALGEBRAIC[9] - STATES[35])/exp(-dt/ALGEBRAIC[25]));
STATES[36] = ALGEBRAIC[10] - ((ALGEBRAIC[10] - STATES[36])/exp(-dt/ALGEBRAIC[27]));
STATES[37] = ALGEBRAIC[26] - ((ALGEBRAIC[26] - STATES[37])/exp(-dt/ALGEBRAIC[38]));
STATES[38] = ALGEBRAIC[11] - ((ALGEBRAIC[11] - STATES[38])/exp(-dt/ALGEBRAIC[28]));
STATES[39] = ALGEBRAIC[86] - ((ALGEBRAIC[86] - STATES[39])/exp(-dt/ALGEBRAIC[92]));
STATES[40] = ALGEBRAIC[87] - ((ALGEBRAIC[87] - STATES[40])/exp(-dt/ALGEBRAIC[93]));
}