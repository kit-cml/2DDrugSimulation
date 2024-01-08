/*
   There are a total of 200 entries in the algebraic variable array.
   There are a total of 49 entries in each of the rate and state variable arrays.
   There are a total of 206 entries in the constant variable array.
 */

#include "ohara_rudy_cipa_v1_2017.hpp"
#include <cmath>
#include <cstdlib>
#include <cstdio>
/*
 * TIME is time in component environment (millisecond).
 * CONSTANTS[celltype] is celltype in component environment (dimensionless).
 * CONSTANTS[nao] is nao in component extracellular (millimolar).
 * CONSTANTS[cao] is cao in component extracellular (millimolar).
 * CONSTANTS[ko] is ko in component extracellular (millimolar).
 * CONSTANTS[R] is R in component physical_constants (joule_per_kilomole_kelvin).
 * CONSTANTS[T] is T in component physical_constants (kelvin).
 * CONSTANTS[F] is F in component physical_constants (coulomb_per_mole).
 * CONSTANTS[zna] is zna in component physical_constants (dimensionless).
 * CONSTANTS[zca] is zca in component physical_constants (dimensionless).
 * CONSTANTS[zk] is zk in component physical_constants (dimensionless).
 * CONSTANTS[L] is L in component cell_geometry (centimeter).
 * CONSTANTS[rad] is rad in component cell_geometry (centimeter).
 * CONSTANTS[vcell] is vcell in component cell_geometry (microliter).
 * CONSTANTS[Ageo] is Ageo in component cell_geometry (centimeter_squared).
 * CONSTANTS[Acap] is Acap in component cell_geometry (centimeter_squared).
 * CONSTANTS[vmyo] is vmyo in component cell_geometry (microliter).
 * CONSTANTS[vnsr] is vnsr in component cell_geometry (microliter).
 * CONSTANTS[vjsr] is vjsr in component cell_geometry (microliter).
 * CONSTANTS[vss] is vss in component cell_geometry (microliter).
 * STATES[V] is v in component membrane (millivolt).
 * ALGEBRAIC[vfrt] is vfrt in component membrane (dimensionless).
 * CONSTANTS[ffrt] is ffrt in component membrane (coulomb_per_mole_millivolt).
 * CONSTANTS[frt] is frt in component membrane (per_millivolt).
 * ALGEBRAIC[INa] is INa in component INa (microA_per_microF).
 * ALGEBRAIC[INaL] is INaL in component INaL (microA_per_microF).
 * ALGEBRAIC[Ito] is Ito in component Ito (microA_per_microF).
 * ALGEBRAIC[ICaL] is ICaL in component ICaL (microA_per_microF).
 * ALGEBRAIC[ICaNa] is ICaNa in component ICaL (microA_per_microF).
 * ALGEBRAIC[ICaK] is ICaK in component ICaL (microA_per_microF).
 * ALGEBRAIC[IKr] is IKr in component IKr (microA_per_microF).
 * ALGEBRAIC[IKs] is IKs in component IKs (microA_per_microF).
 * ALGEBRAIC[IK1] is IK1 in component IK1 (microA_per_microF).
 * ALGEBRAIC[INaCa_i] is INaCa_i in component INaCa_i (microA_per_microF).
 * ALGEBRAIC[INaCa_ss] is INaCa_ss in component INaCa_i (microA_per_microF).
 * ALGEBRAIC[INaK] is INaK in component INaK (microA_per_microF).
 * ALGEBRAIC[INab] is INab in component INab (microA_per_microF).
 * ALGEBRAIC[IKb] is IKb in component IKb (microA_per_microF).
 * ALGEBRAIC[IpCa] is IpCa in component IpCa (microA_per_microF).
 * ALGEBRAIC[ICab] is ICab in component ICab (microA_per_microF).
 * ALGEBRAIC[Istim] is Istim in component membrane (microA_per_microF).
 * CONSTANTS[stim_start] is stim_start in component membrane (millisecond).
 * CONSTANTS[stim_end] is stim_end in component membrane (millisecond).
 * CONSTANTS[amp] is amp in component membrane (microA_per_microF).
 * CONSTANTS[stim_period] is stim_period in component membrane (millisecond).
 * CONSTANTS[duration] is duration in component membrane (millisecond).
 * CONSTANTS[KmCaMK] is KmCaMK in component CaMK (millimolar).
 * CONSTANTS[aCaMK] is aCaMK in component CaMK (per_millimolar_per_millisecond).
 * CONSTANTS[bCaMK] is bCaMK in component CaMK (per_millisecond).
 * CONSTANTS[CaMKo] is CaMKo in component CaMK (dimensionless).
 * CONSTANTS[KmCaM] is KmCaM in component CaMK (millimolar).
 * ALGEBRAIC[CaMKb] is CaMKb in component CaMK (millimolar).
 * ALGEBRAIC[CaMKa] is CaMKa in component CaMK (millimolar).
 * STATES[CaMKt] is CaMKt in component CaMK (millimolar).
 * STATES[cass] is cass in component intracellular_ions (millimolar).
 * CONSTANTS[cmdnmax_b] is cmdnmax_b in component intracellular_ions (millimolar).
 * CONSTANTS[cmdnmax] is cmdnmax in component intracellular_ions (millimolar).
 * CONSTANTS[kmcmdn] is kmcmdn in component intracellular_ions (millimolar).
 * CONSTANTS[trpnmax] is trpnmax in component intracellular_ions (millimolar).
 * CONSTANTS[kmtrpn] is kmtrpn in component intracellular_ions (millimolar).
 * CONSTANTS[BSRmax] is BSRmax in component intracellular_ions (millimolar).
 * CONSTANTS[KmBSR] is KmBSR in component intracellular_ions (millimolar).
 * CONSTANTS[BSLmax] is BSLmax in component intracellular_ions (millimolar).
 * CONSTANTS[KmBSL] is KmBSL in component intracellular_ions (millimolar).
 * CONSTANTS[csqnmax] is csqnmax in component intracellular_ions (millimolar).
 * CONSTANTS[kmcsqn] is kmcsqn in component intracellular_ions (millimolar).
 * STATES[nai] is nai in component intracellular_ions (millimolar).
 * STATES[nass] is nass in component intracellular_ions (millimolar).
 * STATES[ki] is ki in component intracellular_ions (millimolar).
 * STATES[kss] is kss in component intracellular_ions (millimolar).
 * STATES[cansr] is cansr in component intracellular_ions (millimolar).
 * STATES[cajsr] is cajsr in component intracellular_ions (millimolar).
 * STATES[cai] is cai in component intracellular_ions (millimolar).
 * ALGEBRAIC[JdiffNa] is JdiffNa in component diff (millimolar_per_millisecond).
 * ALGEBRAIC[Jdiff] is Jdiff in component diff (millimolar_per_millisecond).
 * ALGEBRAIC[Jup] is Jup in component SERCA (millimolar_per_millisecond).
 * ALGEBRAIC[JdiffK] is JdiffK in component diff (millimolar_per_millisecond).
 * ALGEBRAIC[Jrel] is Jrel in component ryr (millimolar_per_millisecond).
 * ALGEBRAIC[Jtr] is Jtr in component trans_flux (millimolar_per_millisecond).
 * ALGEBRAIC[Bcai] is Bcai in component intracellular_ions (dimensionless).
 * ALGEBRAIC[Bcajsr] is Bcajsr in component intracellular_ions (dimensionless).
 * ALGEBRAIC[Bcass] is Bcass in component intracellular_ions (dimensionless).
 * CONSTANTS[cm] is cm in component intracellular_ions (microF_per_centimeter_squared).
 * CONSTANTS[PKNa] is PKNa in component reversal_potentials (dimensionless).
 * ALGEBRAIC[ENa] is ENa in component reversal_potentials (millivolt).
 * ALGEBRAIC[EK] is EK in component reversal_potentials (millivolt).
 * ALGEBRAIC[EKs] is EKs in component reversal_potentials (millivolt).
 * ALGEBRAIC[mss] is mss in component INa (dimensionless).
 * ALGEBRAIC[tm] is tm in component INa (millisecond).
 * CONSTANTS[mssV1] is mssV1 in component INa (millivolt).
 * CONSTANTS[mssV2] is mssV2 in component INa (millivolt).
 * CONSTANTS[mtV1] is mtV1 in component INa (millivolt).
 * CONSTANTS[mtV2] is mtV2 in component INa (millivolt).
 * CONSTANTS[mtD1] is mtD1 in component INa (dimensionless).
 * CONSTANTS[mtD2] is mtD2 in component INa (dimensionless).
 * CONSTANTS[mtV3] is mtV3 in component INa (millivolt).
 * CONSTANTS[mtV4] is mtV4 in component INa (millivolt).
 * STATES[m] is m in component INa (dimensionless).
 * ALGEBRAIC[hss] is hss in component INa (dimensionless).
 * ALGEBRAIC[thf] is thf in component INa (millisecond).
 * ALGEBRAIC[ths] is ths in component INa (millisecond).
 * CONSTANTS[hssV1] is hssV1 in component INa (millivolt).
 * CONSTANTS[hssV2] is hssV2 in component INa (millivolt).
 * CONSTANTS[Ahs] is Ahs in component INa (dimensionless).
 * CONSTANTS[Ahf] is Ahf in component INa (dimensionless).
 * STATES[hf] is hf in component INa (dimensionless).
 * STATES[hs] is hs in component INa (dimensionless).
 * ALGEBRAIC[h] is h in component INa (dimensionless).
 * CONSTANTS[GNa] is GNa in component INa (milliS_per_microF).
 * CONSTANTS[shift_INa_inact] is shift_INa_inact in component INa (millivolt).
 * ALGEBRAIC[jss] is jss in component INa (dimensionless).
 * ALGEBRAIC[tj] is tj in component INa (millisecond).
 * STATES[j] is j in component INa (dimensionless).
 * ALGEBRAIC[hssp] is hssp in component INa (dimensionless).
 * ALGEBRAIC[thsp] is thsp in component INa (millisecond).
 * STATES[hsp] is hsp in component INa (dimensionless).
 * ALGEBRAIC[hp] is hp in component INa (dimensionless).
 * ALGEBRAIC[tjp] is tjp in component INa (millisecond).
 * STATES[jp] is jp in component INa (dimensionless).
 * ALGEBRAIC[fINap] is fINap in component INa (dimensionless).
 * ALGEBRAIC[mLss] is mLss in component INaL (dimensionless).
 * ALGEBRAIC[tmL] is tmL in component INaL (millisecond).
 * STATES[mL] is mL in component INaL (dimensionless).
 * CONSTANTS[thL] is thL in component INaL (millisecond).
 * ALGEBRAIC[hLss] is hLss in component INaL (dimensionless).
 * STATES[hL] is hL in component INaL (dimensionless).
 * ALGEBRAIC[hLssp] is hLssp in component INaL (dimensionless).
 * CONSTANTS[thLp] is thLp in component INaL (millisecond).
 * STATES[hLp] is hLp in component INaL (dimensionless).
 * CONSTANTS[GNaL_b] is GNaL_b in component INaL (milliS_per_microF).
 * CONSTANTS[GNaL] is GNaL in component INaL (milliS_per_microF).
 * ALGEBRAIC[fINaLp] is fINaLp in component INaL (dimensionless).
 * CONSTANTS[Gto_b] is Gto_b in component Ito (milliS_per_microF).
 * ALGEBRAIC[ass] is ass in component Ito (dimensionless).
 * ALGEBRAIC[ta] is ta in component Ito (millisecond).
 * STATES[a] is a in component Ito (dimensionless).
 * ALGEBRAIC[iss] is iss in component Ito (dimensionless).
 * ALGEBRAIC[delta_epi] is delta_epi in component Ito (dimensionless).
 * ALGEBRAIC[tiF_b] is tiF_b in component Ito (millisecond).
 * ALGEBRAIC[tiS_b] is tiS_b in component Ito (millisecond).
 * ALGEBRAIC[tiF] is tiF in component Ito (millisecond).
 * ALGEBRAIC[tiS] is tiS in component Ito (millisecond).
 * ALGEBRAIC[AiF] is AiF in component Ito (dimensionless).
 * ALGEBRAIC[AiS] is AiS in component Ito (dimensionless).
 * STATES[iF] is iF in component Ito (dimensionless).
 * STATES[iS] is iS in component Ito (dimensionless).
 * ALGEBRAIC[i] is i in component Ito (dimensionless).
 * ALGEBRAIC[assp] is assp in component Ito (dimensionless).
 * STATES[ap] is ap in component Ito (dimensionless).
 * ALGEBRAIC[dti_develop] is dti_develop in component Ito (dimensionless).
 * ALGEBRAIC[dti_recover] is dti_recover in component Ito (dimensionless).
 * ALGEBRAIC[tiFp] is tiFp in component Ito (millisecond).
 * ALGEBRAIC[tiSp] is tiSp in component Ito (millisecond).
 * STATES[iFp] is iFp in component Ito (dimensionless).
 * STATES[iSp] is iSp in component Ito (dimensionless).
 * ALGEBRAIC[ip] is ip in component Ito (dimensionless).
 * CONSTANTS[Gto] is Gto in component Ito (milliS_per_microF).
 * ALGEBRAIC[fItop] is fItop in component Ito (dimensionless).
 * CONSTANTS[Kmn] is Kmn in component ICaL (millimolar).
 * CONSTANTS[k2n] is k2n in component ICaL (per_millisecond).
 * CONSTANTS[PCa_b] is PCa_b in component ICaL (dimensionless).
 * ALGEBRAIC[dss] is dss in component ICaL (dimensionless).
 * STATES[d] is d in component ICaL (dimensionless).
 * ALGEBRAIC[fss] is fss in component ICaL (dimensionless).
 * CONSTANTS[Aff] is Aff in component ICaL (dimensionless).
 * CONSTANTS[Afs] is Afs in component ICaL (dimensionless).
 * STATES[ff] is ff in component ICaL (dimensionless).
 * STATES[fs] is fs in component ICaL (dimensionless).
 * ALGEBRAIC[f] is f in component ICaL (dimensionless).
 * ALGEBRAIC[fcass] is fcass in component ICaL (dimensionless).
 * ALGEBRAIC[Afcaf] is Afcaf in component ICaL (dimensionless).
 * ALGEBRAIC[Afcas] is Afcas in component ICaL (dimensionless).
 * STATES[fcaf] is fcaf in component ICaL (dimensionless).
 * STATES[fcas] is fcas in component ICaL (dimensionless).
 * ALGEBRAIC[fca] is fca in component ICaL (dimensionless).
 * STATES[jca] is jca in component ICaL (dimensionless).
 * STATES[ffp] is ffp in component ICaL (dimensionless).
 * ALGEBRAIC[fp] is fp in component ICaL (dimensionless).
 * STATES[fcafp] is fcafp in component ICaL (dimensionless).
 * ALGEBRAIC[fcap] is fcap in component ICaL (dimensionless).
 * ALGEBRAIC[km2n] is km2n in component ICaL (per_millisecond).
 * ALGEBRAIC[anca] is anca in component ICaL (dimensionless).
 * STATES[nca] is nca in component ICaL (dimensionless).
 * ALGEBRAIC[PhiCaL] is PhiCaL in component ICaL (dimensionless).
 * ALGEBRAIC[PhiCaNa] is PhiCaNa in component ICaL (dimensionless).
 * ALGEBRAIC[PhiCaK] is PhiCaK in component ICaL (dimensionless).
 * CONSTANTS[PCa] is PCa in component ICaL (dimensionless).
 * CONSTANTS[PCap] is PCap in component ICaL (dimensionless).
 * CONSTANTS[PCaNa] is PCaNa in component ICaL (dimensionless).
 * CONSTANTS[PCaK] is PCaK in component ICaL (dimensionless).
 * CONSTANTS[PCaNap] is PCaNap in component ICaL (dimensionless).
 * CONSTANTS[PCaKp] is PCaKp in component ICaL (dimensionless).
 * ALGEBRAIC[fICaLp] is fICaLp in component ICaL (dimensionless).
 * ALGEBRAIC[td] is td in component ICaL (millisecond).
 * ALGEBRAIC[tff] is tff in component ICaL (millisecond).
 * ALGEBRAIC[tfs] is tfs in component ICaL (millisecond).
 * ALGEBRAIC[tfcaf] is tfcaf in component ICaL (millisecond).
 * ALGEBRAIC[tfcas] is tfcas in component ICaL (millisecond).
 * CONSTANTS[tjca] is tjca in component ICaL (millisecond).
 * ALGEBRAIC[tffp] is tffp in component ICaL (millisecond).
 * ALGEBRAIC[tfcafp] is tfcafp in component ICaL (millisecond).
 * CONSTANTS[v0] is v0 in component ICaL (millivolt).
 * ALGEBRAIC[A_1] is A_1 in component ICaL (dimensionless).
 * CONSTANTS[B_1] is B_1 in component ICaL (per_millivolt).
 * ALGEBRAIC[U_1] is U_1 in component ICaL (dimensionless).
 * ALGEBRAIC[A_2] is A_2 in component ICaL (dimensionless).
 * CONSTANTS[B_2] is B_2 in component ICaL (per_millivolt).
 * ALGEBRAIC[U_2] is U_2 in component ICaL (dimensionless).
 * ALGEBRAIC[A_3] is A_3 in component ICaL (dimensionless).
 * CONSTANTS[B_3] is B_3 in component ICaL (per_millivolt).
 * ALGEBRAIC[U_3] is U_3 in component ICaL (dimensionless).
 * CONSTANTS[GKr_b] is GKr_b in component IKr (milliS_per_microF).
 * STATES[IC1] is IC1 in component IKr (dimensionless).
 * STATES[IC2] is IC2 in component IKr (dimensionless).
 * STATES[C1] is C1 in component IKr (dimensionless).
 * STATES[C2] is C2 in component IKr (dimensionless).
 * STATES[O] is O in component IKr (dimensionless).
 * STATES[IO] is IO in component IKr (dimensionless).
 * STATES[IObound] is IObound in component IKr (dimensionless).
 * STATES[Obound] is Obound in component IKr (dimensionless).
 * STATES[Cbound] is Cbound in component IKr (dimensionless).
 * STATES[D] is D in component IKr (dimensionless).
 * CONSTANTS[GKr] is GKr in component IKr (milliS_per_microF).
 * CONSTANTS[A1] is A1 in component IKr (per_millisecond).
 * CONSTANTS[B1] is B1 in component IKr (per_millivolt).
 * CONSTANTS[q1] is q1 in component IKr (dimensionless).
 * CONSTANTS[A2] is A2 in component IKr (per_millisecond).
 * CONSTANTS[B2] is B2 in component IKr (per_millivolt).
 * CONSTANTS[q2] is q2 in component IKr (dimensionless).
 * CONSTANTS[A3] is A3 in component IKr (per_millisecond).
 * CONSTANTS[B3] is B3 in component IKr (per_millivolt).
 * CONSTANTS[q3] is q3 in component IKr (dimensionless).
 * CONSTANTS[A4] is A4 in component IKr (per_millisecond).
 * CONSTANTS[B4] is B4 in component IKr (per_millivolt).
 * CONSTANTS[q4] is q4 in component IKr (dimensionless).
 * CONSTANTS[A11] is A11 in component IKr (per_millisecond).
 * CONSTANTS[B11] is B11 in component IKr (per_millivolt).
 * CONSTANTS[q11] is q11 in component IKr (dimensionless).
 * CONSTANTS[A21] is A21 in component IKr (per_millisecond).
 * CONSTANTS[B21] is B21 in component IKr (per_millivolt).
 * CONSTANTS[q21] is q21 in component IKr (dimensionless).
 * CONSTANTS[A31] is A31 in component IKr (per_millisecond).
 * CONSTANTS[B31] is B31 in component IKr (per_millivolt).
 * CONSTANTS[q31] is q31 in component IKr (dimensionless).
 * CONSTANTS[A41] is A41 in component IKr (per_millisecond).
 * CONSTANTS[B41] is B41 in component IKr (per_millivolt).
 * CONSTANTS[q41] is q41 in component IKr (dimensionless).
 * CONSTANTS[A51] is A51 in component IKr (per_millisecond).
 * CONSTANTS[B51] is B51 in component IKr (per_millivolt).
 * CONSTANTS[q51] is q51 in component IKr (dimensionless).
 * CONSTANTS[A52] is A52 in component IKr (per_millisecond).
 * CONSTANTS[B52] is B52 in component IKr (per_millivolt).
 * CONSTANTS[q52] is q52 in component IKr (dimensionless).
 * CONSTANTS[A53] is A53 in component IKr (per_millisecond).
 * CONSTANTS[B53] is B53 in component IKr (per_millivolt).
 * CONSTANTS[q53] is q53 in component IKr (dimensionless).
 * CONSTANTS[A61] is A61 in component IKr (per_millisecond).
 * CONSTANTS[B61] is B61 in component IKr (per_millivolt).
 * CONSTANTS[q61] is q61 in component IKr (dimensionless).
 * CONSTANTS[A62] is A62 in component IKr (per_millisecond).
 * CONSTANTS[B62] is B62 in component IKr (per_millivolt).
 * CONSTANTS[q62] is q62 in component IKr (dimensionless).
 * CONSTANTS[A63] is A63 in component IKr (per_millisecond).
 * CONSTANTS[B63] is B63 in component IKr (per_millivolt).
 * CONSTANTS[q63] is q63 in component IKr (dimensionless).
 * CONSTANTS[Kmax] is Kmax in component IKr (dimensionless).
 * CONSTANTS[Ku] is Ku in component IKr (per_millisecond).
 * CONSTANTS[n] is n in component IKr (dimensionless).
 * CONSTANTS[halfmax] is halfmax in component IKr (dimensionless).
 * CONSTANTS[Kt] is Kt in component IKr (per_millisecond).
 * CONSTANTS[Vhalf] is Vhalf in component IKr (millivolt).
 * CONSTANTS[Temp] is Temp in component IKr (dimensionless).
 * CONSTANTS[GKs_b] is GKs_b in component IKs (milliS_per_microF).
 * CONSTANTS[GKs] is GKs in component IKs (milliS_per_microF).
 * ALGEBRAIC[xs1ss] is xs1ss in component IKs (dimensionless).
 * ALGEBRAIC[xs2ss] is xs2ss in component IKs (dimensionless).
 * ALGEBRAIC[txs1] is txs1 in component IKs (millisecond).
 * CONSTANTS[txs1_max] is txs1_max in component IKs (millisecond).
 * STATES[xs1] is xs1 in component IKs (dimensionless).
 * STATES[xs2] is xs2 in component IKs (dimensionless).
 * ALGEBRAIC[KsCa] is KsCa in component IKs (dimensionless).
 * ALGEBRAIC[txs2] is txs2 in component IKs (millisecond).
 * CONSTANTS[GK1] is GK1 in component IK1 (milliS_per_microF).
 * CONSTANTS[GK1_b] is GK1_b in component IK1 (milliS_per_microF).
 * ALGEBRAIC[xk1ss] is xk1ss in component IK1 (dimensionless).
 * ALGEBRAIC[txk1] is txk1 in component IK1 (millisecond).
 * STATES[xk1] is xk1 in component IK1 (dimensionless).
 * ALGEBRAIC[rk1] is rk1 in component IK1 (millisecond).
 * CONSTANTS[kna1] is kna1 in component INaCa_i (per_millisecond).
 * CONSTANTS[kna2] is kna2 in component INaCa_i (per_millisecond).
 * CONSTANTS[kna3] is kna3 in component INaCa_i (per_millisecond).
 * CONSTANTS[kasymm] is kasymm in component INaCa_i (dimensionless).
 * CONSTANTS[wna] is wna in component INaCa_i (dimensionless).
 * CONSTANTS[wca] is wca in component INaCa_i (dimensionless).
 * CONSTANTS[wnaca] is wnaca in component INaCa_i (dimensionless).
 * CONSTANTS[kcaon] is kcaon in component INaCa_i (per_millisecond).
 * CONSTANTS[kcaoff] is kcaoff in component INaCa_i (per_millisecond).
 * CONSTANTS[qna] is qna in component INaCa_i (dimensionless).
 * CONSTANTS[qca] is qca in component INaCa_i (dimensionless).
 * ALGEBRAIC[hna] is hna in component INaCa_i (dimensionless).
 * ALGEBRAIC[hca] is hca in component INaCa_i (dimensionless).
 * CONSTANTS[KmCaAct] is KmCaAct in component INaCa_i (millimolar).
 * CONSTANTS[Gncx_b] is Gncx_b in component INaCa_i (milliS_per_microF).
 * CONSTANTS[Gncx] is Gncx in component INaCa_i (milliS_per_microF).
 * ALGEBRAIC[h1_i] is h1_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[h2_i] is h2_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[h3_i] is h3_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[h4_i] is h4_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[h5_i] is h5_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[h6_i] is h6_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[h7_i] is h7_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[h8_i] is h8_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[h9_i] is h9_i in component INaCa_i (dimensionless).
 * CONSTANTS[h10_i] is h10_i in component INaCa_i (dimensionless).
 * CONSTANTS[h11_i] is h11_i in component INaCa_i (dimensionless).
 * CONSTANTS[h12_i] is h12_i in component INaCa_i (dimensionless).
 * CONSTANTS[k1_i] is k1_i in component INaCa_i (dimensionless).
 * CONSTANTS[k2_i] is k2_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[k3p_i] is k3p_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[k3pp_i] is k3pp_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[k3_i] is k3_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[k4_i] is k4_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[k4p_i] is k4p_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[k4pp_i] is k4pp_i in component INaCa_i (dimensionless).
 * CONSTANTS[k5_i] is k5_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[k6_i] is k6_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[k7_i] is k7_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[k8_i] is k8_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[x1_i] is x1_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[x2_i] is x2_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[x3_i] is x3_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[x4_i] is x4_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[E1_i] is E1_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[E2_i] is E2_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[E3_i] is E3_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[E4_i] is E4_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[allo_i] is allo_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[JncxNa_i] is JncxNa_i in component INaCa_i (millimolar_per_millisecond).
 * ALGEBRAIC[JncxCa_i] is JncxCa_i in component INaCa_i (millimolar_per_millisecond).
 * ALGEBRAIC[h1_ss] is h1_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[h2_ss] is h2_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[h3_ss] is h3_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[h4_ss] is h4_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[h5_ss] is h5_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[h6_ss] is h6_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[h7_ss] is h7_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[h8_ss] is h8_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[h9_ss] is h9_ss in component INaCa_i (dimensionless).
 * CONSTANTS[h10_ss] is h10_ss in component INaCa_i (dimensionless).
 * CONSTANTS[h11_ss] is h11_ss in component INaCa_i (dimensionless).
 * CONSTANTS[h12_ss] is h12_ss in component INaCa_i (dimensionless).
 * CONSTANTS[k1_ss] is k1_ss in component INaCa_i (dimensionless).
 * CONSTANTS[k2_ss] is k2_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[k3p_ss] is k3p_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[k3pp_ss] is k3pp_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[k3_ss] is k3_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[k4_ss] is k4_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[k4p_ss] is k4p_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[k4pp_ss] is k4pp_ss in component INaCa_i (dimensionless).
 * CONSTANTS[k5_ss] is k5_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[k6_ss] is k6_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[k7_ss] is k7_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[k8_ss] is k8_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[x1_ss] is x1_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[x2_ss] is x2_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[x3_ss] is x3_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[x4_ss] is x4_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[E1_ss] is E1_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[E2_ss] is E2_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[E3_ss] is E3_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[E4_ss] is E4_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[allo_ss] is allo_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[JncxNa_ss] is JncxNa_ss in component INaCa_i (millimolar_per_millisecond).
 * ALGEBRAIC[JncxCa_ss] is JncxCa_ss in component INaCa_i (millimolar_per_millisecond).
 * CONSTANTS[k1p] is k1p in component INaK (per_millisecond).
 * CONSTANTS[k1m] is k1m in component INaK (per_millisecond).
 * CONSTANTS[k2p] is k2p in component INaK (per_millisecond).
 * CONSTANTS[k2m] is k2m in component INaK (per_millisecond).
 * CONSTANTS[k3p] is k3p in component INaK (per_millisecond).
 * CONSTANTS[k3m] is k3m in component INaK (per_millisecond).
 * CONSTANTS[k4p] is k4p in component INaK (per_millisecond).
 * CONSTANTS[k4m] is k4m in component INaK (per_millisecond).
 * CONSTANTS[Knai0] is Knai0 in component INaK (millimolar).
 * CONSTANTS[Knao0] is Knao0 in component INaK (millimolar).
 * CONSTANTS[delta] is delta in component INaK (millivolt).
 * CONSTANTS[Kki] is Kki in component INaK (per_millisecond).
 * CONSTANTS[Kko] is Kko in component INaK (per_millisecond).
 * CONSTANTS[MgADP] is MgADP in component INaK (millimolar).
 * CONSTANTS[MgATP] is MgATP in component INaK (millimolar).
 * CONSTANTS[Kmgatp] is Kmgatp in component INaK (millimolar).
 * CONSTANTS[H] is H in component INaK (millimolar).
 * CONSTANTS[eP] is eP in component INaK (dimensionless).
 * CONSTANTS[Khp] is Khp in component INaK (millimolar).
 * CONSTANTS[Knap] is Knap in component INaK (millimolar).
 * CONSTANTS[Kxkur] is Kxkur in component INaK (millimolar).
 * CONSTANTS[Pnak_b] is Pnak_b in component INaK (milliS_per_microF).
 * CONSTANTS[Pnak] is Pnak in component INaK (milliS_per_microF).
 * ALGEBRAIC[Knai] is Knai in component INaK (millimolar).
 * ALGEBRAIC[Knao] is Knao in component INaK (millimolar).
 * ALGEBRAIC[P] is P in component INaK (dimensionless).
 * ALGEBRAIC[a1] is a1 in component INaK (dimensionless).
 * CONSTANTS[b1] is b1 in component INaK (dimensionless).
 * CONSTANTS[a2] is a2 in component INaK (dimensionless).
 * ALGEBRAIC[b2] is b2 in component INaK (dimensionless).
 * ALGEBRAIC[a3] is a3 in component INaK (dimensionless).
 * ALGEBRAIC[b3] is b3 in component INaK (dimensionless).
 * CONSTANTS[a4] is a4 in component INaK (dimensionless).
 * ALGEBRAIC[b4] is b4 in component INaK (dimensionless).
 * ALGEBRAIC[x1] is x1 in component INaK (dimensionless).
 * ALGEBRAIC[x2] is x2 in component INaK (dimensionless).
 * ALGEBRAIC[x3] is x3 in component INaK (dimensionless).
 * ALGEBRAIC[x4] is x4 in component INaK (dimensionless).
 * ALGEBRAIC[E1] is E1 in component INaK (dimensionless).
 * ALGEBRAIC[E2] is E2 in component INaK (dimensionless).
 * ALGEBRAIC[E3] is E3 in component INaK (dimensionless).
 * ALGEBRAIC[E4] is E4 in component INaK (dimensionless).
 * ALGEBRAIC[JnakNa] is JnakNa in component INaK (millimolar_per_millisecond).
 * ALGEBRAIC[JnakK] is JnakK in component INaK (millimolar_per_millisecond).
 * ALGEBRAIC[xkb] is xkb in component IKb (dimensionless).
 * CONSTANTS[GKb_b] is GKb_b in component IKb (milliS_per_microF).
 * CONSTANTS[GKb] is GKb in component IKb (milliS_per_microF).
 * CONSTANTS[PNab] is PNab in component INab (milliS_per_microF).
 * ALGEBRAIC[A_INab] is A_INab in component INab (microA_per_microF).
 * CONSTANTS[B_INab] is B_INab in component INab (per_millivolt).
 * CONSTANTS[v0_INab] is v0_INab in component INab (millivolt).
 * ALGEBRAIC[U_INab] is U_INab in component INab (dimensionless).
 * CONSTANTS[PCab] is PCab in component ICab (milliS_per_microF).
 * ALGEBRAIC[A_ICab] is A_ICab in component ICab (microA_per_microF).
 * CONSTANTS[B_ICab] is B_ICab in component ICab (per_millivolt).
 * CONSTANTS[v0_ICab] is v0_ICab in component ICab (millivolt).
 * ALGEBRAIC[U_ICab] is U_ICab in component ICab (dimensionless).
 * CONSTANTS[GpCa] is GpCa in component IpCa (milliS_per_microF).
 * CONSTANTS[KmCap] is KmCap in component IpCa (millimolar).
 * CONSTANTS[bt] is bt in component ryr (millisecond).
 * CONSTANTS[a_rel] is a_rel in component ryr (millisecond).
 * ALGEBRAIC[Jrel_inf] is Jrel_inf in component ryr (dimensionless).
 * ALGEBRAIC[tau_rel] is tau_rel in component ryr (millisecond).
 * ALGEBRAIC[Jrel_infp] is Jrel_infp in component ryr (dimensionless).
 * ALGEBRAIC[Jrel_temp] is Jrel_temp in component ryr (dimensionless).
 * ALGEBRAIC[tau_relp] is tau_relp in component ryr (millisecond).
 * STATES[Jrelnp] is Jrelnp in component ryr (dimensionless).
 * STATES[Jrelp] is Jrelp in component ryr (dimensionless).
 * CONSTANTS[btp] is btp in component ryr (millisecond).
 * CONSTANTS[a_relp] is a_relp in component ryr (millisecond).
 * ALGEBRAIC[Jrel_inf_temp] is Jrel_inf_temp in component ryr (dimensionless).
 * ALGEBRAIC[fJrelp] is fJrelp in component ryr (dimensionless).
 * CONSTANTS[Jrel_scaling_factor] is Jrel_scaling_factor in component ryr (dimensionless).
 * ALGEBRAIC[tau_rel_temp] is tau_rel_temp in component ryr (millisecond).
 * ALGEBRAIC[tau_relp_temp] is tau_relp_temp in component ryr (millisecond).
 * CONSTANTS[upScale] is upScale in component SERCA (dimensionless).
 * ALGEBRAIC[Jupnp] is Jupnp in component SERCA (millimolar_per_millisecond).
 * ALGEBRAIC[Jupp] is Jupp in component SERCA (millimolar_per_millisecond).
 * ALGEBRAIC[fJupp] is fJupp in component SERCA (dimensionless).
 * ALGEBRAIC[Jleak] is Jleak in component SERCA (millimolar_per_millisecond).
 * CONSTANTS[Jup_b] is Jup_b in component SERCA (dimensionless).
 * RATES[V] is d/dt v in component membrane (millivolt).
 * RATES[CaMKt] is d/dt CaMKt in component CaMK (millimolar).
 * RATES[nai] is d/dt nai in component intracellular_ions (millimolar).
 * RATES[nass] is d/dt nass in component intracellular_ions (millimolar).
 * RATES[ki] is d/dt ki in component intracellular_ions (millimolar).
 * RATES[kss] is d/dt kss in component intracellular_ions (millimolar).
 * RATES[cai] is d/dt cai in component intracellular_ions (millimolar).
 * RATES[cass] is d/dt cass in component intracellular_ions (millimolar).
 * RATES[cansr] is d/dt cansr in component intracellular_ions (millimolar).
 * RATES[cajsr] is d/dt cajsr in component intracellular_ions (millimolar).
 * RATES[m] is d/dt m in component INa (dimensionless).
 * RATES[hf] is d/dt hf in component INa (dimensionless).
 * RATES[hs] is d/dt hs in component INa (dimensionless).
 * RATES[j] is d/dt j in component INa (dimensionless).
 * RATES[hsp] is d/dt hsp in component INa (dimensionless).
 * RATES[jp] is d/dt jp in component INa (dimensionless).
 * RATES[mL] is d/dt mL in component INaL (dimensionless).
 * RATES[hL] is d/dt hL in component INaL (dimensionless).
 * RATES[hLp] is d/dt hLp in component INaL (dimensionless).
 * RATES[a] is d/dt a in component Ito (dimensionless).
 * RATES[iF] is d/dt iF in component Ito (dimensionless).
 * RATES[iS] is d/dt iS in component Ito (dimensionless).
 * RATES[ap] is d/dt ap in component Ito (dimensionless).
 * RATES[iFp] is d/dt iFp in component Ito (dimensionless).
 * RATES[iSp] is d/dt iSp in component Ito (dimensionless).
 * RATES[d] is d/dt d in component ICaL (dimensionless).
 * RATES[ff] is d/dt ff in component ICaL (dimensionless).
 * RATES[fs] is d/dt fs in component ICaL (dimensionless).
 * RATES[fcaf] is d/dt fcaf in component ICaL (dimensionless).
 * RATES[fcas] is d/dt fcas in component ICaL (dimensionless).
 * RATES[jca] is d/dt jca in component ICaL (dimensionless).
 * RATES[ffp] is d/dt ffp in component ICaL (dimensionless).
 * RATES[fcafp] is d/dt fcafp in component ICaL (dimensionless).
 * RATES[nca] is d/dt nca in component ICaL (dimensionless).
 * RATES[IC1] is d/dt IC1 in component IKr (dimensionless).
 * RATES[IC2] is d/dt IC2 in component IKr (dimensionless).
 * RATES[C1] is d/dt C1 in component IKr (dimensionless).
 * RATES[C2] is d/dt C2 in component IKr (dimensionless).
 * RATES[O] is d/dt O in component IKr (dimensionless).
 * RATES[IO] is d/dt IO in component IKr (dimensionless).
 * RATES[IObound] is d/dt IObound in component IKr (dimensionless).
 * RATES[Obound] is d/dt Obound in component IKr (dimensionless).
 * RATES[Cbound] is d/dt Cbound in component IKr (dimensionless).
 * RATES[D] is d/dt D in component IKr (dimensionless).
 * RATES[xs1] is d/dt xs1 in component IKs (dimensionless).
 * RATES[xs2] is d/dt xs2 in component IKs (dimensionless).
 * RATES[xk1] is d/dt xk1 in component IK1 (dimensionless).
 * RATES[Jrelnp] is d/dt Jrelnp in component ryr (dimensionless).
 * RATES[Jrelp] is d/dt Jrelp in component ryr (dimensionless).
 */


ohara_rudy_cipa_v1_2017::ohara_rudy_cipa_v1_2017()
{
algebraic_size = 200;
constants_size = 206+3;
states_size = 49;
ALGEBRAIC = new double[algebraic_size];
CONSTANTS = new double[constants_size];
RATES = new double[states_size];
STATES = new double[states_size];
}

ohara_rudy_cipa_v1_2017::~ohara_rudy_cipa_v1_2017()
{
delete []ALGEBRAIC;
delete []CONSTANTS;
delete []RATES;
delete []STATES;
printf("DEALLOCATE P_CELL SUCCESS!!\n");
}

void ohara_rudy_cipa_v1_2017::___initConsts()
{
/*
STATES[V] = -88.00190465;
CONSTANTS[R] = 8314;
CONSTANTS[T] = 310;
CONSTANTS[F] = 96485;
CONSTANTS[cm] = 1;
CONSTANTS[rad] = 0.0011;
CONSTANTS[L] = 0.01;
CONSTANTS[vcell] =  1000.00*3.14000*CONSTANTS[rad]*CONSTANTS[rad]*CONSTANTS[L];
CONSTANTS[stim_start] = 10;
CONSTANTS[stim_end] = 100000000000000000;
CONSTANTS[amp] = -80;
CONSTANTS[stim_period] = 1000;
CONSTANTS[duration] = 0.5;
CONSTANTS[clamp_start] = 10.;
CONSTANTS[clamp_duration] = 10.;
CONSTANTS[clamp_period] = 1000;
CONSTANTS[ffrt] =  CONSTANTS[F]*CONSTANTS[frt];
CONSTANTS[frt] = CONSTANTS[F]/( CONSTANTS[R]*CONSTANTS[T]);
#if defined CAB_ORUDY2011 
CONSTANTS[PCab] = 2.5e-8;
CONSTANTS[frt] = CONSTANTS[F]/( CONSTANTS[R]*CONSTANTS[T]);
CONSTANTS[ffrt] =  CONSTANTS[F]*CONSTANTS[frt];
CONSTANTS[v0] = 0.000000;
CONSTANTS[B_ICab] =  2.00000*CONSTANTS[frt];
CONSTANTS[cao] = 1.8;
CONSTANTS[vss] =  0.0200000*CONSTANTS[vcell];
CONSTANTS[Ageo] =  2.00000*3.14000*CONSTANTS[rad]*CONSTANTS[rad]+ 2.00000*3.14000*CONSTANTS[rad]*CONSTANTS[L];
CONSTANTS[Acap] =  2.00000*CONSTANTS[Ageo];
CONSTANTS[vmyo] =  0.680000*CONSTANTS[vcell];
CONSTANTS[vnsr] =  0.0552000*CONSTANTS[vcell];
CONSTANTS[cmdnmax_b] = 0.05;
CONSTANTS[cmdnmax] = (CONSTANTS[celltype]==1.00000 ?  CONSTANTS[cmdnmax_b]*1.30000 : CONSTANTS[cmdnmax_b]);
CONSTANTS[kmcmdn] = 0.00238;
CONSTANTS[trpnmax] = 0.07;
CONSTANTS[kmtrpn] = 0.0005;
CONSTANTS[KmCaMK] = 0.15;
CONSTANTS[upScale] = (CONSTANTS[celltype]==1.00000 ? 1.30000 : 1.00000);
CONSTANTS[CaMKo] = 0.05;
CONSTANTS[KmCaM] = 0.0015;
CONSTANTS[Jup_b] = 1.0;
CONSTANTS[Jrel_scaling_factor] = 1.0;
CONSTANTS[aCaMK] = 0.05;
CONSTANTS[bCaMK] = 0.00068;
CONSTANTS[vjsr] =  0.00480000*CONSTANTS[vcell];
CONSTANTS[csqnmax] = 10;
CONSTANTS[kmcsqn] = 0.8;
CONSTANTS[bt] = 4.75;
CONSTANTS[btp] =  1.25000*CONSTANTS[bt];
CONSTANTS[a_relp] =  0.500000*CONSTANTS[btp];
CONSTANTS[a_rel] =  0.500000*CONSTANTS[bt];
CONSTANTS[BSLmax] = 1.124;
CONSTANTS[BSRmax] = 0.047;
#endif
#if defined CAK_ORUDY2011 
CONSTANTS[bt] = 4.75;
CONSTANTS[btp] =  1.25000*CONSTANTS[bt];
CONSTANTS[a_relp] =  0.500000*CONSTANTS[btp];
CONSTANTS[a_rel] =  0.500000*CONSTANTS[bt];
CONSTANTS[frt] = CONSTANTS[F]/( CONSTANTS[R]*CONSTANTS[T]);
CONSTANTS[ffrt] =  CONSTANTS[F]*CONSTANTS[frt];
CONSTANTS[v0] = 0.000000;
CONSTANTS[B_3] = CONSTANTS[frt];
CONSTANTS[Jrel_scaling_factor] = 1.0;
CONSTANTS[Jup_b] = 1.0;
CONSTANTS[upScale] = (CONSTANTS[celltype]==1.00000 ? 1.30000 : 1.00000);
CONSTANTS[csqnmax] = 10;
CONSTANTS[vjsr] =  0.00480000*CONSTANTS[vcell];
CONSTANTS[kmcsqn] = 0.8;
CONSTANTS[cmdnmax_b] = 0.05;
CONSTANTS[cmdnmax] = (CONSTANTS[celltype]==1.00000 ?  CONSTANTS[cmdnmax_b]*1.30000 : CONSTANTS[cmdnmax_b]);
CONSTANTS[kmcmdn] = 0.00238;
CONSTANTS[trpnmax] = 0.07;
CONSTANTS[kmtrpn] = 0.0005;
CONSTANTS[vmyo] =  0.680000*CONSTANTS[vcell];
CONSTANTS[vnsr] =  0.0552000*CONSTANTS[vcell];
CONSTANTS[BSRmax] = 0.047;
CONSTANTS[KmBSR] = 0.00087;
CONSTANTS[BSLmax] = 1.124;
CONSTANTS[KmBSL] = 0.0087;
CONSTANTS[aCaMK] = 0.05;
CONSTANTS[bCaMK] = 0.00068;
CONSTANTS[CaMKo] = 0.05;
CONSTANTS[KmCaM] = 0.0015;
CONSTANTS[Kmn] = 0.002;
CONSTANTS[Ageo] =  2.00000*3.14000*CONSTANTS[rad]*CONSTANTS[rad]+ 2.00000*3.14000*CONSTANTS[rad]*CONSTANTS[L];
CONSTANTS[Acap] =  2.00000*CONSTANTS[Ageo];
CONSTANTS[vss] =  0.0200000*CONSTANTS[vcell];
CONSTANTS[vjsr] =  0.00480000*CONSTANTS[vcell];
CONSTANTS[k2n] = 1000;
CONSTANTS[tjca] = 75.0000;
CONSTANTS[Aff] = 0.600000;
CONSTANTS[Afs] = 1.00000 - CONSTANTS[Aff];
CONSTANTS[cao] = 1.8;
CONSTANTS[KmCaMK] = 0.15;
CONSTANTS[PCa_b] = 0.0001007;
CONSTANTS[PCa] = (CONSTANTS[celltype]==1.00000 ?  CONSTANTS[PCa_b]*1.20000 : CONSTANTS[celltype]==2.00000 ?  CONSTANTS[PCa_b]*2.50000 : CONSTANTS[PCa_b]);
CONSTANTS[PCap] =  1.10000*CONSTANTS[PCa];
CONSTANTS[PCaK] =  0.000357400*CONSTANTS[PCa];
CONSTANTS[PCaKp] =  0.000357400*CONSTANTS[PCap];
CONSTANTS[ko] = 5.4;
#endif
#if defined CAL_ORUDY2011 
CONSTANTS[bt] = 4.75;
CONSTANTS[btp] =  1.25000*CONSTANTS[bt];
CONSTANTS[frt] = CONSTANTS[F]/( CONSTANTS[R]*CONSTANTS[T]);
CONSTANTS[ffrt] =  CONSTANTS[F]*CONSTANTS[frt];
CONSTANTS[v0] = 0.000000;
CONSTANTS[B_1] =  2.00000*CONSTANTS[frt];
CONSTANTS[a_relp] =  0.500000*CONSTANTS[btp];
CONSTANTS[a_rel] =  0.500000*CONSTANTS[bt];
CONSTANTS[upScale] = (CONSTANTS[celltype]==1.00000 ? 1.30000 : 1.00000);
CONSTANTS[csqnmax] = 10;
CONSTANTS[vjsr] =  0.00480000*CONSTANTS[vcell];
CONSTANTS[kmcsqn] = 0.8;
CONSTANTS[cmdnmax_b] = 0.05;
CONSTANTS[cmdnmax] = (CONSTANTS[celltype]==1.00000 ?  CONSTANTS[cmdnmax_b]*1.30000 : CONSTANTS[cmdnmax_b]);
CONSTANTS[kmcmdn] = 0.00238;
CONSTANTS[trpnmax] = 0.07;
CONSTANTS[kmtrpn] = 0.0005;
CONSTANTS[vmyo] =  0.680000*CONSTANTS[vcell];
CONSTANTS[vnsr] =  0.0552000*CONSTANTS[vcell];
CONSTANTS[BSRmax] = 0.047;
CONSTANTS[KmBSR] = 0.00087;
CONSTANTS[BSLmax] = 1.124;
CONSTANTS[KmBSL] = 0.0087;
CONSTANTS[aCaMK] = 0.05;
CONSTANTS[bCaMK] = 0.00068;
CONSTANTS[CaMKo] = 0.05;
CONSTANTS[KmCaM] = 0.0015;
CONSTANTS[Kmn] = 0.002;
CONSTANTS[Ageo] =  2.00000*3.14000*CONSTANTS[rad]*CONSTANTS[rad]+ 2.00000*3.14000*CONSTANTS[rad]*CONSTANTS[L];
CONSTANTS[Acap] =  2.00000*CONSTANTS[Ageo];
CONSTANTS[vss] =  0.0200000*CONSTANTS[vcell];
CONSTANTS[vjsr] =  0.00480000*CONSTANTS[vcell];
CONSTANTS[Jrel_scaling_factor] = 1.0;
CONSTANTS[Jup_b] = 1.0;
CONSTANTS[k2n] = 1000;
CONSTANTS[tjca] = 75.0000;
CONSTANTS[Aff] = 0.600000;
CONSTANTS[Afs] = 1.00000 - CONSTANTS[Aff];
CONSTANTS[cao] = 1.8;
CONSTANTS[KmCaMK] = 0.15;
CONSTANTS[PCa_b] = 0.0001007;
CONSTANTS[PCa] = (CONSTANTS[celltype]==1.00000 ?  CONSTANTS[PCa_b]*1.20000 : CONSTANTS[celltype]==2.00000 ?  CONSTANTS[PCa_b]*2.50000 : CONSTANTS[PCa_b]);
CONSTANTS[PCap] =  1.10000*CONSTANTS[PCa];
#endif
#if defined CANA_ORUDY2011 
CONSTANTS[bt] = 4.75;
CONSTANTS[btp] =  1.25000*CONSTANTS[bt];
CONSTANTS[a_relp] =  0.500000*CONSTANTS[btp];
CONSTANTS[a_rel] =  0.500000*CONSTANTS[bt];
CONSTANTS[frt] = CONSTANTS[F]/( CONSTANTS[R]*CONSTANTS[T]);
CONSTANTS[ffrt] =  CONSTANTS[F]*CONSTANTS[frt];
CONSTANTS[v0] = 0.000000;
CONSTANTS[B_2] = CONSTANTS[frt];
CONSTANTS[Jrel_scaling_factor] = 1.0;
CONSTANTS[Jup_b] = 1.0;
CONSTANTS[upScale] = (CONSTANTS[celltype]==1.00000 ? 1.30000 : 1.00000);
CONSTANTS[csqnmax] = 10;
CONSTANTS[vjsr] =  0.00480000*CONSTANTS[vcell];
CONSTANTS[kmcsqn] = 0.8;
CONSTANTS[cmdnmax_b] = 0.05;
CONSTANTS[cmdnmax] = (CONSTANTS[celltype]==1.00000 ?  CONSTANTS[cmdnmax_b]*1.30000 : CONSTANTS[cmdnmax_b]);
CONSTANTS[kmcmdn] = 0.00238;
CONSTANTS[trpnmax] = 0.07;
CONSTANTS[kmtrpn] = 0.0005;
CONSTANTS[vmyo] =  0.680000*CONSTANTS[vcell];
CONSTANTS[vnsr] =  0.0552000*CONSTANTS[vcell];
CONSTANTS[BSRmax] = 0.047;
CONSTANTS[KmBSR] = 0.00087;
CONSTANTS[BSLmax] = 1.124;
CONSTANTS[KmBSL] = 0.0087;
CONSTANTS[aCaMK] = 0.05;
CONSTANTS[bCaMK] = 0.00068;
CONSTANTS[CaMKo] = 0.05;
CONSTANTS[KmCaM] = 0.0015;
CONSTANTS[Kmn] = 0.002;
CONSTANTS[Ageo] =  2.00000*3.14000*CONSTANTS[rad]*CONSTANTS[rad]+ 2.00000*3.14000*CONSTANTS[rad]*CONSTANTS[L];
CONSTANTS[Acap] =  2.00000*CONSTANTS[Ageo];
CONSTANTS[vss] =  0.0200000*CONSTANTS[vcell];
CONSTANTS[vjsr] =  0.00480000*CONSTANTS[vcell];
CONSTANTS[k2n] = 1000;
CONSTANTS[tjca] = 75.0000;
CONSTANTS[Aff] = 0.600000;
CONSTANTS[Afs] = 1.00000 - CONSTANTS[Aff];
CONSTANTS[cao] = 1.8;
CONSTANTS[KmCaMK] = 0.15;
CONSTANTS[PCa_b] = 0.0001007;
CONSTANTS[PCa] = (CONSTANTS[celltype]==1.00000 ?  CONSTANTS[PCa_b]*1.20000 : CONSTANTS[celltype]==2.00000 ?  CONSTANTS[PCa_b]*2.50000 : CONSTANTS[PCa_b]);
CONSTANTS[PCap] =  1.10000*CONSTANTS[PCa];
CONSTANTS[PCaNa] =  0.00125000*CONSTANTS[PCa];
CONSTANTS[PCaNap] =  0.00125000*CONSTANTS[PCap];
CONSTANTS[nao] = 140;
#endif
#if defined K1_ORUDY2011 
CONSTANTS[GK1_b] = 0.3239783999999998;
CONSTANTS[GK1] = (CONSTANTS[celltype]==1.00000 ?  CONSTANTS[GK1_b]*1.20000 : CONSTANTS[celltype]==2.00000 ?  CONSTANTS[GK1_b]*1.30000 : CONSTANTS[GK1_b]);
CONSTANTS[ko] = 5.4;
CONSTANTS[Ageo] =  2.00000*3.14000*CONSTANTS[rad]*CONSTANTS[rad]+ 2.00000*3.14000*CONSTANTS[rad]*CONSTANTS[L];
CONSTANTS[Acap] =  2.00000*CONSTANTS[Ageo];
CONSTANTS[vmyo] =  0.680000*CONSTANTS[vcell];
CONSTANTS[vss] =  0.0200000*CONSTANTS[vcell];
#endif
#if defined KB_ORUDY2011 
CONSTANTS[GKb_b] = 0.003;
CONSTANTS[GKb] = (CONSTANTS[celltype]==1.00000 ?  CONSTANTS[GKb_b]*0.600000 : CONSTANTS[GKb_b]);
CONSTANTS[ko] = 5.4;
CONSTANTS[Acap] =  2.00000*CONSTANTS[Ageo];
CONSTANTS[Ageo] =  2.00000*3.14000*CONSTANTS[rad]*CONSTANTS[rad]+ 2.00000*3.14000*CONSTANTS[rad]*CONSTANTS[L];
CONSTANTS[vmyo] =  0.680000*CONSTANTS[vcell];
CONSTANTS[vss] =  0.0200000*CONSTANTS[vcell];
#endif
#if defined KR_ORUDY2011 
CONSTANTS[GKr_b] = 0.04658545454545456;
CONSTANTS[GKr] = (CONSTANTS[celltype]==1.00000 ?  CONSTANTS[GKr_b]*1.30000 : CONSTANTS[celltype]==2.00000 ?  CONSTANTS[GKr_b]*0.800000 : CONSTANTS[GKr_b]);
CONSTANTS[Kmax] = 0;
CONSTANTS[Ku] = 0;
CONSTANTS[n] = 1;
CONSTANTS[halfmax] = 1;
CONSTANTS[ko] = 5.4;
CONSTANTS[Kt] = 0;
CONSTANTS[Vhalf] = 1;
CONSTANTS[Temp] = 37;
CONSTANTS[Ageo] =  2.00000*3.14000*CONSTANTS[rad]*CONSTANTS[rad]+ 2.00000*3.14000*CONSTANTS[rad]*CONSTANTS[L];
CONSTANTS[Acap] =  2.00000*CONSTANTS[Ageo];
CONSTANTS[vcell] =  1000.00*3.14000*CONSTANTS[rad]*CONSTANTS[rad]*CONSTANTS[L];
CONSTANTS[vmyo] =  0.680000*CONSTANTS[vcell];
CONSTANTS[vss] =  0.0200000*CONSTANTS[vcell];
CONSTANTS[A1] = 0.0264;
CONSTANTS[B1] = 4.631E-05;
CONSTANTS[q1] = 4.843;
CONSTANTS[A2] = 4.986E-06;
CONSTANTS[B2] = -0.004226;
CONSTANTS[q2] = 4.23;
CONSTANTS[A3] = 0.001214;
CONSTANTS[B3] = 0.008516;
CONSTANTS[q3] = 4.962;
CONSTANTS[A4] = 1.854E-05;
CONSTANTS[B4] = -0.04641;
CONSTANTS[q4] = 3.769;
CONSTANTS[A11] = 0.0007868;
CONSTANTS[B11] = 1.535E-08;
CONSTANTS[q11] = 4.942;
CONSTANTS[A21] = 5.455E-06;
CONSTANTS[B21] = -0.1688;
CONSTANTS[q21] = 4.156;
CONSTANTS[A31] = 0.005509;
CONSTANTS[B31] = 7.771E-09;
CONSTANTS[q31] = 4.22;
CONSTANTS[A41] = 0.001416;
CONSTANTS[B41] = -0.02877;
CONSTANTS[q41] = 1.459;
CONSTANTS[A51] = 0.4492;
CONSTANTS[B51] = 0.008595;
CONSTANTS[q51] = 5;
CONSTANTS[A52] = 0.3181;
CONSTANTS[B52] = 3.613E-08;
CONSTANTS[q52] = 4.663;
CONSTANTS[A53] = 0.149;
CONSTANTS[B53] = 0.004668;
CONSTANTS[q53] = 2.412;
CONSTANTS[A61] = 0.01241;
CONSTANTS[B61] = 0.1725;
CONSTANTS[q61] = 5.568;
CONSTANTS[A62] = 0.3226;
CONSTANTS[B62] = -0.0006575;
CONSTANTS[q62] = 5;
CONSTANTS[A63] = 0.008978;
CONSTANTS[B63] = -0.02215;
CONSTANTS[q63] = 5.682;
#endif
#if defined KS_ORUDY2011 
CONSTANTS[GKs_b] = 0.006358000000000001;
CONSTANTS[GKs] = (CONSTANTS[celltype]==1.00000 ?  CONSTANTS[GKs_b]*1.40000 : CONSTANTS[GKs_b]);
CONSTANTS[PKNa] = 0.01833;
CONSTANTS[bt] = 4.75;
CONSTANTS[btp] =  1.25000*CONSTANTS[bt];
CONSTANTS[a_relp] =  0.500000*CONSTANTS[btp];
CONSTANTS[a_rel] =  0.500000*CONSTANTS[bt];
CONSTANTS[txs1_max] = 817.3;
CONSTANTS[Jup_b] = 1.0;
CONSTANTS[ko] = 5.4;
CONSTANTS[Ageo] =  2.00000*3.14000*CONSTANTS[rad]*CONSTANTS[rad]+ 2.00000*3.14000*CONSTANTS[rad]*CONSTANTS[L];
CONSTANTS[Acap] =  2.00000*CONSTANTS[Ageo];
CONSTANTS[vmyo] =  0.680000*CONSTANTS[vcell];
CONSTANTS[vss] =  0.0200000*CONSTANTS[vcell];
CONSTANTS[KmCaM] = 0.0015;
CONSTANTS[CaMKo] = 0.05;
CONSTANTS[aCaMK] = 0.05;
CONSTANTS[bCaMK] = 0.00068;
#endif
#if defined NA_ORUDY2011 
CONSTANTS[CaMKo] = 0.05;
CONSTANTS[KmCaM] = 0.0015;
CONSTANTS[nao] = 140;
CONSTANTS[KmCaMK] = 0.15;
CONSTANTS[GNa] = 75;
CONSTANTS[mssV1] = 39.57;
CONSTANTS[mssV2] = 9.871;
CONSTANTS[mtD1] = 6.765;
CONSTANTS[mtD2] = 8.552;
CONSTANTS[mtV1] = 11.64;
CONSTANTS[mtV2] = 34.77;
CONSTANTS[mtV3] = 77.42;
CONSTANTS[mtV4] = 5.955;
CONSTANTS[hssV1] = 82.9;
CONSTANTS[hssV2] = 6.086;
CONSTANTS[aCaMK] = 0.05;
CONSTANTS[bCaMK] = 0.00068;
CONSTANTS[shift_INa_inact] = 0;
CONSTANTS[nao] = 140;
CONSTANTS[Ahf] = 0.99;
CONSTANTS[Ahs] = 1.00000 - CONSTANTS[Ahf];
CONSTANTS[Ageo] =  2.00000*3.14000*CONSTANTS[rad]*CONSTANTS[rad]+ 2.00000*3.14000*CONSTANTS[rad]*CONSTANTS[L];
CONSTANTS[Acap] =  2.00000*CONSTANTS[Ageo];
CONSTANTS[vss] =  0.0200000*CONSTANTS[vcell];
CONSTANTS[vmyo] =  0.680000*CONSTANTS[vcell];
#endif
#if defined NAB_ORUDY2011 
CONSTANTS[v0] = 0.000000;
CONSTANTS[B_INab] = CONSTANTS[frt];
CONSTANTS[PNab] = 3.75e-10;
CONSTANTS[ffrt] =  CONSTANTS[F]*CONSTANTS[frt];
CONSTANTS[v0_INab] = 0.000000;
CONSTANTS[Ageo] =  2.00000*3.14000*CONSTANTS[rad]*CONSTANTS[rad]+ 2.00000*3.14000*CONSTANTS[rad]*CONSTANTS[L];
CONSTANTS[Acap] =  2.00000*CONSTANTS[Ageo];
CONSTANTS[vmyo] =  0.680000*CONSTANTS[vcell];
CONSTANTS[vss] =  0.0200000*CONSTANTS[vcell];
#endif
#if defined NACA_I_ORUDY2011 
CONSTANTS[Gncx_b] = 0.0008;
CONSTANTS[Gncx] = (CONSTANTS[celltype]==1.00000 ?  CONSTANTS[Gncx_b]*1.10000 : CONSTANTS[celltype]==2.00000 ?  CONSTANTS[Gncx_b]*1.40000 : CONSTANTS[Gncx_b]);
CONSTANTS[zna] = 1;
CONSTANTS[zca] = 2;
CONSTANTS[bt] = 4.75;
CONSTANTS[btp] =  1.25000*CONSTANTS[bt];
CONSTANTS[a_relp] =  0.500000*CONSTANTS[btp];
CONSTANTS[a_rel] =  0.500000*CONSTANTS[bt];
CONSTANTS[upScale] = (CONSTANTS[celltype]==1.00000 ? 1.30000 : 1.00000);
CONSTANTS[Jup_b] = 1.0;
CONSTANTS[csqnmax] = 10;
CONSTANTS[kmcsqn] = 0.8;
CONSTANTS[vjsr] =  0.00480000*CONSTANTS[vcell];
CONSTANTS[vnsr] =  0.0552000*CONSTANTS[vcell];
CONSTANTS[cmdnmax_b] = 0.05;
CONSTANTS[cmdnmax] = (CONSTANTS[celltype]==1.00000 ?  CONSTANTS[cmdnmax_b]*1.30000 : CONSTANTS[cmdnmax_b]);
CONSTANTS[kmcmdn] = 0.00238;
CONSTANTS[trpnmax] = 0.07;
CONSTANTS[kmtrpn] = 0.0005;
CONSTANTS[vmyo] =  0.680000*CONSTANTS[vcell];
CONSTANTS[KmCaMK] = 0.15;
CONSTANTS[aCaMK] = 0.05;
CONSTANTS[bCaMK] = 0.00068;
CONSTANTS[CaMKo] = 0.05;
CONSTANTS[KmCaM] = 0.0015;
CONSTANTS[BSRmax] = 0.047;
CONSTANTS[KmBSR] = 0.00087;
CONSTANTS[BSLmax] = 1.124;
CONSTANTS[KmBSL] = 0.0087;
CONSTANTS[kna1] = 15;
CONSTANTS[kna2] = 5;
CONSTANTS[kna3] = 88.12;
CONSTANTS[Ageo] =  2.00000*3.14000*CONSTANTS[rad]*CONSTANTS[rad]+ 2.00000*3.14000*CONSTANTS[rad]*CONSTANTS[L];
CONSTANTS[Acap] =  2.00000*CONSTANTS[Ageo];
CONSTANTS[vss] =  0.0200000*CONSTANTS[vcell];
CONSTANTS[cao] = 1.8;
CONSTANTS[nao] = 140;
CONSTANTS[qna] = 0.5224;
CONSTANTS[qca] = 0.167;
CONSTANTS[kasymm] = 12.5;
CONSTANTS[wnaca] = 5e3;
CONSTANTS[wna] = 6e4;
CONSTANTS[wca] = 6e4;
CONSTANTS[kcaon] = 1.5e6;
CONSTANTS[kcaoff] = 5e3;
CONSTANTS[h10_i] = CONSTANTS[kasymm]+1.00000+ (CONSTANTS[nao]/CONSTANTS[kna1])*(1.00000+CONSTANTS[nao]/CONSTANTS[kna2]);
CONSTANTS[h11_i] = ( CONSTANTS[nao]*CONSTANTS[nao])/( CONSTANTS[h10_i]*CONSTANTS[kna1]*CONSTANTS[kna2]);
CONSTANTS[h12_i] = 1.00000/CONSTANTS[h10_i];
CONSTANTS[k1_i] =  CONSTANTS[h12_i]*CONSTANTS[cao]*CONSTANTS[kcaon];
CONSTANTS[k2_i] = CONSTANTS[kcaoff];
CONSTANTS[k5_i] = CONSTANTS[kcaoff];
#endif
#if defined NACA_SS_ORUDY2011 
CONSTANTS[Gncx_b] = 0.0008;
CONSTANTS[Gncx] = (CONSTANTS[celltype]==1.00000 ?  CONSTANTS[Gncx_b]*1.10000 : CONSTANTS[celltype]==2.00000 ?  CONSTANTS[Gncx_b]*1.40000 : CONSTANTS[Gncx_b]);
CONSTANTS[zna] = 1;
CONSTANTS[zca] = 2;
CONSTANTS[bt] = 4.75;
CONSTANTS[btp] =  1.25000*CONSTANTS[bt];
CONSTANTS[a_relp] =  0.500000*CONSTANTS[btp];
CONSTANTS[a_rel] =  0.500000*CONSTANTS[bt];
CONSTANTS[upScale] = (CONSTANTS[celltype]==1.00000 ? 1.30000 : 1.00000);
CONSTANTS[csqnmax] = 10;
CONSTANTS[kmcsqn] = 0.8;
CONSTANTS[vjsr] =  0.00480000*CONSTANTS[vcell];
CONSTANTS[vnsr] =  0.0552000*CONSTANTS[vcell];
CONSTANTS[Jup_b] = 1.0;
CONSTANTS[cmdnmax_b] = 0.05;
CONSTANTS[cmdnmax] = (CONSTANTS[celltype]==1.00000 ?  CONSTANTS[cmdnmax_b]*1.30000 : CONSTANTS[cmdnmax_b]);
CONSTANTS[kmcmdn] = 0.00238;
CONSTANTS[trpnmax] = 0.07;
CONSTANTS[kmtrpn] = 0.0005;
CONSTANTS[vmyo] =  0.680000*CONSTANTS[vcell];
CONSTANTS[KmCaMK] = 0.15;
CONSTANTS[aCaMK] = 0.05;
CONSTANTS[bCaMK] = 0.00068;
CONSTANTS[CaMKo] = 0.05;
CONSTANTS[KmCaM] = 0.0015;
CONSTANTS[BSRmax] = 0.047;
CONSTANTS[KmBSR] = 0.00087;
CONSTANTS[BSLmax] = 1.124;
CONSTANTS[KmBSL] = 0.0087;
CONSTANTS[kna1] = 15;
CONSTANTS[kna2] = 5;
CONSTANTS[kna3] = 88.12;
CONSTANTS[Ageo] =  2.00000*3.14000*CONSTANTS[rad]*CONSTANTS[rad]+ 2.00000*3.14000*CONSTANTS[rad]*CONSTANTS[L];
CONSTANTS[Acap] =  2.00000*CONSTANTS[Ageo];
CONSTANTS[vss] =  0.0200000*CONSTANTS[vcell];
CONSTANTS[cao] = 1.8;
CONSTANTS[nao] = 140;
CONSTANTS[qna] = 0.5224;
CONSTANTS[qca] = 0.167;
CONSTANTS[kasymm] = 12.5;
CONSTANTS[wnaca] = 5e3;
CONSTANTS[wna] = 6e4;
CONSTANTS[wca] = 6e4;
CONSTANTS[kcaon] = 1.5e6;
CONSTANTS[kcaoff] = 5e3;
CONSTANTS[h10_ss] = CONSTANTS[kasymm]+1.00000+ (CONSTANTS[nao]/CONSTANTS[kna1])*(1.00000+CONSTANTS[nao]/CONSTANTS[kna2]);
CONSTANTS[h11_ss] = ( CONSTANTS[nao]*CONSTANTS[nao])/( CONSTANTS[h10_ss]*CONSTANTS[kna1]*CONSTANTS[kna2]);
CONSTANTS[h12_ss] = 1.00000/CONSTANTS[h10_ss];
CONSTANTS[k1_ss] =  CONSTANTS[h12_ss]*CONSTANTS[cao]*CONSTANTS[kcaon];
CONSTANTS[k2_ss] = CONSTANTS[kcaoff];
CONSTANTS[k5_ss] = CONSTANTS[kcaoff];
#endif
#if defined NAK_ORUDY2011 
CONSTANTS[Pnak_b] = 30;
CONSTANTS[Pnak] = (CONSTANTS[celltype]==1.00000 ?  CONSTANTS[Pnak_b]*0.900000 : CONSTANTS[celltype]==2.00000 ?  CONSTANTS[Pnak_b]*0.700000 : CONSTANTS[Pnak_b]);
CONSTANTS[zna] = 1;
CONSTANTS[zk] = 1;
CONSTANTS[k1m] = 182.4;
CONSTANTS[MgADP] = 0.05;
CONSTANTS[b1] =  CONSTANTS[k1m]*CONSTANTS[MgADP];
CONSTANTS[k3p] = 1899;
CONSTANTS[ko] = 5.4;
CONSTANTS[Kko] = 0.3582;
CONSTANTS[nao] = 140;
CONSTANTS[k3m] = 79300;
CONSTANTS[H] = 1e-7;
CONSTANTS[MgATP] = 9.8;
CONSTANTS[Kmgatp] = 1.698e-7;
CONSTANTS[k1p] = 949.5;
CONSTANTS[Kki] = 0.5;
CONSTANTS[k4p] = 639;
CONSTANTS[k2p] = 687.2;
CONSTANTS[k3p] = 1899;
CONSTANTS[a4] = (( CONSTANTS[k4p]*CONSTANTS[MgATP])/CONSTANTS[Kmgatp])/(1.00000+CONSTANTS[MgATP]/CONSTANTS[Kmgatp]);
CONSTANTS[a2] = CONSTANTS[k2p];
CONSTANTS[b1] =  CONSTANTS[k1m]*CONSTANTS[MgADP];
CONSTANTS[Knao0] = 27.78;
CONSTANTS[delta] = -0.155;
CONSTANTS[eP] = 4.2;
CONSTANTS[Khp] = 1.698e-7;
CONSTANTS[Knap] = 224;
CONSTANTS[Kxkur] = 292;
CONSTANTS[Ageo] =  2.00000*3.14000*CONSTANTS[rad]*CONSTANTS[rad]+ 2.00000*3.14000*CONSTANTS[rad]*CONSTANTS[L];
CONSTANTS[Acap] =  2.00000*CONSTANTS[Ageo];
CONSTANTS[vmyo] =  0.680000*CONSTANTS[vcell];
CONSTANTS[vss] =  0.0200000*CONSTANTS[vcell];
CONSTANTS[Knai0] = 9.073;
CONSTANTS[k2m] = 39.4;
CONSTANTS[k4m] = 40;
#endif
#if defined NAL_ORUDY2011 
CONSTANTS[GNaL_b] = 0.019957499999999975;
CONSTANTS[GNaL] = (CONSTANTS[celltype]==1.00000 ?  CONSTANTS[GNaL_b]*0.600000 : CONSTANTS[GNaL_b]);
CONSTANTS[CaMKo] = 0.05;
CONSTANTS[KmCaM] = 0.0015;
CONSTANTS[KmCaMK] = 0.15;
CONSTANTS[nao] = 140;
CONSTANTS[Ageo] =  2.00000*3.14000*CONSTANTS[rad]*CONSTANTS[rad]+ 2.00000*3.14000*CONSTANTS[rad]*CONSTANTS[L];
CONSTANTS[Acap] =  2.00000*CONSTANTS[Ageo];
CONSTANTS[vmyo] =  0.680000*CONSTANTS[vcell];
CONSTANTS[vss] =  0.0200000*CONSTANTS[vcell];
CONSTANTS[mtD1] = 6.765;
CONSTANTS[mtD2] = 8.552;
CONSTANTS[mtV1] = 11.64;
CONSTANTS[mtV2] = 34.77;
CONSTANTS[mtV3] = 77.42;
CONSTANTS[mtV4] = 5.955;
CONSTANTS[thL] = 200;
CONSTANTS[thLp] =  3.00000*CONSTANTS[thL];
#endif
#if defined PCA_ORUDY2011 
CONSTANTS[GpCa] = 0.0005;
CONSTANTS[KmCap] = 0.0005;
CONSTANTS[vss] =  0.0200000*CONSTANTS[vcell];
CONSTANTS[Ageo] =  2.00000*3.14000*CONSTANTS[rad]*CONSTANTS[rad]+ 2.00000*3.14000*CONSTANTS[rad]*CONSTANTS[L];
CONSTANTS[Acap] =  2.00000*CONSTANTS[Ageo];
CONSTANTS[vmyo] =  0.680000*CONSTANTS[vcell];
CONSTANTS[vnsr] =  0.0552000*CONSTANTS[vcell];
CONSTANTS[Jup_b] = 1.0;
CONSTANTS[cmdnmax_b] = 0.05;
CONSTANTS[cmdnmax] = (CONSTANTS[celltype]==1.00000 ?  CONSTANTS[cmdnmax_b]*1.30000 : CONSTANTS[cmdnmax_b]);
CONSTANTS[cmdnmax] = (CONSTANTS[celltype]==1.00000 ?  CONSTANTS[cmdnmax_b]*1.30000 : CONSTANTS[cmdnmax_b]);
CONSTANTS[kmcmdn] = 0.00238;
CONSTANTS[trpnmax] = 0.07;
CONSTANTS[kmtrpn] = 0.0005;
CONSTANTS[KmCaMK] = 0.15;
CONSTANTS[upScale] = (CONSTANTS[celltype]==1.00000 ? 1.30000 : 1.00000);
CONSTANTS[CaMKo] = 0.05;
CONSTANTS[KmCaM] = 0.0015;
CONSTANTS[aCaMK] = 0.05;
CONSTANTS[bCaMK] = 0.00068;
CONSTANTS[vjsr] =  0.00480000*CONSTANTS[vcell];
CONSTANTS[csqnmax] = 10;
CONSTANTS[kmcsqn] = 0.8;
CONSTANTS[bt] = 4.75;
CONSTANTS[btp] =  1.25000*CONSTANTS[bt];
CONSTANTS[a_relp] =  0.500000*CONSTANTS[btp];
CONSTANTS[a_rel] =  0.500000*CONSTANTS[bt];
CONSTANTS[BSLmax] = 1.124;
CONSTANTS[BSRmax] = 0.047;
#endif
#if defined TO_ORUDY2011 
CONSTANTS[Gto_b] = 0.02;
CONSTANTS[Gto] = (CONSTANTS[celltype]==1.00000 ?  CONSTANTS[Gto_b]*4.00000 : CONSTANTS[celltype]==2.00000 ?  CONSTANTS[Gto_b]*4.00000 : CONSTANTS[Gto_b]);
CONSTANTS[CaMKo] = 0.05;
CONSTANTS[KmCaM] = 0.0015;
CONSTANTS[KmCaMK] = 0.15;
CONSTANTS[ko] = 5.4;
CONSTANTS[Ageo] =  2.00000*3.14000*CONSTANTS[rad]*CONSTANTS[rad]+ 2.00000*3.14000*CONSTANTS[rad]*CONSTANTS[L];
CONSTANTS[Acap] =  2.00000*CONSTANTS[Ageo];
CONSTANTS[vss] =  0.0200000*CONSTANTS[vcell];
CONSTANTS[vmyo] =  0.680000*CONSTANTS[vcell];
#endif
#if defined CAB_ORUDY2011 
STATES[cai] = 8.6e-05;
STATES[cansr] = 1.619574538;
STATES[CaMKt] = 0.0125840447;
STATES[cajsr] = 1.571234014;
STATES[Jrelnp] = 2.5e-7;
STATES[Jrelp] = 3.12e-7;
STATES[cass] = 8.49e-05;
#endif
#if defined CAK_ORUDY2011 
STATES[cass] = 8.49e-05;
STATES[d] = 2.34e-9;
STATES[CaMKt] = 0.0125840447;
STATES[ff] = 0.9999999909;
STATES[fs] = 0.9102412777;
STATES[fcaf] = 0.9999999909;
STATES[nca] = 0.002749414044;
STATES[jca] = 0.9999738312;
STATES[fcas] = 0.9998046777;
STATES[ffp] = 0.9999999909;
STATES[fcafp] = 0.9999999909;
STATES[Jrelnp] = 2.5e-7;
STATES[Jrelp] = 3.12e-7;
STATES[cai] = 8.6e-05;
STATES[cansr] = 1.619574538;
STATES[cajsr] = 1.571234014;
STATES[ki] = 144.6555918;
STATES[kss] = 144.6555651;
#endif
#if defined CAL_ORUDY2011 
STATES[cass] = 8.49e-05;
STATES[d] = 2.34e-9;
STATES[CaMKt] = 0.0125840447;
STATES[ff] = 0.9999999909;
STATES[fs] = 0.9102412777;
STATES[fcaf] = 0.9999999909;
STATES[jca] = 0.9999738312;
STATES[nca] = 0.002749414044;
STATES[fcas] = 0.9998046777;
STATES[ffp] = 0.9999999909;
STATES[fcafp] = 0.9999999909;
STATES[Jrelnp] = 2.5e-7;
STATES[Jrelp] = 3.12e-7;
STATES[cai] = 8.6e-05;
STATES[cansr] = 1.619574538;
STATES[cajsr] = 1.571234014;
#endif
#if defined CANA_ORUDY2011 
STATES[d] = 2.34e-9;
STATES[cass] = 8.49e-05;
STATES[ff] = 0.9999999909;
STATES[fs] = 0.9102412777;
STATES[nca] = 0.002749414044;
STATES[jca] = 0.9999738312;
STATES[fcaf] = 0.9999999909;
STATES[ffp] = 0.9999999909;
STATES[fcafp] = 0.9999999909;
STATES[fcas] = 0.9998046777;
STATES[CaMKt] = 0.0125840447;
STATES[nass] = 7.268089977;
STATES[nai] = 7.268004498;
STATES[cansr] = 1.619574538;
STATES[cajsr] = 1.571234014;
STATES[cai] = 8.6e-05;
STATES[Jrelnp] = 2.5e-7;
STATES[Jrelp] = 3.12e-7;
#endif
#if defined K1_ORUDY2011 
STATES[xk1] = 0.9967597594;
STATES[ki] = 144.6555918;
STATES[kss] = 144.6555651;
#endif
#if defined KB_ORUDY2011 
STATES[ki] = 144.6555918;
STATES[kss] = 144.6555651;
#endif
#if defined KR_ORUDY2011 
STATES[ki] = 144.6555918;
STATES[kss] = 144.6555651;
STATES[IC1] = 0.999637;
STATES[IC2] = 6.83208e-05;
STATES[C1] = 1.80145e-08;
STATES[C2] = 8.26619e-05;
STATES[O] = 0.00015551;
STATES[IO] = 5.67623e-05;
STATES[D] = 0;
STATES[IObound] = 0;
STATES[Obound] = 0;
STATES[Cbound] = 0;
#endif
#if defined KS_ORUDY2011 
STATES[xs1] = 0.2707758025;
STATES[xs2] = 0.0001928503426;
STATES[ki] = 144.6555918;
STATES[kss] = 144.6555651;
STATES[CaMKt] = 0.0125840447;
#endif
#if defined NA_ORUDY2011 
STATES[m] = 0.007344121102;
STATES[j] = 0.6979908432;
STATES[jp] = 0.6979245865;
STATES[hf] = 0.6981071913;
STATES[hs] = 0.6980895801;
STATES[hsp] = 0.4549485525;
STATES[nai] = 7.268004498;
STATES[nass] = 7.268089977;
STATES[CaMKt] = 0.0125840447;
#endif
#if defined NAB_ORUDY2011 
STATES[nass] = 7.268089977;
STATES[nai] = 7.268004498;
#endif
#if defined NACA_I_ORUDY2011 
STATES[cai] = 8.6e-05;
STATES[nai] = 7.268004498;
STATES[cansr] = 1.619574538;
STATES[cass] = 8.49e-05;
STATES[nass] = 7.268089977;
STATES[CaMKt] = 0.0125840447;
STATES[cajsr] = 1.571234014;
STATES[Jrelnp] = 2.5e-7;
STATES[Jrelp] = 3.12e-7;
#endif
#if defined NACA_SS_ORUDY2011 
STATES[cai] = 8.6e-05;
STATES[nai] = 7.268004498;
STATES[cansr] = 1.619574538;
STATES[cass] = 8.49e-05;
STATES[nass] = 7.268089977;
STATES[CaMKt] = 0.0125840447;
STATES[cajsr] = 1.571234014;
STATES[Jrelnp] = 2.5e-7;
STATES[Jrelp] = 3.12e-7;
#endif
#if defined NAK_ORUDY2011 
STATES[nai] = 7.268004498;
STATES[ki] = 144.6555918;
STATES[kss] = 144.6555651;
STATES[nass] = 7.268089977;
#endif
#if defined NAL_ORUDY2011 
STATES[mL] = 0.0001882617273;
STATES[hL] = 0.5008548855;
STATES[hLp] = 0.2693065357;
STATES[CaMKt] = 0.0125840447;
STATES[nai] = 7.268004498;
STATES[nass] = 7.268089977;
#endif
#if defined PCA_ORUDY2011 
STATES[cai] = 8.6e-05;
STATES[cansr] = 1.619574538;
STATES[CaMKt] = 0.0125840447;
STATES[cajsr] = 1.571234014;
STATES[Jrelnp] = 2.5e-7;
STATES[Jrelp] = 3.12e-7;
STATES[cass] = 8.49e-05;
#endif
#if defined TO_ORUDY2011 
STATES[a] = 0.001001097687;
STATES[ap] = 0.0005100862934;
STATES[ki] = 144.6555918;
STATES[iF] = 0.9995541745;
STATES[iS] = 0.5865061736;
STATES[kss] = 144.6555651;
STATES[iFp] = 0.9995541823;
STATES[iSp] = 0.6393399482;
STATES[CaMKt] = 0.0125840447;
#endif*/
CONSTANTS[0] = 0;
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
STATES[0] = -88.00190465;
CONSTANTS[12] = 10;
CONSTANTS[13] = 100000000000000000;
CONSTANTS[14] = -80;
CONSTANTS[15] = 1000;
CONSTANTS[16] = 0.5;
CONSTANTS[17] = 0.15;
CONSTANTS[18] = 0.05;
CONSTANTS[19] = 0.00068;
CONSTANTS[20] = 0.05;
CONSTANTS[21] = 0.0015;
STATES[1] = 0.0125840447;
STATES[2] = 8.49e-05;
CONSTANTS[22] = 0.05;
CONSTANTS[23] = 0.00238;
CONSTANTS[24] = 0.07;
CONSTANTS[25] = 0.0005;
CONSTANTS[26] = 0.047;
CONSTANTS[27] = 0.00087;
CONSTANTS[28] = 1.124;
CONSTANTS[29] = 0.0087;
CONSTANTS[30] = 10;
CONSTANTS[31] = 0.8;
STATES[3] = 7.268004498;
STATES[4] = 7.268089977;
STATES[5] = 144.6555918;
STATES[6] = 144.6555651;
STATES[7] = 1.619574538;
STATES[8] = 1.571234014;
STATES[9] = 8.6e-05;
CONSTANTS[32] = 1;
CONSTANTS[33] = 0.01833;
CONSTANTS[34] = 39.57;
CONSTANTS[35] = 9.871;
CONSTANTS[36] = 11.64;
CONSTANTS[37] = 34.77;
CONSTANTS[38] = 6.765;
CONSTANTS[39] = 8.552;
CONSTANTS[40] = 77.42;
CONSTANTS[41] = 5.955;
STATES[10] = 0.007344121102;
CONSTANTS[42] = 82.9;
CONSTANTS[43] = 6.086;
CONSTANTS[44] = 0.99;
STATES[11] = 0.6981071913;
STATES[12] = 0.6980895801;
CONSTANTS[45] = 75;
CONSTANTS[46] = 0;
STATES[13] = 0.6979908432;
STATES[14] = 0.4549485525;
STATES[15] = 0.6979245865;
STATES[16] = 0.0001882617273;
CONSTANTS[47] = 200;
STATES[17] = 0.5008548855;
STATES[18] = 0.2693065357;
CONSTANTS[48] = 0.019957499999999975;
CONSTANTS[49] = 0.02;
STATES[19] = 0.001001097687;
STATES[20] = 0.9995541745;
STATES[21] = 0.5865061736;
STATES[22] = 0.0005100862934;
STATES[23] = 0.9995541823;
STATES[24] = 0.6393399482;
CONSTANTS[50] = 0.002;
CONSTANTS[51] = 1000;
CONSTANTS[52] = 0.0001007;
STATES[25] = 2.34e-9;
STATES[26] = 0.9999999909;
STATES[27] = 0.9102412777;
STATES[28] = 0.9999999909;
STATES[29] = 0.9998046777;
STATES[30] = 0.9999738312;
STATES[31] = 0.9999999909;
STATES[32] = 0.9999999909;
STATES[33] = 0.002749414044;
CONSTANTS[53] = 0.04658545454545456;
STATES[34] = 0.999637;
STATES[35] = 6.83208e-05;
STATES[36] = 1.80145e-08;
STATES[37] = 8.26619e-05;
STATES[38] = 0.00015551;
STATES[39] = 5.67623e-05;
STATES[40] = 0;
STATES[41] = 0;
STATES[42] = 0;
STATES[43] = 0;
CONSTANTS[54] = 0.0264;
CONSTANTS[55] = 4.631E-05;
CONSTANTS[56] = 4.843;
CONSTANTS[57] = 4.986E-06;
CONSTANTS[58] = -0.004226;
CONSTANTS[59] = 4.23;
CONSTANTS[60] = 0.001214;
CONSTANTS[61] = 0.008516;
CONSTANTS[62] = 4.962;
CONSTANTS[63] = 1.854E-05;
CONSTANTS[64] = -0.04641;
CONSTANTS[65] = 3.769;
CONSTANTS[66] = 0.0007868;
CONSTANTS[67] = 1.535E-08;
CONSTANTS[68] = 4.942;
CONSTANTS[69] = 5.455E-06;
CONSTANTS[70] = -0.1688;
CONSTANTS[71] = 4.156;
CONSTANTS[72] = 0.005509;
CONSTANTS[73] = 7.771E-09;
CONSTANTS[74] = 4.22;
CONSTANTS[75] = 0.001416;
CONSTANTS[76] = -0.02877;
CONSTANTS[77] = 1.459;
CONSTANTS[78] = 0.4492;
CONSTANTS[79] = 0.008595;
CONSTANTS[80] = 5;
CONSTANTS[81] = 0.3181;
CONSTANTS[82] = 3.613E-08;
CONSTANTS[83] = 4.663;
CONSTANTS[84] = 0.149;
CONSTANTS[85] = 0.004668;
CONSTANTS[86] = 2.412;
CONSTANTS[87] = 0.01241;
CONSTANTS[88] = 0.1725;
CONSTANTS[89] = 5.568;
CONSTANTS[90] = 0.3226;
CONSTANTS[91] = -0.0006575;
CONSTANTS[92] = 5;
CONSTANTS[93] = 0.008978;
CONSTANTS[94] = -0.02215;
CONSTANTS[95] = 5.682;
CONSTANTS[96] = 0;
CONSTANTS[97] = 0;
CONSTANTS[98] = 1;
CONSTANTS[99] = 1;
CONSTANTS[100] = 0;
CONSTANTS[101] = 1;
CONSTANTS[102] = 37;
CONSTANTS[103] = 0.006358000000000001;
CONSTANTS[104] = 817.3;
STATES[44] = 0.2707758025;
STATES[45] = 0.0001928503426;
CONSTANTS[105] = 0.3239783999999998;
STATES[46] = 0.9967597594;
CONSTANTS[106] = 15;
CONSTANTS[107] = 5;
CONSTANTS[108] = 88.12;
CONSTANTS[109] = 12.5;
CONSTANTS[110] = 6e4;
CONSTANTS[111] = 6e4;
CONSTANTS[112] = 5e3;
CONSTANTS[113] = 1.5e6;
CONSTANTS[114] = 5e3;
CONSTANTS[115] = 0.5224;
CONSTANTS[116] = 0.167;
CONSTANTS[117] = 150e-6;
CONSTANTS[118] = 0.0008;
CONSTANTS[119] = 949.5;
CONSTANTS[120] = 182.4;
CONSTANTS[121] = 687.2;
CONSTANTS[122] = 39.4;
CONSTANTS[123] = 1899;
CONSTANTS[124] = 79300;
CONSTANTS[125] = 639;
CONSTANTS[126] = 40;
CONSTANTS[127] = 9.073;
CONSTANTS[128] = 27.78;
CONSTANTS[129] = -0.155;
CONSTANTS[130] = 0.5;
CONSTANTS[131] = 0.3582;
CONSTANTS[132] = 0.05;
CONSTANTS[133] = 9.8;
CONSTANTS[134] = 1.698e-7;
CONSTANTS[135] = 1e-7;
CONSTANTS[136] = 4.2;
CONSTANTS[137] = 1.698e-7;
CONSTANTS[138] = 224;
CONSTANTS[139] = 292;
CONSTANTS[140] = 30;
CONSTANTS[141] = 0.003;
CONSTANTS[142] = 3.75e-10;
CONSTANTS[143] = 2.5e-8;
CONSTANTS[144] = 0.0005;
CONSTANTS[145] = 0.0005;
CONSTANTS[146] = 4.75;
STATES[47] = 2.5e-7;
STATES[48] = 3.12e-7;
CONSTANTS[147] = 1.0;
CONSTANTS[148] = 1.0;
CONSTANTS[149] = CONSTANTS[6]/( CONSTANTS[4]*CONSTANTS[5]);
CONSTANTS[150] = (CONSTANTS[0]==1.00000 ?  CONSTANTS[22]*1.30000 : CONSTANTS[22]);
CONSTANTS[151] = 1.00000 - CONSTANTS[44];
CONSTANTS[152] =  3.00000*CONSTANTS[47];
CONSTANTS[153] = (CONSTANTS[0]==1.00000 ?  CONSTANTS[48]*0.600000 : CONSTANTS[48]);
CONSTANTS[154] = (CONSTANTS[0]==1.00000 ?  CONSTANTS[49]*4.00000 : CONSTANTS[0]==2.00000 ?  CONSTANTS[49]*4.00000 : CONSTANTS[49]);
CONSTANTS[155] = 0.600000;
CONSTANTS[156] = (CONSTANTS[0]==1.00000 ?  CONSTANTS[52]*1.20000 : CONSTANTS[0]==2.00000 ?  CONSTANTS[52]*2.50000 : CONSTANTS[52]);
CONSTANTS[157] = 75.0000;
CONSTANTS[158] = 0.000000;
CONSTANTS[159] = (CONSTANTS[0]==1.00000 ?  CONSTANTS[53]*1.30000 : CONSTANTS[0]==2.00000 ?  CONSTANTS[53]*0.800000 : CONSTANTS[53]);
CONSTANTS[160] = (CONSTANTS[0]==1.00000 ?  CONSTANTS[103]*1.40000 : CONSTANTS[103]);
CONSTANTS[161] = (CONSTANTS[0]==1.00000 ?  CONSTANTS[105]*1.20000 : CONSTANTS[0]==2.00000 ?  CONSTANTS[105]*1.30000 : CONSTANTS[105]);
CONSTANTS[162] =  1000.00*3.14000*CONSTANTS[11]*CONSTANTS[11]*CONSTANTS[10];
CONSTANTS[163] = (CONSTANTS[0]==1.00000 ?  CONSTANTS[141]*0.600000 : CONSTANTS[141]);
CONSTANTS[164] = 0.000000;
CONSTANTS[165] = 0.000000;
CONSTANTS[166] =  0.500000*CONSTANTS[146];
CONSTANTS[167] =  1.25000*CONSTANTS[146];
CONSTANTS[168] = (CONSTANTS[0]==1.00000 ? 1.30000 : 1.00000);
CONSTANTS[205] = 0.000000;
CONSTANTS[169] =  CONSTANTS[6]*CONSTANTS[149];
CONSTANTS[170] = 1.00000 - CONSTANTS[155];
CONSTANTS[171] =  1.10000*CONSTANTS[156];
CONSTANTS[172] =  0.00125000*CONSTANTS[156];
CONSTANTS[173] =  0.000357400*CONSTANTS[156];
CONSTANTS[174] =  2.00000*CONSTANTS[149];
CONSTANTS[175] = CONSTANTS[149];
CONSTANTS[176] = CONSTANTS[149];
CONSTANTS[177] =  2.00000*3.14000*CONSTANTS[11]*CONSTANTS[11]+ 2.00000*3.14000*CONSTANTS[11]*CONSTANTS[10];
CONSTANTS[178] = CONSTANTS[149];
CONSTANTS[179] =  2.00000*CONSTANTS[149];
CONSTANTS[180] =  0.500000*CONSTANTS[167];
CONSTANTS[181] =  0.00125000*CONSTANTS[171];
CONSTANTS[182] =  0.000357400*CONSTANTS[171];
CONSTANTS[183] =  2.00000*CONSTANTS[177];
CONSTANTS[184] =  0.680000*CONSTANTS[162];
CONSTANTS[185] =  0.0552000*CONSTANTS[162];
CONSTANTS[186] =  0.00480000*CONSTANTS[162];
CONSTANTS[187] =  0.0200000*CONSTANTS[162];
CONSTANTS[188] = CONSTANTS[109]+1.00000+ (CONSTANTS[1]/CONSTANTS[106])*(1.00000+CONSTANTS[1]/CONSTANTS[107]);
CONSTANTS[189] = ( CONSTANTS[1]*CONSTANTS[1])/( CONSTANTS[188]*CONSTANTS[106]*CONSTANTS[107]);
CONSTANTS[190] = 1.00000/CONSTANTS[188];
CONSTANTS[191] =  CONSTANTS[190]*CONSTANTS[2]*CONSTANTS[113];
CONSTANTS[192] = CONSTANTS[114];
CONSTANTS[193] = CONSTANTS[114];
CONSTANTS[194] = (CONSTANTS[0]==1.00000 ?  CONSTANTS[118]*1.10000 : CONSTANTS[0]==2.00000 ?  CONSTANTS[118]*1.40000 : CONSTANTS[118]);
CONSTANTS[195] = CONSTANTS[109]+1.00000+ (CONSTANTS[1]/CONSTANTS[106])*(1.00000+CONSTANTS[1]/CONSTANTS[107]);
CONSTANTS[196] = ( CONSTANTS[1]*CONSTANTS[1])/( CONSTANTS[195]*CONSTANTS[106]*CONSTANTS[107]);
CONSTANTS[197] = 1.00000/CONSTANTS[195];
CONSTANTS[198] =  CONSTANTS[197]*CONSTANTS[2]*CONSTANTS[113];
CONSTANTS[199] = CONSTANTS[114];
CONSTANTS[200] = CONSTANTS[114];
CONSTANTS[201] =  CONSTANTS[120]*CONSTANTS[132];
CONSTANTS[202] = CONSTANTS[121];
CONSTANTS[203] = (( CONSTANTS[125]*CONSTANTS[133])/CONSTANTS[134])/(1.00000+CONSTANTS[133]/CONSTANTS[134]);
CONSTANTS[204] = (CONSTANTS[0]==1.00000 ?  CONSTANTS[140]*0.900000 : CONSTANTS[0]==2.00000 ?  CONSTANTS[140]*0.700000 : CONSTANTS[140]);
}

void ohara_rudy_cipa_v1_2017::initConsts()
{
CONSTANTS[celltype] = 0;
___initConsts();
}

void ohara_rudy_cipa_v1_2017::initConsts(double type)
{
CONSTANTS[celltype] = type;
___initConsts();
}

void ohara_rudy_cipa_v1_2017::initConsts(double type, double conc, double *hill)
{
initConsts(type);
STATES[D] = conc;
CONSTANTS[PCa]  *= (hill[0] > 10E-14 && hill[1] > 10E-14) ?
                  pow(1.+pow(STATES[D]/hill[0], hill[1]), -1) : 1.;
CONSTANTS[GK1] *= (hill[2] > 10E-14 && hill[3] > 10E-14) ?
                  pow(1.+pow(STATES[D]/hill[2], hill[3]), -1) : 1.;
CONSTANTS[GKs] *= (hill[4] > 10E-14 && hill[5] > 10E-14) ?
                  pow(1.+pow(STATES[D]/hill[4], hill[5]), -1) : 1.;
CONSTANTS[GNa]  *= (hill[6] > 10E-14 && hill[7] > 10E-14) ?
                  pow(1.+pow(STATES[D]/hill[6], hill[7]), -1) : 1.;
CONSTANTS[GNaL]  *= (hill[8] > 10E-14 && hill[9] > 10E-14) ?
                  pow(1.+pow(STATES[D]/hill[8], hill[9]), -1) : 1.;
CONSTANTS[Gto]  *= (hill[10] > 10E-14 && hill[11] > 10E-14) ?
                  pow(1.+pow(STATES[D]/hill[10], hill[11]), -1) : 1.;
CONSTANTS[GKr] *= (hill[12] > 10E-14 && hill[13] > 10E-14) ?
                  pow(1.+pow(STATES[D]/hill[12], hill[13]), -1) : 1.;
CONSTANTS[Kmax] = (hill[14] > 10E-14) ? hill[14] : CONSTANTS[Kmax];
CONSTANTS[Ku] = (hill[15] > 10E-14) ? hill[15] : CONSTANTS[Ku];
CONSTANTS[n] = (hill[16] > 10E-14) ? hill[16] : CONSTANTS[n];
CONSTANTS[halfmax] = (hill[17] > 10E-14) ? hill[17] : CONSTANTS[halfmax];
CONSTANTS[Vhalf] = (hill[18] > 10E-14) ? hill[18] : CONSTANTS[Vhalf];
}


void ohara_rudy_cipa_v1_2017::computeRates( double TIME, double *CONSTANTS, double *RATES, double *STATES, double *ALGEBRAIC )
{
/*
#if defined SINGLE_CELL
ALGEBRAIC[Istim] = (TIME>=CONSTANTS[stim_start]&&TIME<=CONSTANTS[stim_end]&&(TIME - CONSTANTS[stim_start]) -  floor((TIME - CONSTANTS[stim_start])/CONSTANTS[stim_period])*CONSTANTS[stim_period]<=CONSTANTS[duration] ? CONSTANTS[amp] : 0.000000);
ALGEBRAIC[vfrt] =  STATES[V]*CONSTANTS[frt];
#else
STATES[V] = (TIME>=CONSTANTS[clamp_start]&&(TIME - CONSTANTS[clamp_start]) -  floor((TIME - CONSTANTS[clamp_start])/CONSTANTS[clamp_period])*CONSTANTS[clamp_period]<=CONSTANTS[clamp_duration] ? 0. : -80.);
#endif
#if defined CAB_ORUDY2011 
ALGEBRAIC[vfrt] =  STATES[V]*CONSTANTS[frt];
ALGEBRAIC[A_ICab] = ( CONSTANTS[PCab]*4.00000*CONSTANTS[ffrt]*( STATES[cai]*exp( 2.00000*ALGEBRAIC[vfrt]) -  0.341000*CONSTANTS[cao]))/CONSTANTS[B_ICab];
ALGEBRAIC[U_ICab] =  CONSTANTS[B_ICab]*(STATES[V] - CONSTANTS[v0_ICab]);
ALGEBRAIC[ICab] = (- 1.00000e-07<=ALGEBRAIC[U_ICab]&&ALGEBRAIC[U_ICab]<=1.00000e-07 ?  ALGEBRAIC[A_ICab]*(1.00000 -  0.500000*ALGEBRAIC[U_ICab]) : ( ALGEBRAIC[A_ICab]*ALGEBRAIC[U_ICab])/(exp(ALGEBRAIC[U_ICab]) - 1.00000));
ALGEBRAIC[Jdiff] = (STATES[cass] - STATES[cai])/0.200000;
ALGEBRAIC[Bcai] = 1.00000/(1.00000+( CONSTANTS[cmdnmax]*CONSTANTS[kmcmdn])/pow(CONSTANTS[kmcmdn]+STATES[cai], 2.00000)+( CONSTANTS[trpnmax]*CONSTANTS[kmtrpn])/pow(CONSTANTS[kmtrpn]+STATES[cai], 2.00000));
ALGEBRAIC[CaMKb] = ( CONSTANTS[CaMKo]*(1.00000 - STATES[CaMKt]))/(1.00000+CONSTANTS[KmCaM]/STATES[cass]);
ALGEBRAIC[CaMKa] = ALGEBRAIC[CaMKb]+STATES[CaMKt];
ALGEBRAIC[fJupp] = 1.00000/(1.00000+CONSTANTS[KmCaMK]/ALGEBRAIC[CaMKa]);
ALGEBRAIC[Jupnp] = ( CONSTANTS[upScale]*0.00437500*STATES[cai])/(STATES[cai]+0.000920000);
ALGEBRAIC[Jupp] = ( CONSTANTS[upScale]*2.75000*0.00437500*STATES[cai])/((STATES[cai]+0.000920000) - 0.000170000);
ALGEBRAIC[Jleak] = ( 0.00393750*STATES[cansr])/15.0000;
ALGEBRAIC[Jup] =  CONSTANTS[Jup_b]*(( (1.00000 - ALGEBRAIC[fJupp])*ALGEBRAIC[Jupnp]+ ALGEBRAIC[fJupp]*ALGEBRAIC[Jupp]) - ALGEBRAIC[Jleak]);
ALGEBRAIC[Bcajsr] = 1.00000/(1.00000+( CONSTANTS[csqnmax]*CONSTANTS[kmcsqn])/pow(CONSTANTS[kmcsqn]+STATES[cajsr], 2.00000));
ALGEBRAIC[Bcass] = 1.00000/(1.00000+( CONSTANTS[BSRmax]*CONSTANTS[KmBSR])/pow(CONSTANTS[KmBSR]+STATES[cass], 2.00000)+( CONSTANTS[BSLmax]*CONSTANTS[KmBSL])/pow(CONSTANTS[KmBSL]+STATES[cass], 2.00000));
ALGEBRAIC[fJrelp] = 1.00000/(1.00000+CONSTANTS[KmCaMK]/ALGEBRAIC[CaMKa]);
ALGEBRAIC[Jrel_inf_temp] = ( CONSTANTS[a_rel]*- ALGEBRAIC[ICaL])/(1.00000+ 1.00000*pow(1.50000/STATES[cajsr], 8.00000));
ALGEBRAIC[tau_relp_temp] = CONSTANTS[btp]/(1.00000+0.0123000/STATES[cajsr]);
ALGEBRAIC[Jrel_temp] = ( CONSTANTS[a_relp]*- ALGEBRAIC[ICaL])/(1.00000+pow(1.50000/STATES[cajsr], 8.00000));
ALGEBRAIC[tau_rel_temp] = CONSTANTS[bt]/(1.00000+0.0123000/STATES[cajsr]);
ALGEBRAIC[tau_rel] = (ALGEBRAIC[tau_rel_temp]<0.00100000 ? 0.00100000 : ALGEBRAIC[tau_rel_temp]);
ALGEBRAIC[Jrel_inf] = (CONSTANTS[celltype]==2.00000 ?  ALGEBRAIC[Jrel_inf_temp]*1.70000 : ALGEBRAIC[Jrel_inf_temp]);
ALGEBRAIC[Jrel_infp] = (CONSTANTS[celltype]==2.00000 ?  ALGEBRAIC[Jrel_temp]*1.70000 : ALGEBRAIC[Jrel_temp]);
ALGEBRAIC[tau_relp] = (ALGEBRAIC[tau_relp_temp]<0.00100000 ? 0.00100000 : ALGEBRAIC[tau_relp_temp]);
ALGEBRAIC[Jrel] =  CONSTANTS[Jrel_scaling_factor]*( (1.00000 - ALGEBRAIC[fJrelp])*STATES[Jrelnp]+ ALGEBRAIC[fJrelp]*STATES[Jrelp]);
ALGEBRAIC[Jtr] = (STATES[cansr] - STATES[cajsr])/100.000;
#endif
#if defined CAK_ORUDY2011 
ALGEBRAIC[CaMKb] = ( CONSTANTS[CaMKo]*(1.00000 - STATES[CaMKt]))/(1.00000+CONSTANTS[KmCaM]/STATES[cass]);
ALGEBRAIC[CaMKa] = ALGEBRAIC[CaMKb]+STATES[CaMKt];
ALGEBRAIC[vfrt] =  STATES[V]*CONSTANTS[frt];
ALGEBRAIC[U_3] =  CONSTANTS[B_3]*(STATES[V] - CONSTANTS[v0]);
ALGEBRAIC[A_3] = ( 0.750000*CONSTANTS[ffrt]*( STATES[kss]*exp(ALGEBRAIC[vfrt]) - CONSTANTS[ko]))/CONSTANTS[B_3];
ALGEBRAIC[fICaLp] = 1.00000/(1.00000+CONSTANTS[KmCaMK]/ALGEBRAIC[CaMKa]);
ALGEBRAIC[PhiCaK] = (- 1.00000e-07<=ALGEBRAIC[U_3]&&ALGEBRAIC[U_3]<=1.00000e-07 ?  ALGEBRAIC[A_3]*(1.00000 -  0.500000*ALGEBRAIC[U_3]) : ( ALGEBRAIC[A_3]*ALGEBRAIC[U_3])/(exp(ALGEBRAIC[U_3]) - 1.00000));
ALGEBRAIC[f] =  CONSTANTS[Aff]*STATES[ff]+ CONSTANTS[Afs]*STATES[fs];
ALGEBRAIC[fp] =  CONSTANTS[Aff]*STATES[ffp]+ CONSTANTS[Afs]*STATES[fs];
ALGEBRAIC[tfcaf] = 7.00000+1.00000/( 0.0400000*exp(- (STATES[V] - 4.00000)/7.00000)+ 0.0400000*exp((STATES[V] - 4.00000)/7.00000));
ALGEBRAIC[tfcafp] =  2.50000*ALGEBRAIC[tfcaf];
ALGEBRAIC[Afcaf] = 0.300000+0.600000/(1.00000+exp((STATES[V] - 10.0000)/10.0000));
ALGEBRAIC[Afcas] = 1.00000 - ALGEBRAIC[Afcaf];
ALGEBRAIC[fca] =  ALGEBRAIC[Afcaf]*STATES[fcaf]+ ALGEBRAIC[Afcas]*STATES[fcas];
ALGEBRAIC[fcap] =  ALGEBRAIC[Afcaf]*STATES[fcafp]+ ALGEBRAIC[Afcas]*STATES[fcas];
ALGEBRAIC[ICaK] =  (1.00000 - ALGEBRAIC[fICaLp])*CONSTANTS[PCaK]*ALGEBRAIC[PhiCaK]*STATES[d]*( ALGEBRAIC[f]*(1.00000 - STATES[nca])+ STATES[jca]*ALGEBRAIC[fca]*STATES[nca])+ ALGEBRAIC[fICaLp]*CONSTANTS[PCaKp]*ALGEBRAIC[PhiCaK]*STATES[d]*( ALGEBRAIC[fp]*(1.00000 - STATES[nca])+ STATES[jca]*ALGEBRAIC[fcap]*STATES[nca]);
ALGEBRAIC[dss] = 1.00000/(1.00000+exp(- (STATES[V]+3.94000)/4.23000));
ALGEBRAIC[km2n] =  STATES[jca]*1.00000;
ALGEBRAIC[anca] = 1.00000/(CONSTANTS[k2n]/ALGEBRAIC[km2n]+pow(1.00000+CONSTANTS[Kmn]/STATES[cass], 4.00000));
ALGEBRAIC[fss] = 1.00000/(1.00000+exp((STATES[V]+19.5800)/3.69600));
ALGEBRAIC[tff] = 7.00000+1.00000/( 0.00450000*exp(- (STATES[V]+20.0000)/10.0000)+ 0.00450000*exp((STATES[V]+20.0000)/10.0000));
ALGEBRAIC[tfs] = 1000.00+1.00000/( 3.50000e-05*exp(- (STATES[V]+5.00000)/4.00000)+ 3.50000e-05*exp((STATES[V]+5.00000)/6.00000));
ALGEBRAIC[td] = 0.600000+1.00000/(exp( - 0.0500000*(STATES[V]+6.00000))+exp( 0.0900000*(STATES[V]+14.0000)));
ALGEBRAIC[ths] = 1.00000/( 0.00979400*exp(- ((STATES[V]+17.9500) - CONSTANTS[shift_INa_inact])/28.0500)+ 0.334300*exp(((STATES[V]+5.73000) - CONSTANTS[shift_INa_inact])/56.6600));
ALGEBRAIC[tffp] =  2.50000*ALGEBRAIC[tff];
ALGEBRAIC[fcass] = ALGEBRAIC[fss];
ALGEBRAIC[tfcas] = 100.000+1.00000/( 0.000120000*exp(- STATES[V]/3.00000)+ 0.000120000*exp(STATES[V]/7.00000));
ALGEBRAIC[Bcass] = 1.00000/(1.00000+( CONSTANTS[BSRmax]*CONSTANTS[KmBSR])/pow(CONSTANTS[KmBSR]+STATES[cass], 2.00000)+( CONSTANTS[BSLmax]*CONSTANTS[KmBSL])/pow(CONSTANTS[KmBSL]+STATES[cass], 2.00000));
ALGEBRAIC[fJrelp] = 1.00000/(1.00000+CONSTANTS[KmCaMK]/ALGEBRAIC[CaMKa]);
ALGEBRAIC[Jrel] =  CONSTANTS[Jrel_scaling_factor]*( (1.00000 - ALGEBRAIC[fJrelp])*STATES[Jrelnp]+ ALGEBRAIC[fJrelp]*STATES[Jrelp]);
ALGEBRAIC[Jdiff] = (STATES[cass] - STATES[cai])/0.200000;
ALGEBRAIC[JdiffK] = (STATES[kss] - STATES[ki])/2.00000;
ALGEBRAIC[Jrel_temp] = ( CONSTANTS[a_relp]*- ALGEBRAIC[ICaL])/(1.00000+pow(1.50000/STATES[cajsr], 8.00000));
ALGEBRAIC[Jrel_inf_temp] = ( CONSTANTS[a_rel]*- ALGEBRAIC[ICaL])/(1.00000+ 1.00000*pow(1.50000/STATES[cajsr], 8.00000));
ALGEBRAIC[tau_rel_temp] = CONSTANTS[bt]/(1.00000+0.0123000/STATES[cajsr]);
ALGEBRAIC[tau_relp_temp] = CONSTANTS[btp]/(1.00000+0.0123000/STATES[cajsr]);
ALGEBRAIC[fJupp] = 1.00000/(1.00000+CONSTANTS[KmCaMK]/ALGEBRAIC[CaMKa]);
ALGEBRAIC[Jupnp] = ( CONSTANTS[upScale]*0.00437500*STATES[cai])/(STATES[cai]+0.000920000);
ALGEBRAIC[Jupp] = ( CONSTANTS[upScale]*2.75000*0.00437500*STATES[cai])/((STATES[cai]+0.000920000) - 0.000170000);
ALGEBRAIC[Jleak] = ( 0.00393750*STATES[cansr])/15.0000;
ALGEBRAIC[Jrel_inf] = (CONSTANTS[celltype]==2.00000 ?  ALGEBRAIC[Jrel_inf_temp]*1.70000 : ALGEBRAIC[Jrel_inf_temp]);
ALGEBRAIC[tau_rel] = (ALGEBRAIC[tau_rel_temp]<0.00100000 ? 0.00100000 : ALGEBRAIC[tau_rel_temp]);
ALGEBRAIC[Jrel_infp] = (CONSTANTS[celltype]==2.00000 ?  ALGEBRAIC[Jrel_temp]*1.70000 : ALGEBRAIC[Jrel_temp]);
ALGEBRAIC[tau_relp] = (ALGEBRAIC[tau_relp_temp]<0.00100000 ? 0.00100000 : ALGEBRAIC[tau_relp_temp]);
ALGEBRAIC[Bcai] = 1.00000/(1.00000+( CONSTANTS[cmdnmax]*CONSTANTS[kmcmdn])/pow(CONSTANTS[kmcmdn]+STATES[cai], 2.00000)+( CONSTANTS[trpnmax]*CONSTANTS[kmtrpn])/pow(CONSTANTS[kmtrpn]+STATES[cai], 2.00000));
ALGEBRAIC[Jup] =  CONSTANTS[Jup_b]*(( (1.00000 - ALGEBRAIC[fJupp])*ALGEBRAIC[Jupnp]+ ALGEBRAIC[fJupp]*ALGEBRAIC[Jupp]) - ALGEBRAIC[Jleak]);
ALGEBRAIC[Bcajsr] = 1.00000/(1.00000+( CONSTANTS[csqnmax]*CONSTANTS[kmcsqn])/pow(CONSTANTS[kmcsqn]+STATES[cajsr], 2.00000));
ALGEBRAIC[Jtr] = (STATES[cansr] - STATES[cajsr])/100.000;
#endif
#if defined CAL_ORUDY2011 
ALGEBRAIC[CaMKb] = ( CONSTANTS[CaMKo]*(1.00000 - STATES[CaMKt]))/(1.00000+CONSTANTS[KmCaM]/STATES[cass]);
ALGEBRAIC[CaMKa] = ALGEBRAIC[CaMKb]+STATES[CaMKt];
ALGEBRAIC[vfrt] =  STATES[V]*CONSTANTS[frt];
ALGEBRAIC[U_1] =  CONSTANTS[B_1]*(STATES[V] - CONSTANTS[v0]);
ALGEBRAIC[A_1] = ( 4.00000*CONSTANTS[ffrt]*( STATES[cass]*exp( 2.00000*ALGEBRAIC[vfrt]) -  0.341000*CONSTANTS[cao]))/CONSTANTS[B_1];
ALGEBRAIC[fICaLp] = 1.00000/(1.00000+CONSTANTS[KmCaMK]/ALGEBRAIC[CaMKa]);
ALGEBRAIC[PhiCaL] = (- 1.00000e-07<=ALGEBRAIC[U_1]&&ALGEBRAIC[U_1]<=1.00000e-07 ?  ALGEBRAIC[A_1]*(1.00000 -  0.500000*ALGEBRAIC[U_1]) : ( ALGEBRAIC[A_1]*ALGEBRAIC[U_1])/(exp(ALGEBRAIC[U_1]) - 1.00000));
ALGEBRAIC[f] =  CONSTANTS[Aff]*STATES[ff]+ CONSTANTS[Afs]*STATES[fs];
ALGEBRAIC[tfcaf] = 7.00000+1.00000/( 0.0400000*exp(- (STATES[V] - 4.00000)/7.00000)+ 0.0400000*exp((STATES[V] - 4.00000)/7.00000));
ALGEBRAIC[tfcafp] =  2.50000*ALGEBRAIC[tfcaf];
ALGEBRAIC[Afcaf] = 0.300000+0.600000/(1.00000+exp((STATES[V] - 10.0000)/10.0000));
ALGEBRAIC[Afcas] = 1.00000 - ALGEBRAIC[Afcaf];
ALGEBRAIC[fca] =  ALGEBRAIC[Afcaf]*STATES[fcaf]+ ALGEBRAIC[Afcas]*STATES[fcas];
ALGEBRAIC[fp] =  CONSTANTS[Aff]*STATES[ffp]+ CONSTANTS[Afs]*STATES[fs];
ALGEBRAIC[fcap] =  ALGEBRAIC[Afcaf]*STATES[fcafp]+ ALGEBRAIC[Afcas]*STATES[fcas];
ALGEBRAIC[ICaL] =  (1.00000 - ALGEBRAIC[fICaLp])*CONSTANTS[PCa]*ALGEBRAIC[PhiCaL]*STATES[d]*( ALGEBRAIC[f]*(1.00000 - STATES[nca])+ STATES[jca]*ALGEBRAIC[fca]*STATES[nca])+ ALGEBRAIC[fICaLp]*CONSTANTS[PCap]*ALGEBRAIC[PhiCaL]*STATES[d]*( ALGEBRAIC[fp]*(1.00000 - STATES[nca])+ STATES[jca]*ALGEBRAIC[fcap]*STATES[nca]);
ALGEBRAIC[dss] = 1.00000/(1.00000+exp(- (STATES[V]+3.94000)/4.23000));
ALGEBRAIC[td] = 0.600000+1.00000/(exp( - 0.0500000*(STATES[V]+6.00000))+exp( 0.0900000*(STATES[V]+14.0000)));
ALGEBRAIC[km2n] =  STATES[jca]*1.00000;
ALGEBRAIC[anca] = 1.00000/(CONSTANTS[k2n]/ALGEBRAIC[km2n]+pow(1.00000+CONSTANTS[Kmn]/STATES[cass], 4.00000));
ALGEBRAIC[fss] = 1.00000/(1.00000+exp((STATES[V]+19.5800)/3.69600));
ALGEBRAIC[fcass] = ALGEBRAIC[fss];
ALGEBRAIC[Jleak] = ( 0.00393750*STATES[cansr])/15.0000;
ALGEBRAIC[Jupp] = ( CONSTANTS[upScale]*2.75000*0.00437500*STATES[cai])/((STATES[cai]+0.000920000) - 0.000170000);
ALGEBRAIC[Jupnp] = ( CONSTANTS[upScale]*0.00437500*STATES[cai])/(STATES[cai]+0.000920000);
ALGEBRAIC[fJupp] = 1.00000/(1.00000+CONSTANTS[KmCaMK]/ALGEBRAIC[CaMKa]);
ALGEBRAIC[Jup] =  CONSTANTS[Jup_b]*(( (1.00000 - ALGEBRAIC[fJupp])*ALGEBRAIC[Jupnp]+ ALGEBRAIC[fJupp]*ALGEBRAIC[Jupp]) - ALGEBRAIC[Jleak]);
ALGEBRAIC[tau_rel_temp] = CONSTANTS[bt]/(1.00000+0.0123000/STATES[cajsr]);
ALGEBRAIC[Jrel_inf_temp] = ( CONSTANTS[a_rel]*- ALGEBRAIC[ICaL])/(1.00000+ 1.00000*pow(1.50000/STATES[cajsr], 8.00000));
ALGEBRAIC[Jrel_inf] = (CONSTANTS[celltype]==2.00000 ?  ALGEBRAIC[Jrel_inf_temp]*1.70000 : ALGEBRAIC[Jrel_inf_temp]);
ALGEBRAIC[tau_rel] = (ALGEBRAIC[tau_rel_temp]<0.00100000 ? 0.00100000 : ALGEBRAIC[tau_rel_temp]);
ALGEBRAIC[Jrel_temp] = ( CONSTANTS[a_relp]*- ALGEBRAIC[ICaL])/(1.00000+pow(1.50000/STATES[cajsr], 8.00000));
ALGEBRAIC[tau_relp_temp] = CONSTANTS[btp]/(1.00000+0.0123000/STATES[cajsr]);
ALGEBRAIC[tau_relp] = (ALGEBRAIC[tau_relp_temp]<0.00100000 ? 0.00100000 : ALGEBRAIC[tau_relp_temp]);
ALGEBRAIC[Jrel_infp] = (CONSTANTS[celltype]==2.00000 ?  ALGEBRAIC[Jrel_temp]*1.70000 : ALGEBRAIC[Jrel_temp]);
ALGEBRAIC[Bcai] = 1.00000/(1.00000+( CONSTANTS[cmdnmax]*CONSTANTS[kmcmdn])/pow(CONSTANTS[kmcmdn]+STATES[cai], 2.00000)+( CONSTANTS[trpnmax]*CONSTANTS[kmtrpn])/pow(CONSTANTS[kmtrpn]+STATES[cai], 2.00000));
ALGEBRAIC[fJrelp] = 1.00000/(1.00000+CONSTANTS[KmCaMK]/ALGEBRAIC[CaMKa]);
ALGEBRAIC[Jrel] =  CONSTANTS[Jrel_scaling_factor]*( (1.00000 - ALGEBRAIC[fJrelp])*STATES[Jrelnp]+ ALGEBRAIC[fJrelp]*STATES[Jrelp]);
ALGEBRAIC[tff] = 7.00000+1.00000/( 0.00450000*exp(- (STATES[V]+20.0000)/10.0000)+ 0.00450000*exp((STATES[V]+20.0000)/10.0000));
ALGEBRAIC[tffp] =  2.50000*ALGEBRAIC[tff];
ALGEBRAIC[tfcas] = 100.000+1.00000/( 0.000120000*exp(- STATES[V]/3.00000)+ 0.000120000*exp(STATES[V]/7.00000));
ALGEBRAIC[tfs] = 1000.00+1.00000/( 3.50000e-05*exp(- (STATES[V]+5.00000)/4.00000)+ 3.50000e-05*exp((STATES[V]+5.00000)/6.00000));
ALGEBRAIC[Jdiff] = (STATES[cass] - STATES[cai])/0.200000;
ALGEBRAIC[Bcass] = 1.00000/(1.00000+( CONSTANTS[BSRmax]*CONSTANTS[KmBSR])/pow(CONSTANTS[KmBSR]+STATES[cass], 2.00000)+( CONSTANTS[BSLmax]*CONSTANTS[KmBSL])/pow(CONSTANTS[KmBSL]+STATES[cass], 2.00000));
ALGEBRAIC[Bcajsr] = 1.00000/(1.00000+( CONSTANTS[csqnmax]*CONSTANTS[kmcsqn])/pow(CONSTANTS[kmcsqn]+STATES[cajsr], 2.00000));
ALGEBRAIC[Jtr] = (STATES[cansr] - STATES[cajsr])/100.000;
#endif
#if defined CANA_ORUDY2011 
ALGEBRAIC[CaMKb] = ( CONSTANTS[CaMKo]*(1.00000 - STATES[CaMKt]))/(1.00000+CONSTANTS[KmCaM]/STATES[cass]);
ALGEBRAIC[CaMKa] = ALGEBRAIC[CaMKb]+STATES[CaMKt];
ALGEBRAIC[fICaLp] = 1.00000/(1.00000+CONSTANTS[KmCaMK]/ALGEBRAIC[CaMKa]);
ALGEBRAIC[vfrt] =  STATES[V]*CONSTANTS[frt];
ALGEBRAIC[U_2] =  CONSTANTS[B_2]*(STATES[V] - CONSTANTS[v0]);
ALGEBRAIC[A_2] = ( 0.750000*CONSTANTS[ffrt]*( STATES[nass]*exp(ALGEBRAIC[vfrt]) - CONSTANTS[nao]))/CONSTANTS[B_2];
ALGEBRAIC[PhiCaNa] = (- 1.00000e-07<=ALGEBRAIC[U_2]&&ALGEBRAIC[U_2]<=1.00000e-07 ?  ALGEBRAIC[A_2]*(1.00000 -  0.500000*ALGEBRAIC[U_2]) : ( ALGEBRAIC[A_2]*ALGEBRAIC[U_2])/(exp(ALGEBRAIC[U_2]) - 1.00000));
ALGEBRAIC[f] =  CONSTANTS[Aff]*STATES[ff]+ CONSTANTS[Afs]*STATES[fs];
ALGEBRAIC[fss] = 1.00000/(1.00000+exp((STATES[V]+19.5800)/3.69600));
ALGEBRAIC[fcass] = ALGEBRAIC[fss];
ALGEBRAIC[tfcaf] = 7.00000+1.00000/( 0.0400000*exp(- (STATES[V] - 4.00000)/7.00000)+ 0.0400000*exp((STATES[V] - 4.00000)/7.00000));
ALGEBRAIC[tfcafp] =  2.50000*ALGEBRAIC[tfcaf];
ALGEBRAIC[Afcaf] = 0.300000+0.600000/(1.00000+exp((STATES[V] - 10.0000)/10.0000));
ALGEBRAIC[Afcas] = 1.00000 - ALGEBRAIC[Afcaf];
ALGEBRAIC[fca] =  ALGEBRAIC[Afcaf]*STATES[fcaf]+ ALGEBRAIC[Afcas]*STATES[fcas];
ALGEBRAIC[fp] =  CONSTANTS[Aff]*STATES[ffp]+ CONSTANTS[Afs]*STATES[fs];
ALGEBRAIC[fcap] =  ALGEBRAIC[Afcaf]*STATES[fcafp]+ ALGEBRAIC[Afcas]*STATES[fcas];
ALGEBRAIC[ICaNa] =  (1.00000 - ALGEBRAIC[fICaLp])*CONSTANTS[PCaNa]*ALGEBRAIC[PhiCaNa]*STATES[d]*( ALGEBRAIC[f]*(1.00000 - STATES[nca])+ STATES[jca]*ALGEBRAIC[fca]*STATES[nca])+ ALGEBRAIC[fICaLp]*CONSTANTS[PCaNap]*ALGEBRAIC[PhiCaNa]*STATES[d]*( ALGEBRAIC[fp]*(1.00000 - STATES[nca])+ STATES[jca]*ALGEBRAIC[fcap]*STATES[nca]);
ALGEBRAIC[dss] = 1.00000/(1.00000+exp(- (STATES[V]+3.94000)/4.23000));
ALGEBRAIC[td] = 0.600000+1.00000/(exp( - 0.0500000*(STATES[V]+6.00000))+exp( 0.0900000*(STATES[V]+14.0000)));
ALGEBRAIC[km2n] =  STATES[jca]*1.00000;
ALGEBRAIC[anca] = 1.00000/(CONSTANTS[k2n]/ALGEBRAIC[km2n]+pow(1.00000+CONSTANTS[Kmn]/STATES[cass], 4.00000));
ALGEBRAIC[Jleak] = ( 0.00393750*STATES[cansr])/15.0000;
ALGEBRAIC[Jupp] = ( CONSTANTS[upScale]*2.75000*0.00437500*STATES[cai])/((STATES[cai]+0.000920000) - 0.000170000);
ALGEBRAIC[Jupnp] = ( CONSTANTS[upScale]*0.00437500*STATES[cai])/(STATES[cai]+0.000920000);
ALGEBRAIC[fJupp] = 1.00000/(1.00000+CONSTANTS[KmCaMK]/ALGEBRAIC[CaMKa]);
ALGEBRAIC[Jup] =  CONSTANTS[Jup_b]*(( (1.00000 - ALGEBRAIC[fJupp])*ALGEBRAIC[Jupnp]+ ALGEBRAIC[fJupp]*ALGEBRAIC[Jupp]) - ALGEBRAIC[Jleak]);
ALGEBRAIC[tau_rel_temp] = CONSTANTS[bt]/(1.00000+0.0123000/STATES[cajsr]);
ALGEBRAIC[Jrel_inf_temp] = ( CONSTANTS[a_rel]*- ALGEBRAIC[ICaL])/(1.00000+ 1.00000*pow(1.50000/STATES[cajsr], 8.00000));
ALGEBRAIC[Jrel_inf] = (CONSTANTS[celltype]==2.00000 ?  ALGEBRAIC[Jrel_inf_temp]*1.70000 : ALGEBRAIC[Jrel_inf_temp]);
ALGEBRAIC[tau_rel] = (ALGEBRAIC[tau_rel_temp]<0.00100000 ? 0.00100000 : ALGEBRAIC[tau_rel_temp]);
ALGEBRAIC[Jrel_temp] = ( CONSTANTS[a_relp]*- ALGEBRAIC[ICaL])/(1.00000+pow(1.50000/STATES[cajsr], 8.00000));
ALGEBRAIC[tau_relp_temp] = CONSTANTS[btp]/(1.00000+0.0123000/STATES[cajsr]);
ALGEBRAIC[tau_relp] = (ALGEBRAIC[tau_relp_temp]<0.00100000 ? 0.00100000 : ALGEBRAIC[tau_relp_temp]);
ALGEBRAIC[Jrel_infp] = (CONSTANTS[celltype]==2.00000 ?  ALGEBRAIC[Jrel_temp]*1.70000 : ALGEBRAIC[Jrel_temp]);
ALGEBRAIC[Bcai] = 1.00000/(1.00000+( CONSTANTS[cmdnmax]*CONSTANTS[kmcmdn])/pow(CONSTANTS[kmcmdn]+STATES[cai], 2.00000)+( CONSTANTS[trpnmax]*CONSTANTS[kmtrpn])/pow(CONSTANTS[kmtrpn]+STATES[cai], 2.00000));
ALGEBRAIC[fJrelp] = 1.00000/(1.00000+CONSTANTS[KmCaMK]/ALGEBRAIC[CaMKa]);
ALGEBRAIC[Jrel] =  CONSTANTS[Jrel_scaling_factor]*( (1.00000 - ALGEBRAIC[fJrelp])*STATES[Jrelnp]+ ALGEBRAIC[fJrelp]*STATES[Jrelp]);
ALGEBRAIC[tff] = 7.00000+1.00000/( 0.00450000*exp(- (STATES[V]+20.0000)/10.0000)+ 0.00450000*exp((STATES[V]+20.0000)/10.0000));
ALGEBRAIC[tffp] =  2.50000*ALGEBRAIC[tff];
ALGEBRAIC[tfcas] = 100.000+1.00000/( 0.000120000*exp(- STATES[V]/3.00000)+ 0.000120000*exp(STATES[V]/7.00000));
ALGEBRAIC[tfs] = 1000.00+1.00000/( 3.50000e-05*exp(- (STATES[V]+5.00000)/4.00000)+ 3.50000e-05*exp((STATES[V]+5.00000)/6.00000));
ALGEBRAIC[Jdiff] = (STATES[cass] - STATES[cai])/0.200000;
ALGEBRAIC[JdiffNa] = (STATES[nass] - STATES[nai])/2.00000;
ALGEBRAIC[Bcass] = 1.00000/(1.00000+( CONSTANTS[BSRmax]*CONSTANTS[KmBSR])/pow(CONSTANTS[KmBSR]+STATES[cass], 2.00000)+( CONSTANTS[BSLmax]*CONSTANTS[KmBSL])/pow(CONSTANTS[KmBSL]+STATES[cass], 2.00000));
ALGEBRAIC[Bcajsr] = 1.00000/(1.00000+( CONSTANTS[csqnmax]*CONSTANTS[kmcsqn])/pow(CONSTANTS[kmcsqn]+STATES[cajsr], 2.00000));
ALGEBRAIC[Jtr] = (STATES[cansr] - STATES[cajsr])/100.000;
#endif
#if defined K1_ORUDY2011 
ALGEBRAIC[EK] =  (( CONSTANTS[R]*CONSTANTS[T])/CONSTANTS[F])*log(CONSTANTS[ko]/STATES[ki]);
ALGEBRAIC[rk1] = 1.00000/(1.00000+exp(((STATES[V]+105.800) -  2.60000*CONSTANTS[ko])/9.49300));
ALGEBRAIC[IK1] =  CONSTANTS[GK1]* pow(CONSTANTS[ko], 1.0 / 2)*ALGEBRAIC[rk1]*STATES[xk1]*(STATES[V] - ALGEBRAIC[EK]);
ALGEBRAIC[xk1ss] = 1.00000/(1.00000+exp(- (STATES[V]+ 2.55380*CONSTANTS[ko]+144.590)/( 1.56920*CONSTANTS[ko]+3.81150)));
ALGEBRAIC[txk1] = 122.200/(exp(- (STATES[V]+127.200)/20.3600)+exp((STATES[V]+236.800)/69.3300));
ALGEBRAIC[JdiffK] = (STATES[kss] - STATES[ki])/2.00000;
#endif
#if defined KB_ORUDY2011 
ALGEBRAIC[xkb] = 1.00000/(1.00000+exp(- (STATES[V] - 14.4800)/18.3400));
ALGEBRAIC[EK] =  (( CONSTANTS[R]*CONSTANTS[T])/CONSTANTS[F])*log(CONSTANTS[ko]/STATES[ki]);
ALGEBRAIC[IKb] =  CONSTANTS[GKb]*ALGEBRAIC[xkb]*(STATES[V] - ALGEBRAIC[EK]);
ALGEBRAIC[JdiffK] = (STATES[kss] - STATES[ki])/2.00000;
#endif
#if defined KR_ORUDY2011 
ALGEBRAIC[EK] =  (( CONSTANTS[R]*CONSTANTS[T])/CONSTANTS[F])*log(CONSTANTS[ko]/STATES[ki]);
ALGEBRAIC[IKr] =  CONSTANTS[GKr]* pow((CONSTANTS[ko]/5.40000), 1.0 / 2)*STATES[O]*(STATES[V] - ALGEBRAIC[EK]);
ALGEBRAIC[JdiffK] = (STATES[kss] - STATES[ki])/2.00000;
#endif
#if defined KS_ORUDY2011 
ALGEBRAIC[CaMKb] = ( CONSTANTS[CaMKo]*(1.00000 - STATES[CaMKt]))/(1.00000+CONSTANTS[KmCaM]/STATES[cass]);
ALGEBRAIC[CaMKa] = ALGEBRAIC[CaMKb]+STATES[CaMKt];
ALGEBRAIC[xs1ss] = 1.00000/(1.00000+exp(- (STATES[V]+11.6000)/8.93200));
ALGEBRAIC[xs2ss] = ALGEBRAIC[xs1ss];
ALGEBRAIC[txs1] = CONSTANTS[txs1_max]+1.00000/( 0.000232600*exp((STATES[V]+48.2800)/17.8000)+ 0.00129200*exp(- (STATES[V]+210.000)/230.000));
ALGEBRAIC[txs2] = 1.00000/( 0.0100000*exp((STATES[V] - 50.0000)/20.0000)+ 0.0193000*exp(- (STATES[V]+66.5400)/31.0000));
ALGEBRAIC[EKs] =  (( CONSTANTS[R]*CONSTANTS[T])/CONSTANTS[F])*log((CONSTANTS[ko]+ CONSTANTS[PKNa]*CONSTANTS[nao])/(STATES[ki]+ CONSTANTS[PKNa]*STATES[nai]));
ALGEBRAIC[KsCa] = 1.00000+0.600000/(1.00000+pow(3.80000e-05/STATES[cai], 1.40000));
ALGEBRAIC[IKs] =  CONSTANTS[GKs]*ALGEBRAIC[KsCa]*STATES[xs1]*STATES[xs2]*(STATES[V] - ALGEBRAIC[EKs]);
#endif
#if defined NA_ORUDY2011 
ALGEBRAIC[CaMKb] = ( CONSTANTS[CaMKo]*(1.00000 - STATES[CaMKt]))/(1.00000+CONSTANTS[KmCaM]/STATES[cass]);
ALGEBRAIC[CaMKa] = ALGEBRAIC[CaMKb]+STATES[CaMKt];
ALGEBRAIC[ENa] =  (( CONSTANTS[R]*CONSTANTS[T])/CONSTANTS[F])*log(CONSTANTS[nao]/STATES[nai]);
ALGEBRAIC[fINap] = 1.00000/(1.00000+CONSTANTS[KmCaMK]/ALGEBRAIC[CaMKa]);
ALGEBRAIC[h] =  CONSTANTS[Ahf]*STATES[hf]+ CONSTANTS[Ahs]*STATES[hs];
ALGEBRAIC[hp] =  CONSTANTS[Ahf]*STATES[hf]+ CONSTANTS[Ahs]*STATES[hsp];
ALGEBRAIC[mss] = 1.00000/(1.00000+exp(- (STATES[V]+CONSTANTS[mssV1])/CONSTANTS[mssV2]));
ALGEBRAIC[tm] = 1.00000/( CONSTANTS[mtD1]*exp((STATES[V]+CONSTANTS[mtV1])/CONSTANTS[mtV2])+ CONSTANTS[mtD2]*exp(- (STATES[V]+CONSTANTS[mtV3])/CONSTANTS[mtV4]));
ALGEBRAIC[ths] = 1.00000/( 0.00979400*exp(- ((STATES[V]+17.9500) - CONSTANTS[shift_INa_inact])/28.0500)+ 0.334300*exp(((STATES[V]+5.73000) - CONSTANTS[shift_INa_inact])/56.6600));
ALGEBRAIC[thsp] =  3.00000*ALGEBRAIC[ths];
ALGEBRAIC[thf] = 1.00000/( 1.43200e-05*exp(- ((STATES[V]+1.19600) - CONSTANTS[shift_INa_inact])/6.28500)+ 6.14900*exp(((STATES[V]+0.509600) - CONSTANTS[shift_INa_inact])/20.2700));
ALGEBRAIC[tj] = 2.03800+1.00000/( 0.0213600*exp(- ((STATES[V]+100.600) - CONSTANTS[shift_INa_inact])/8.28100)+ 0.305200*exp(((STATES[V]+0.994100) - CONSTANTS[shift_INa_inact])/38.4500));
ALGEBRAIC[hss] = 1.00000/(1.00000+exp(((STATES[V]+CONSTANTS[hssV1]) - CONSTANTS[shift_INa_inact])/CONSTANTS[hssV2]));
ALGEBRAIC[jss] = ALGEBRAIC[hss];
ALGEBRAIC[JdiffNa] = (STATES[nass] - STATES[nai])/2.00000;
ALGEBRAIC[tjp] =  1.46000*ALGEBRAIC[tj];
ALGEBRAIC[hssp] = 1.00000/(1.00000+exp(((STATES[V]+89.1000) - CONSTANTS[shift_INa_inact])/6.08600));
ALGEBRAIC[INa] =  CONSTANTS[GNa]*(STATES[V] - ALGEBRAIC[ENa])*pow(STATES[m], 3.00000)*( (1.00000 - ALGEBRAIC[fINap])*ALGEBRAIC[h]*STATES[j]+ ALGEBRAIC[fINap]*ALGEBRAIC[hp]*STATES[jp]);
#endif
#if defined NAB_ORUDY2011 
ALGEBRAIC[A_INab] = ( CONSTANTS[PNab]*CONSTANTS[ffrt]*( STATES[nai]*exp(ALGEBRAIC[vfrt]) - CONSTANTS[nao]))/CONSTANTS[B_INab];
ALGEBRAIC[vfrt] =  STATES[V]*CONSTANTS[frt];
ALGEBRAIC[U_INab] =  CONSTANTS[B_INab]*(STATES[V] - CONSTANTS[v0_INab]);
ALGEBRAIC[INab] = (- 1.00000e-07<=ALGEBRAIC[U_INab]&&ALGEBRAIC[U_INab]<=1.00000e-07 ?  ALGEBRAIC[A_INab]*(1.00000 -  0.500000*ALGEBRAIC[U_INab]) : ( ALGEBRAIC[A_INab]*ALGEBRAIC[U_INab])/(exp(ALGEBRAIC[U_INab]) - 1.00000));
ALGEBRAIC[JdiffNa] = (STATES[nass] - STATES[nai])/2.00000;
#endif
#if defined NACA_I_ORUDY2011 
ALGEBRAIC[CaMKb] = ( CONSTANTS[CaMKo]*(1.00000 - STATES[CaMKt]))/(1.00000+CONSTANTS[KmCaM]/STATES[cass]);
ALGEBRAIC[CaMKa] = ALGEBRAIC[CaMKb]+STATES[CaMKt];
ALGEBRAIC[hna] = exp(( CONSTANTS[qna]*STATES[V]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[T]));
ALGEBRAIC[hca] = exp(( CONSTANTS[qca]*STATES[V]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[T]));
ALGEBRAIC[h1_i] = 1.00000+ (STATES[nai]/CONSTANTS[kna3])*(1.00000+ALGEBRAIC[hna]);
ALGEBRAIC[h2_i] = ( STATES[nai]*ALGEBRAIC[hna])/( CONSTANTS[kna3]*ALGEBRAIC[h1_i]);
ALGEBRAIC[h3_i] = 1.00000/ALGEBRAIC[h1_i];
ALGEBRAIC[h4_i] = 1.00000+ (STATES[nai]/CONSTANTS[kna1])*(1.00000+STATES[nai]/CONSTANTS[kna2]);
ALGEBRAIC[h5_i] = ( STATES[nai]*STATES[nai])/( ALGEBRAIC[h4_i]*CONSTANTS[kna1]*CONSTANTS[kna2]);
ALGEBRAIC[h6_i] = 1.00000/ALGEBRAIC[h4_i];
ALGEBRAIC[h7_i] = 1.00000+ (CONSTANTS[nao]/CONSTANTS[kna3])*(1.00000+1.00000/ALGEBRAIC[hna]);
ALGEBRAIC[h8_i] = CONSTANTS[nao]/( CONSTANTS[kna3]*ALGEBRAIC[hna]*ALGEBRAIC[h7_i]);
ALGEBRAIC[h9_i] = 1.00000/ALGEBRAIC[h7_i];
ALGEBRAIC[k3p_i] =  ALGEBRAIC[h9_i]*CONSTANTS[wca];
ALGEBRAIC[k3pp_i] =  ALGEBRAIC[h8_i]*CONSTANTS[wnaca];
ALGEBRAIC[k3_i] = ALGEBRAIC[k3p_i]+ALGEBRAIC[k3pp_i];
ALGEBRAIC[k4p_i] = ( ALGEBRAIC[h3_i]*CONSTANTS[wca])/ALGEBRAIC[hca];
ALGEBRAIC[k4pp_i] =  ALGEBRAIC[h2_i]*CONSTANTS[wnaca];
ALGEBRAIC[k4_i] = ALGEBRAIC[k4p_i]+ALGEBRAIC[k4pp_i];
ALGEBRAIC[k6_i] =  ALGEBRAIC[h6_i]*STATES[cai]*CONSTANTS[kcaon];
ALGEBRAIC[k7_i] =  ALGEBRAIC[h5_i]*ALGEBRAIC[h2_i]*CONSTANTS[wna];
ALGEBRAIC[k8_i] =  ALGEBRAIC[h8_i]*CONSTANTS[h11_i]*CONSTANTS[wna];
ALGEBRAIC[x1_i] =  CONSTANTS[k2_i]*ALGEBRAIC[k4_i]*(ALGEBRAIC[k7_i]+ALGEBRAIC[k6_i])+ CONSTANTS[k5_i]*ALGEBRAIC[k7_i]*(CONSTANTS[k2_i]+ALGEBRAIC[k3_i]);
ALGEBRAIC[x2_i] =  CONSTANTS[k1_i]*ALGEBRAIC[k7_i]*(ALGEBRAIC[k4_i]+CONSTANTS[k5_i])+ ALGEBRAIC[k4_i]*ALGEBRAIC[k6_i]*(CONSTANTS[k1_i]+ALGEBRAIC[k8_i]);
ALGEBRAIC[x3_i] =  CONSTANTS[k1_i]*ALGEBRAIC[k3_i]*(ALGEBRAIC[k7_i]+ALGEBRAIC[k6_i])+ ALGEBRAIC[k8_i]*ALGEBRAIC[k6_i]*(CONSTANTS[k2_i]+ALGEBRAIC[k3_i]);
ALGEBRAIC[x4_i] =  CONSTANTS[k2_i]*ALGEBRAIC[k8_i]*(ALGEBRAIC[k4_i]+CONSTANTS[k5_i])+ ALGEBRAIC[k3_i]*CONSTANTS[k5_i]*(CONSTANTS[k1_i]+ALGEBRAIC[k8_i]);
ALGEBRAIC[E1_i] = ALGEBRAIC[x1_i]/(ALGEBRAIC[x1_i]+ALGEBRAIC[x2_i]+ALGEBRAIC[x3_i]+ALGEBRAIC[x4_i]);
ALGEBRAIC[E2_i] = ALGEBRAIC[x2_i]/(ALGEBRAIC[x1_i]+ALGEBRAIC[x2_i]+ALGEBRAIC[x3_i]+ALGEBRAIC[x4_i]);
ALGEBRAIC[E3_i] = ALGEBRAIC[x3_i]/(ALGEBRAIC[x1_i]+ALGEBRAIC[x2_i]+ALGEBRAIC[x3_i]+ALGEBRAIC[x4_i]);
ALGEBRAIC[E4_i] = ALGEBRAIC[x4_i]/(ALGEBRAIC[x1_i]+ALGEBRAIC[x2_i]+ALGEBRAIC[x3_i]+ALGEBRAIC[x4_i]);
ALGEBRAIC[JncxCa_i] =  ALGEBRAIC[E2_i]*CONSTANTS[k2_i] -  ALGEBRAIC[E1_i]*CONSTANTS[k1_i];
ALGEBRAIC[JncxNa_i] = ( 3.00000*( ALGEBRAIC[E4_i]*ALGEBRAIC[k7_i] -  ALGEBRAIC[E1_i]*ALGEBRAIC[k8_i])+ ALGEBRAIC[E3_i]*ALGEBRAIC[k4pp_i]) -  ALGEBRAIC[E2_i]*ALGEBRAIC[k3pp_i];
ALGEBRAIC[allo_i] = 1.00000/(1.00000+pow(CONSTANTS[KmCaAct]/STATES[cai], 2.00000));
ALGEBRAIC[INaCa_i] =  0.800000*CONSTANTS[Gncx]*ALGEBRAIC[allo_i]*( CONSTANTS[zna]*ALGEBRAIC[JncxNa_i]+ CONSTANTS[zca]*ALGEBRAIC[JncxCa_i]);
ALGEBRAIC[Jleak] = ( 0.00393750*STATES[cansr])/15.0000;
ALGEBRAIC[Jupp] = ( CONSTANTS[upScale]*2.75000*0.00437500*STATES[cai])/((STATES[cai]+0.000920000) - 0.000170000);
ALGEBRAIC[Jupnp] = ( CONSTANTS[upScale]*0.00437500*STATES[cai])/(STATES[cai]+0.000920000);
ALGEBRAIC[fJupp] = 1.00000/(1.00000+CONSTANTS[KmCaMK]/ALGEBRAIC[CaMKa]);
ALGEBRAIC[Jup] =  CONSTANTS[Jup_b]*(( (1.00000 - ALGEBRAIC[fJupp])*ALGEBRAIC[Jupnp]+ ALGEBRAIC[fJupp]*ALGEBRAIC[Jupp]) - ALGEBRAIC[Jleak]);
ALGEBRAIC[tau_rel_temp] = CONSTANTS[bt]/(1.00000+0.0123000/STATES[cajsr]);
ALGEBRAIC[Jrel_inf_temp] = ( CONSTANTS[a_rel]*- ALGEBRAIC[ICaL])/(1.00000+ 1.00000*pow(1.50000/STATES[cajsr], 8.00000));
ALGEBRAIC[Jrel_inf] = (CONSTANTS[celltype]==2.00000 ?  ALGEBRAIC[Jrel_inf_temp]*1.70000 : ALGEBRAIC[Jrel_inf_temp]);
ALGEBRAIC[tau_rel] = (ALGEBRAIC[tau_rel_temp]<0.00100000 ? 0.00100000 : ALGEBRAIC[tau_rel_temp]);
ALGEBRAIC[Jrel_temp] = ( CONSTANTS[a_relp]*- ALGEBRAIC[ICaL])/(1.00000+pow(1.50000/STATES[cajsr], 8.00000));
ALGEBRAIC[tau_relp_temp] = CONSTANTS[btp]/(1.00000+0.0123000/STATES[cajsr]);
ALGEBRAIC[tau_relp] = (ALGEBRAIC[tau_relp_temp]<0.00100000 ? 0.00100000 : ALGEBRAIC[tau_relp_temp]);
ALGEBRAIC[Jrel_infp] = (CONSTANTS[celltype]==2.00000 ?  ALGEBRAIC[Jrel_temp]*1.70000 : ALGEBRAIC[Jrel_temp]);
ALGEBRAIC[Bcai] = 1.00000/(1.00000+( CONSTANTS[cmdnmax]*CONSTANTS[kmcmdn])/pow(CONSTANTS[kmcmdn]+STATES[cai], 2.00000)+( CONSTANTS[trpnmax]*CONSTANTS[kmtrpn])/pow(CONSTANTS[kmtrpn]+STATES[cai], 2.00000));
ALGEBRAIC[fJrelp] = 1.00000/(1.00000+CONSTANTS[KmCaMK]/ALGEBRAIC[CaMKa]);
ALGEBRAIC[Jrel] =  CONSTANTS[Jrel_scaling_factor]*( (1.00000 - ALGEBRAIC[fJrelp])*STATES[Jrelnp]+ ALGEBRAIC[fJrelp]*STATES[Jrelp]);
ALGEBRAIC[Jtr] = (STATES[cansr] - STATES[cajsr])/100.000;
ALGEBRAIC[Bcajsr] = 1.00000/(1.00000+( CONSTANTS[csqnmax]*CONSTANTS[kmcsqn])/pow(CONSTANTS[kmcsqn]+STATES[cajsr], 2.00000));
ALGEBRAIC[JdiffNa] = (STATES[nass] - STATES[nai])/2.00000;
ALGEBRAIC[Jdiff] = (STATES[cass] - STATES[cai])/0.200000;
ALGEBRAIC[Bcass] = 1.00000/(1.00000+( CONSTANTS[BSRmax]*CONSTANTS[KmBSR])/pow(CONSTANTS[KmBSR]+STATES[cass], 2.00000)+( CONSTANTS[BSLmax]*CONSTANTS[KmBSL])/pow(CONSTANTS[KmBSL]+STATES[cass], 2.00000));
#endif
#if defined NACA_SS_ORUDY2011 
ALGEBRAIC[CaMKb] = ( CONSTANTS[CaMKo]*(1.00000 - STATES[CaMKt]))/(1.00000+CONSTANTS[KmCaM]/STATES[cass]);
ALGEBRAIC[CaMKa] = ALGEBRAIC[CaMKb]+STATES[CaMKt];
ALGEBRAIC[hna] = exp(( CONSTANTS[qna]*STATES[V]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[T]));
ALGEBRAIC[hca] = exp(( CONSTANTS[qca]*STATES[V]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[T]));
ALGEBRAIC[h1_ss] = 1.00000+ (STATES[nass]/CONSTANTS[kna3])*(1.00000+ALGEBRAIC[hna]);
ALGEBRAIC[h2_ss] = ( STATES[nass]*ALGEBRAIC[hna])/( CONSTANTS[kna3]*ALGEBRAIC[h1_ss]);
ALGEBRAIC[h3_ss] = 1.00000/ALGEBRAIC[h1_ss];
ALGEBRAIC[h4_ss] = 1.00000+ (STATES[nass]/CONSTANTS[kna1])*(1.00000+STATES[nass]/CONSTANTS[kna2]);
ALGEBRAIC[h5_ss] = ( STATES[nass]*STATES[nass])/( ALGEBRAIC[h4_ss]*CONSTANTS[kna1]*CONSTANTS[kna2]);
ALGEBRAIC[h6_ss] = 1.00000/ALGEBRAIC[h4_ss];
ALGEBRAIC[h7_ss] = 1.00000+ (CONSTANTS[nao]/CONSTANTS[kna3])*(1.00000+1.00000/ALGEBRAIC[hna]);
ALGEBRAIC[h8_ss] = CONSTANTS[nao]/( CONSTANTS[kna3]*ALGEBRAIC[hna]*ALGEBRAIC[h7_ss]);
ALGEBRAIC[h9_ss] = 1.00000/ALGEBRAIC[h7_ss];
ALGEBRAIC[k3p_ss] =  ALGEBRAIC[h9_ss]*CONSTANTS[wca];
ALGEBRAIC[k3pp_ss] =  ALGEBRAIC[h8_ss]*CONSTANTS[wnaca];
ALGEBRAIC[k3_ss] = ALGEBRAIC[k3p_ss]+ALGEBRAIC[k3pp_ss];
ALGEBRAIC[k4p_ss] = ( ALGEBRAIC[h3_ss]*CONSTANTS[wca])/ALGEBRAIC[hca];
ALGEBRAIC[k4pp_ss] =  ALGEBRAIC[h2_ss]*CONSTANTS[wnaca];
ALGEBRAIC[k4_ss] = ALGEBRAIC[k4p_ss]+ALGEBRAIC[k4pp_ss];
ALGEBRAIC[k6_ss] =  ALGEBRAIC[h6_ss]*STATES[cass]*CONSTANTS[kcaon];
ALGEBRAIC[k7_ss] =  ALGEBRAIC[h5_ss]*ALGEBRAIC[h2_ss]*CONSTANTS[wna];
ALGEBRAIC[k8_ss] =  ALGEBRAIC[h8_ss]*CONSTANTS[h11_ss]*CONSTANTS[wna];
ALGEBRAIC[x1_ss] =  CONSTANTS[k2_ss]*ALGEBRAIC[k4_ss]*(ALGEBRAIC[k7_ss]+ALGEBRAIC[k6_ss])+ CONSTANTS[k5_ss]*ALGEBRAIC[k7_ss]*(CONSTANTS[k2_ss]+ALGEBRAIC[k3_ss]);
ALGEBRAIC[x2_ss] =  CONSTANTS[k1_ss]*ALGEBRAIC[k7_ss]*(ALGEBRAIC[k4_ss]+CONSTANTS[k5_ss])+ ALGEBRAIC[k4_ss]*ALGEBRAIC[k6_ss]*(CONSTANTS[k1_ss]+ALGEBRAIC[k8_ss]);
ALGEBRAIC[x3_ss] =  CONSTANTS[k1_ss]*ALGEBRAIC[k3_ss]*(ALGEBRAIC[k7_ss]+ALGEBRAIC[k6_ss])+ ALGEBRAIC[k8_ss]*ALGEBRAIC[k6_ss]*(CONSTANTS[k2_ss]+ALGEBRAIC[k3_ss]);
ALGEBRAIC[x4_ss] =  CONSTANTS[k2_ss]*ALGEBRAIC[k8_ss]*(ALGEBRAIC[k4_ss]+CONSTANTS[k5_ss])+ ALGEBRAIC[k3_ss]*CONSTANTS[k5_ss]*(CONSTANTS[k1_ss]+ALGEBRAIC[k8_ss]);
ALGEBRAIC[E1_ss] = ALGEBRAIC[x1_ss]/(ALGEBRAIC[x1_ss]+ALGEBRAIC[x2_ss]+ALGEBRAIC[x3_ss]+ALGEBRAIC[x4_ss]);
ALGEBRAIC[E2_ss] = ALGEBRAIC[x2_ss]/(ALGEBRAIC[x1_ss]+ALGEBRAIC[x2_ss]+ALGEBRAIC[x3_ss]+ALGEBRAIC[x4_ss]);
ALGEBRAIC[E3_ss] = ALGEBRAIC[x3_ss]/(ALGEBRAIC[x1_ss]+ALGEBRAIC[x2_ss]+ALGEBRAIC[x3_ss]+ALGEBRAIC[x4_ss]);
ALGEBRAIC[E4_ss] = ALGEBRAIC[x4_ss]/(ALGEBRAIC[x1_ss]+ALGEBRAIC[x2_ss]+ALGEBRAIC[x3_ss]+ALGEBRAIC[x4_ss]);
ALGEBRAIC[JncxCa_ss] =  ALGEBRAIC[E2_ss]*CONSTANTS[k2_ss] -  ALGEBRAIC[E1_ss]*CONSTANTS[k1_ss];
ALGEBRAIC[JncxNa_ss] = ( 3.00000*( ALGEBRAIC[E4_ss]*ALGEBRAIC[k7_ss] -  ALGEBRAIC[E1_ss]*ALGEBRAIC[k8_ss])+ ALGEBRAIC[E3_ss]*ALGEBRAIC[k4pp_ss]) -  ALGEBRAIC[E2_ss]*ALGEBRAIC[k3pp_ss];
ALGEBRAIC[allo_ss] = 1.00000/(1.00000+pow(CONSTANTS[KmCaAct]/STATES[cass], 2.00000));
ALGEBRAIC[INaCa_ss] =  0.200000*CONSTANTS[Gncx]*ALGEBRAIC[allo_ss]*( CONSTANTS[zna]*ALGEBRAIC[JncxNa_ss]+ CONSTANTS[zca]*ALGEBRAIC[JncxCa_ss]);
ALGEBRAIC[Jleak] = ( 0.00393750*STATES[cansr])/15.0000;
ALGEBRAIC[Jupp] = ( CONSTANTS[upScale]*2.75000*0.00437500*STATES[cai])/((STATES[cai]+0.000920000) - 0.000170000);
ALGEBRAIC[Jupnp] = ( CONSTANTS[upScale]*0.00437500*STATES[cai])/(STATES[cai]+0.000920000);
ALGEBRAIC[fJupp] = 1.00000/(1.00000+CONSTANTS[KmCaMK]/ALGEBRAIC[CaMKa]);
ALGEBRAIC[Jup] =  CONSTANTS[Jup_b]*(( (1.00000 - ALGEBRAIC[fJupp])*ALGEBRAIC[Jupnp]+ ALGEBRAIC[fJupp]*ALGEBRAIC[Jupp]) - ALGEBRAIC[Jleak]);
ALGEBRAIC[tau_rel_temp] = CONSTANTS[bt]/(1.00000+0.0123000/STATES[cajsr]);
ALGEBRAIC[Jrel_inf_temp] = ( CONSTANTS[a_rel]*- ALGEBRAIC[ICaL])/(1.00000+ 1.00000*pow(1.50000/STATES[cajsr], 8.00000));
ALGEBRAIC[Jrel_inf] = (CONSTANTS[celltype]==2.00000 ?  ALGEBRAIC[Jrel_inf_temp]*1.70000 : ALGEBRAIC[Jrel_inf_temp]);
ALGEBRAIC[tau_rel] = (ALGEBRAIC[tau_rel_temp]<0.00100000 ? 0.00100000 : ALGEBRAIC[tau_rel_temp]);
ALGEBRAIC[Jrel_temp] = ( CONSTANTS[a_relp]*- ALGEBRAIC[ICaL])/(1.00000+pow(1.50000/STATES[cajsr], 8.00000));
ALGEBRAIC[tau_relp_temp] = CONSTANTS[btp]/(1.00000+0.0123000/STATES[cajsr]);
ALGEBRAIC[tau_relp] = (ALGEBRAIC[tau_relp_temp]<0.00100000 ? 0.00100000 : ALGEBRAIC[tau_relp_temp]);
ALGEBRAIC[Jrel_infp] = (CONSTANTS[celltype]==2.00000 ?  ALGEBRAIC[Jrel_temp]*1.70000 : ALGEBRAIC[Jrel_temp]);
ALGEBRAIC[Bcai] = 1.00000/(1.00000+( CONSTANTS[cmdnmax]*CONSTANTS[kmcmdn])/pow(CONSTANTS[kmcmdn]+STATES[cai], 2.00000)+( CONSTANTS[trpnmax]*CONSTANTS[kmtrpn])/pow(CONSTANTS[kmtrpn]+STATES[cai], 2.00000));
ALGEBRAIC[fJrelp] = 1.00000/(1.00000+CONSTANTS[KmCaMK]/ALGEBRAIC[CaMKa]);
ALGEBRAIC[Jrel] =  CONSTANTS[Jrel_scaling_factor]*( (1.00000 - ALGEBRAIC[fJrelp])*STATES[Jrelnp]+ ALGEBRAIC[fJrelp]*STATES[Jrelp]);
ALGEBRAIC[Jtr] = (STATES[cansr] - STATES[cajsr])/100.000;
ALGEBRAIC[Bcajsr] = 1.00000/(1.00000+( CONSTANTS[csqnmax]*CONSTANTS[kmcsqn])/pow(CONSTANTS[kmcsqn]+STATES[cajsr], 2.00000));
ALGEBRAIC[JdiffNa] = (STATES[nass] - STATES[nai])/2.00000;
ALGEBRAIC[Jdiff] = (STATES[cass] - STATES[cai])/0.200000;
ALGEBRAIC[Bcass] = 1.00000/(1.00000+( CONSTANTS[BSRmax]*CONSTANTS[KmBSR])/pow(CONSTANTS[KmBSR]+STATES[cass], 2.00000)+( CONSTANTS[BSLmax]*CONSTANTS[KmBSL])/pow(CONSTANTS[KmBSL]+STATES[cass], 2.00000));
#endif
#if defined NAK_ORUDY2011 
ALGEBRAIC[Knao] =  CONSTANTS[Knao0]*exp(( (1.00000 - CONSTANTS[delta])*STATES[V]*CONSTANTS[F])/( 3.00000*CONSTANTS[R]*CONSTANTS[T]));
ALGEBRAIC[Knai] =  CONSTANTS[Knai0]*exp(( CONSTANTS[delta]*STATES[V]*CONSTANTS[F])/( 3.00000*CONSTANTS[R]*CONSTANTS[T]));
ALGEBRAIC[b2] = ( CONSTANTS[k2m]*pow(CONSTANTS[nao]/ALGEBRAIC[Knao], 3.00000))/((pow(1.00000+CONSTANTS[nao]/ALGEBRAIC[Knao], 3.00000)+pow(1.00000+CONSTANTS[ko]/CONSTANTS[Kko], 2.00000)) - 1.00000);
ALGEBRAIC[b4] = ( CONSTANTS[k4m]*pow(STATES[ki]/CONSTANTS[Kki], 2.00000))/((pow(1.00000+STATES[nai]/ALGEBRAIC[Knai], 3.00000)+pow(1.00000+STATES[ki]/CONSTANTS[Kki], 2.00000)) - 1.00000);
ALGEBRAIC[a1] = ( CONSTANTS[k1p]*pow(STATES[nai]/ALGEBRAIC[Knai], 3.00000))/((pow(1.00000+STATES[nai]/ALGEBRAIC[Knai], 3.00000)+pow(1.00000+STATES[ki]/CONSTANTS[Kki], 2.00000)) - 1.00000);
ALGEBRAIC[a3] = ( CONSTANTS[k3p]*pow(CONSTANTS[ko]/CONSTANTS[Kko], 2.00000))/((pow(1.00000+CONSTANTS[nao]/ALGEBRAIC[Knao], 3.00000)+pow(1.00000+CONSTANTS[ko]/CONSTANTS[Kko], 2.00000)) - 1.00000);
ALGEBRAIC[b3] = ( CONSTANTS[k3m]*ALGEBRAIC[P]*CONSTANTS[H])/(1.00000+CONSTANTS[MgATP]/CONSTANTS[Kmgatp]);
ALGEBRAIC[x1] =  CONSTANTS[a4]*ALGEBRAIC[a1]*CONSTANTS[a2]+ ALGEBRAIC[b2]*ALGEBRAIC[b4]*ALGEBRAIC[b3]+ CONSTANTS[a2]*ALGEBRAIC[b4]*ALGEBRAIC[b3]+ ALGEBRAIC[b3]*ALGEBRAIC[a1]*CONSTANTS[a2];
ALGEBRAIC[x2] =  ALGEBRAIC[b2]*CONSTANTS[b1]*ALGEBRAIC[b4]+ ALGEBRAIC[a1]*CONSTANTS[a2]*ALGEBRAIC[a3]+ ALGEBRAIC[a3]*CONSTANTS[b1]*ALGEBRAIC[b4]+ CONSTANTS[a2]*ALGEBRAIC[a3]*ALGEBRAIC[b4];
ALGEBRAIC[x3] =  CONSTANTS[a2]*ALGEBRAIC[a3]*CONSTANTS[a4]+ ALGEBRAIC[b3]*ALGEBRAIC[b2]*CONSTANTS[b1]+ ALGEBRAIC[b2]*CONSTANTS[b1]*CONSTANTS[a4]+ ALGEBRAIC[a3]*CONSTANTS[a4]*CONSTANTS[b1];
ALGEBRAIC[x4] =  ALGEBRAIC[b4]*ALGEBRAIC[b3]*ALGEBRAIC[b2]+ ALGEBRAIC[a3]*CONSTANTS[a4]*ALGEBRAIC[a1]+ ALGEBRAIC[b2]*CONSTANTS[a4]*ALGEBRAIC[a1]+ ALGEBRAIC[b3]*ALGEBRAIC[b2]*ALGEBRAIC[a1];
ALGEBRAIC[P] = CONSTANTS[eP]/(1.00000+CONSTANTS[H]/CONSTANTS[Khp]+STATES[nai]/CONSTANTS[Knap]+STATES[ki]/CONSTANTS[Kxkur]);
ALGEBRAIC[E1] = ALGEBRAIC[x1]/(ALGEBRAIC[x1]+ALGEBRAIC[x2]+ALGEBRAIC[x3]+ALGEBRAIC[x4]);
ALGEBRAIC[a3] = ( CONSTANTS[k3p]*pow(CONSTANTS[ko]/CONSTANTS[Kko], 2.00000))/((pow(1.00000+CONSTANTS[nao]/ALGEBRAIC[Knao], 3.00000)+pow(1.00000+CONSTANTS[ko]/CONSTANTS[Kko], 2.00000)) - 1.00000);
ALGEBRAIC[E2] = ALGEBRAIC[x2]/(ALGEBRAIC[x1]+ALGEBRAIC[x2]+ALGEBRAIC[x3]+ALGEBRAIC[x4]);
ALGEBRAIC[E4] = ALGEBRAIC[x4]/(ALGEBRAIC[x1]+ALGEBRAIC[x2]+ALGEBRAIC[x3]+ALGEBRAIC[x4]);
ALGEBRAIC[E3] = ALGEBRAIC[x3]/(ALGEBRAIC[x1]+ALGEBRAIC[x2]+ALGEBRAIC[x3]+ALGEBRAIC[x4]);
ALGEBRAIC[JnakNa] =  3.00000*( ALGEBRAIC[E1]*ALGEBRAIC[a3] -  ALGEBRAIC[E2]*ALGEBRAIC[b3]);
ALGEBRAIC[JnakK] =  2.00000*( ALGEBRAIC[E4]*CONSTANTS[b1] -  ALGEBRAIC[E3]*ALGEBRAIC[a1]);
ALGEBRAIC[INaK] =  CONSTANTS[Pnak]*( CONSTANTS[zna]*ALGEBRAIC[JnakNa]+ CONSTANTS[zk]*ALGEBRAIC[JnakK]);
ALGEBRAIC[JdiffK] = (STATES[kss] - STATES[ki])/2.00000;
ALGEBRAIC[JdiffNa] = (STATES[nass] - STATES[nai])/2.00000;
#endif
#if defined NAL_ORUDY2011 
ALGEBRAIC[CaMKb] = ( CONSTANTS[CaMKo]*(1.00000 - STATES[CaMKt]))/(1.00000+CONSTANTS[KmCaM]/STATES[cass]);
ALGEBRAIC[CaMKa] = ALGEBRAIC[CaMKb]+STATES[CaMKt];
ALGEBRAIC[ENa] =  (( CONSTANTS[R]*CONSTANTS[T])/CONSTANTS[F])*log(CONSTANTS[nao]/STATES[nai]);
ALGEBRAIC[INaL] =  CONSTANTS[GNaL]*(STATES[V] - ALGEBRAIC[ENa])*STATES[mL]*( (1.00000 - ALGEBRAIC[fINaLp])*STATES[hL]+ ALGEBRAIC[fINaLp]*STATES[hLp]);
ALGEBRAIC[mLss] = 1.00000/(1.00000+exp(- (STATES[V]+42.8500)/5.26400));
ALGEBRAIC[tm] = 1.00000/( CONSTANTS[mtD1]*exp((STATES[V]+CONSTANTS[mtV1])/CONSTANTS[mtV2])+ CONSTANTS[mtD2]*exp(- (STATES[V]+CONSTANTS[mtV3])/CONSTANTS[mtV4]));
ALGEBRAIC[tmL] = ALGEBRAIC[tm];
ALGEBRAIC[hLss] = 1.00000/(1.00000+exp((STATES[V]+87.6100)/7.48800));
ALGEBRAIC[hLssp] = 1.00000/(1.00000+exp((STATES[V]+93.8100)/7.48800));
ALGEBRAIC[JdiffNa] = (STATES[nass] - STATES[nai])/2.00000;
#endif
#if defined PCA_ORUDY2011 
ALGEBRAIC[IpCa] = ( CONSTANTS[GpCa]*STATES[cai])/(CONSTANTS[KmCap]+STATES[cai]);
ALGEBRAIC[Bcai] = 1.00000/(1.00000+( CONSTANTS[cmdnmax]*CONSTANTS[kmcmdn])/pow(CONSTANTS[kmcmdn]+STATES[cai], 2.00000)+( CONSTANTS[trpnmax]*CONSTANTS[kmtrpn])/pow(CONSTANTS[kmtrpn]+STATES[cai], 2.00000));
ALGEBRAIC[CaMKb] = ( CONSTANTS[CaMKo]*(1.00000 - STATES[CaMKt]))/(1.00000+CONSTANTS[KmCaM]/STATES[cass]);
ALGEBRAIC[CaMKa] = ALGEBRAIC[CaMKb]+STATES[CaMKt];
ALGEBRAIC[fJupp] = 1.00000/(1.00000+CONSTANTS[KmCaMK]/ALGEBRAIC[CaMKa]);
ALGEBRAIC[Jupnp] = ( CONSTANTS[upScale]*0.00437500*STATES[cai])/(STATES[cai]+0.000920000);
ALGEBRAIC[Jupp] = ( CONSTANTS[upScale]*2.75000*0.00437500*STATES[cai])/((STATES[cai]+0.000920000) - 0.000170000);
ALGEBRAIC[Jleak] = ( 0.00393750*STATES[cansr])/15.0000;
ALGEBRAIC[Jup] =  CONSTANTS[Jup_b]*(( (1.00000 - ALGEBRAIC[fJupp])*ALGEBRAIC[Jupnp]+ ALGEBRAIC[fJupp]*ALGEBRAIC[Jupp]) - ALGEBRAIC[Jleak]);
ALGEBRAIC[Bcajsr] = 1.00000/(1.00000+( CONSTANTS[csqnmax]*CONSTANTS[kmcsqn])/pow(CONSTANTS[kmcsqn]+STATES[cajsr], 2.00000));
ALGEBRAIC[fJrelp] = 1.00000/(1.00000+CONSTANTS[KmCaMK]/ALGEBRAIC[CaMKa]);
ALGEBRAIC[Jrel_inf_temp] = ( CONSTANTS[a_rel]*- ALGEBRAIC[ICaL])/(1.00000+ 1.00000*pow(1.50000/STATES[cajsr], 8.00000));
ALGEBRAIC[tau_rel_temp] = CONSTANTS[bt]/(1.00000+0.0123000/STATES[cajsr]);
ALGEBRAIC[Jrel_temp] = ( CONSTANTS[a_relp]*- ALGEBRAIC[ICaL])/(1.00000+pow(1.50000/STATES[cajsr], 8.00000));
ALGEBRAIC[tau_rel] = (ALGEBRAIC[tau_rel_temp]<0.00100000 ? 0.00100000 : ALGEBRAIC[tau_rel_temp]);
ALGEBRAIC[Jrel_inf] = (CONSTANTS[celltype]==2.00000 ?  ALGEBRAIC[Jrel_inf_temp]*1.70000 : ALGEBRAIC[Jrel_inf_temp]);
ALGEBRAIC[Jrel_infp] = (CONSTANTS[celltype]==2.00000 ?  ALGEBRAIC[Jrel_temp]*1.70000 : ALGEBRAIC[Jrel_temp]);
ALGEBRAIC[tau_relp] = (ALGEBRAIC[tau_relp_temp]<0.00100000 ? 0.00100000 : ALGEBRAIC[tau_relp_temp]);
ALGEBRAIC[Jrel] =  CONSTANTS[Jrel_scaling_factor]*( (1.00000 - ALGEBRAIC[fJrelp])*STATES[Jrelnp]+ ALGEBRAIC[fJrelp]*STATES[Jrelp]);
ALGEBRAIC[Jtr] = (STATES[cansr] - STATES[cajsr])/100.000;
#endif
#if defined TO_ORUDY2011 
ALGEBRAIC[EK] =  (( CONSTANTS[R]*CONSTANTS[T])/CONSTANTS[F])*log(CONSTANTS[ko]/STATES[ki]);
ALGEBRAIC[CaMKb] = ( CONSTANTS[CaMKo]*(1.00000 - STATES[CaMKt]))/(1.00000+CONSTANTS[KmCaM]/STATES[cass]);
ALGEBRAIC[CaMKa] = ALGEBRAIC[CaMKb]+STATES[CaMKt];
ALGEBRAIC[fItop] = 1.00000/(1.00000+CONSTANTS[KmCaMK]/ALGEBRAIC[CaMKa]);
ALGEBRAIC[AiF] = 1.00000/(1.00000+exp((STATES[V] - 213.600)/151.200));
ALGEBRAIC[AiS] = 1.00000 - ALGEBRAIC[AiF];
ALGEBRAIC[i] =  ALGEBRAIC[AiF]*STATES[iF]+ ALGEBRAIC[AiS]*STATES[iS];
ALGEBRAIC[ip] =  ALGEBRAIC[AiF]*STATES[iFp]+ ALGEBRAIC[AiS]*STATES[iSp];
ALGEBRAIC[Ito] =  CONSTANTS[Gto]*(STATES[V] - ALGEBRAIC[EK])*( (1.00000 - ALGEBRAIC[fItop])*STATES[a]*ALGEBRAIC[i]+ ALGEBRAIC[fItop]*STATES[ap]*ALGEBRAIC[ip]);
ALGEBRAIC[assp] = 1.00000/(1.00000+exp(- (STATES[V] - 24.3400)/14.8200));
ALGEBRAIC[ass] = 1.00000/(1.00000+exp(- (STATES[V] - 14.3400)/14.8200));
ALGEBRAIC[ta] = 1.05150/(1.00000/( 1.20890*(1.00000+exp(- (STATES[V] - 18.4099)/29.3814)))+3.50000/(1.00000+exp((STATES[V]+100.000)/29.3814)));
ALGEBRAIC[JdiffK] = (STATES[kss] - STATES[ki])/2.00000;
ALGEBRAIC[iss] = 1.00000/(1.00000+exp((STATES[V]+43.9400)/5.71100));
ALGEBRAIC[dti_recover] = 1.00000 - 0.500000/(1.00000+exp((STATES[V]+70.0000)/20.0000));
ALGEBRAIC[dti_develop] = 1.35400+0.000100000/(exp((STATES[V] - 167.400)/15.8900)+exp(- (STATES[V] - 12.2300)/0.215400));
ALGEBRAIC[tiF_b] = 4.56200+1.00000/( 0.393300*exp(- (STATES[V]+100.000)/100.000)+ 0.0800400*exp((STATES[V]+50.0000)/16.5900));
ALGEBRAIC[delta_epi] = (CONSTANTS[celltype]==1.00000 ? 1.00000 - 0.950000/(1.00000+exp((STATES[V]+70.0000)/5.00000)) : 1.00000);
ALGEBRAIC[tiS_b] = 23.6200+1.00000/( 0.00141600*exp(- (STATES[V]+96.5200)/59.0500)+ 1.78000e-08*exp((STATES[V]+114.100)/8.07900));
ALGEBRAIC[tiF] =  ALGEBRAIC[tiF_b]*ALGEBRAIC[delta_epi];
ALGEBRAIC[tiS] =  ALGEBRAIC[tiS_b]*ALGEBRAIC[delta_epi];
ALGEBRAIC[tiFp] =  ALGEBRAIC[dti_develop]*ALGEBRAIC[dti_recover]*ALGEBRAIC[tiF];
ALGEBRAIC[tiSp] =  ALGEBRAIC[dti_develop]*ALGEBRAIC[dti_recover]*ALGEBRAIC[tiS];
#endif


#if defined SINGLE_CELL
RATES[V] = - (ALGEBRAIC[INa]+ALGEBRAIC[INaL]+ALGEBRAIC[Ito]+ALGEBRAIC[ICaL]+ALGEBRAIC[ICaNa]+ALGEBRAIC[ICaK]+ALGEBRAIC[IKr]+ALGEBRAIC[IKs]+ALGEBRAIC[IK1]+ALGEBRAIC[INaCa_i]+ALGEBRAIC[INaCa_ss]+ALGEBRAIC[INaK]+ALGEBRAIC[INab]+ALGEBRAIC[IKb]+ALGEBRAIC[IpCa]+ALGEBRAIC[ICab]+ALGEBRAIC[Istim]);
#else
RATES[V] = 0.;
#endif
#if defined CAB_ORUDY2011 
RATES[cai] =  ALGEBRAIC[Bcai]*((( - ((ALGEBRAIC[IpCa]+ALGEBRAIC[ICab]) -  2.00000*ALGEBRAIC[INaCa_i])*CONSTANTS[cm]*CONSTANTS[Acap])/( 2.00000*CONSTANTS[F]*CONSTANTS[vmyo]) - ( ALGEBRAIC[Jup]*CONSTANTS[vnsr])/CONSTANTS[vmyo])+( ALGEBRAIC[Jdiff]*CONSTANTS[vss])/CONSTANTS[vmyo]);
RATES[cansr] = ALGEBRAIC[Jup] - ( ALGEBRAIC[Jtr]*CONSTANTS[vjsr])/CONSTANTS[vnsr];
RATES[CaMKt] =  CONSTANTS[aCaMK]*ALGEBRAIC[CaMKb]*(ALGEBRAIC[CaMKb]+STATES[CaMKt]) -  CONSTANTS[bCaMK]*STATES[CaMKt];
RATES[cajsr] =  ALGEBRAIC[Bcajsr]*(ALGEBRAIC[Jtr] - ALGEBRAIC[Jrel]);
RATES[Jrelnp] = (ALGEBRAIC[Jrel_inf] - STATES[Jrelnp])/ALGEBRAIC[tau_rel];
RATES[Jrelp] = (ALGEBRAIC[Jrel_infp] - STATES[Jrelp])/ALGEBRAIC[tau_relp];
RATES[cass] =  ALGEBRAIC[Bcass]*((( - (ALGEBRAIC[ICaL] -  2.00000*ALGEBRAIC[INaCa_ss])*CONSTANTS[cm]*CONSTANTS[Acap])/( 2.00000*CONSTANTS[F]*CONSTANTS[vss])+( ALGEBRAIC[Jrel]*CONSTANTS[vjsr])/CONSTANTS[vss]) - ALGEBRAIC[Jdiff]);
#endif
#if defined CAK_ORUDY2011 
RATES[cass] =  ALGEBRAIC[Bcass]*((( - (ALGEBRAIC[ICaL] -  2.00000*ALGEBRAIC[INaCa_ss])*CONSTANTS[cm]*CONSTANTS[Acap])/( 2.00000*CONSTANTS[F]*CONSTANTS[vss])+( ALGEBRAIC[Jrel]*CONSTANTS[vjsr])/CONSTANTS[vss]) - ALGEBRAIC[Jdiff]);
RATES[d] = (ALGEBRAIC[dss] - STATES[d])/ALGEBRAIC[td];
RATES[CaMKt] =  CONSTANTS[aCaMK]*ALGEBRAIC[CaMKb]*(ALGEBRAIC[CaMKb]+STATES[CaMKt]) -  CONSTANTS[bCaMK]*STATES[CaMKt];
RATES[ff] = (ALGEBRAIC[fss] - STATES[ff])/ALGEBRAIC[tff];
RATES[fs] = (ALGEBRAIC[fss] - STATES[fs])/ALGEBRAIC[tfs];
RATES[fcaf] = (ALGEBRAIC[fcass] - STATES[fcaf])/ALGEBRAIC[tfcaf];
RATES[nca] =  ALGEBRAIC[anca]*CONSTANTS[k2n] -  STATES[nca]*ALGEBRAIC[km2n];
RATES[jca] = (ALGEBRAIC[fcass] - STATES[jca])/CONSTANTS[tjca];
RATES[fcas] = (ALGEBRAIC[fcass] - STATES[fcas])/ALGEBRAIC[tfcas];
RATES[ffp] = (ALGEBRAIC[fss] - STATES[ffp])/ALGEBRAIC[tffp];
RATES[fcafp] = (ALGEBRAIC[fcass] - STATES[fcafp])/ALGEBRAIC[tfcafp];
RATES[Jrelnp] = (ALGEBRAIC[Jrel_inf] - STATES[Jrelnp])/ALGEBRAIC[tau_rel];
RATES[Jrelp] = (ALGEBRAIC[Jrel_infp] - STATES[Jrelp])/ALGEBRAIC[tau_relp];
RATES[cai] =  ALGEBRAIC[Bcai]*((( - ((ALGEBRAIC[IpCa]+ALGEBRAIC[ICab]) -  2.00000*ALGEBRAIC[INaCa_i])*CONSTANTS[cm]*CONSTANTS[Acap])/( 2.00000*CONSTANTS[F]*CONSTANTS[vmyo]) - ( ALGEBRAIC[Jup]*CONSTANTS[vnsr])/CONSTANTS[vmyo])+( ALGEBRAIC[Jdiff]*CONSTANTS[vss])/CONSTANTS[vmyo]);
RATES[cansr] = ALGEBRAIC[Jup] - ( ALGEBRAIC[Jtr]*CONSTANTS[vjsr])/CONSTANTS[vnsr];
RATES[cajsr] =  ALGEBRAIC[Bcajsr]*(ALGEBRAIC[Jtr] - ALGEBRAIC[Jrel]);
RATES[ki] = ( - ((ALGEBRAIC[Ito]+ALGEBRAIC[IKr]+ALGEBRAIC[IKs]+ALGEBRAIC[IK1]+ALGEBRAIC[IKb]+ALGEBRAIC[Istim]) -  2.00000*ALGEBRAIC[INaK])*CONSTANTS[cm]*CONSTANTS[Acap])/( CONSTANTS[F]*CONSTANTS[vmyo])+( ALGEBRAIC[JdiffK]*CONSTANTS[vss])/CONSTANTS[vmyo];
RATES[kss] = ( - ALGEBRAIC[ICaK]*CONSTANTS[cm]*CONSTANTS[Acap])/( CONSTANTS[F]*CONSTANTS[vss]) - ALGEBRAIC[JdiffK];
#endif
#if defined CAL_ORUDY2011 
RATES[cass] =  ALGEBRAIC[Bcass]*((( - (ALGEBRAIC[ICaL] -  2.00000*ALGEBRAIC[INaCa_ss])*CONSTANTS[cm]*CONSTANTS[Acap])/( 2.00000*CONSTANTS[F]*CONSTANTS[vss])+( ALGEBRAIC[Jrel]*CONSTANTS[vjsr])/CONSTANTS[vss]) - ALGEBRAIC[Jdiff]);
RATES[d] = (ALGEBRAIC[dss] - STATES[d])/ALGEBRAIC[td];
RATES[CaMKt] =  CONSTANTS[aCaMK]*ALGEBRAIC[CaMKb]*(ALGEBRAIC[CaMKb]+STATES[CaMKt]) -  CONSTANTS[bCaMK]*STATES[CaMKt];
RATES[ff] = (ALGEBRAIC[fss] - STATES[ff])/ALGEBRAIC[tff];
RATES[fs] = (ALGEBRAIC[fss] - STATES[fs])/ALGEBRAIC[tfs];
RATES[fcaf] = (ALGEBRAIC[fcass] - STATES[fcaf])/ALGEBRAIC[tfcaf];
RATES[jca] = (ALGEBRAIC[fcass] - STATES[jca])/CONSTANTS[tjca];
RATES[nca] =  ALGEBRAIC[anca]*CONSTANTS[k2n] -  STATES[nca]*ALGEBRAIC[km2n];
RATES[fcas] = (ALGEBRAIC[fcass] - STATES[fcas])/ALGEBRAIC[tfcas];
RATES[ffp] = (ALGEBRAIC[fss] - STATES[ffp])/ALGEBRAIC[tffp];
RATES[fcafp] = (ALGEBRAIC[fcass] - STATES[fcafp])/ALGEBRAIC[tfcafp];
RATES[Jrelnp] = (ALGEBRAIC[Jrel_inf] - STATES[Jrelnp])/ALGEBRAIC[tau_rel];
RATES[Jrelp] = (ALGEBRAIC[Jrel_infp] - STATES[Jrelp])/ALGEBRAIC[tau_relp];
RATES[cai] =  ALGEBRAIC[Bcai]*((( - ((ALGEBRAIC[IpCa]+ALGEBRAIC[ICab]) -  2.00000*ALGEBRAIC[INaCa_i])*CONSTANTS[cm]*CONSTANTS[Acap])/( 2.00000*CONSTANTS[F]*CONSTANTS[vmyo]) - ( ALGEBRAIC[Jup]*CONSTANTS[vnsr])/CONSTANTS[vmyo])+( ALGEBRAIC[Jdiff]*CONSTANTS[vss])/CONSTANTS[vmyo]);
RATES[cansr] = ALGEBRAIC[Jup] - ( ALGEBRAIC[Jtr]*CONSTANTS[vjsr])/CONSTANTS[vnsr];
RATES[cajsr] =  ALGEBRAIC[Bcajsr]*(ALGEBRAIC[Jtr] - ALGEBRAIC[Jrel]);
#endif
#if defined CANA_ORUDY2011 
RATES[d] = (ALGEBRAIC[dss] - STATES[d])/ALGEBRAIC[td];
RATES[cass] =  ALGEBRAIC[Bcass]*((( - (ALGEBRAIC[ICaL] -  2.00000*ALGEBRAIC[INaCa_ss])*CONSTANTS[cm]*CONSTANTS[Acap])/( 2.00000*CONSTANTS[F]*CONSTANTS[vss])+( ALGEBRAIC[Jrel]*CONSTANTS[vjsr])/CONSTANTS[vss]) - ALGEBRAIC[Jdiff]);
RATES[ff] = (ALGEBRAIC[fss] - STATES[ff])/ALGEBRAIC[tff];
RATES[fs] = (ALGEBRAIC[fss] - STATES[fs])/ALGEBRAIC[tfs];
RATES[nca] =  ALGEBRAIC[anca]*CONSTANTS[k2n] -  STATES[nca]*ALGEBRAIC[km2n];
RATES[jca] = (ALGEBRAIC[fcass] - STATES[jca])/CONSTANTS[tjca];
RATES[fcaf] = (ALGEBRAIC[fcass] - STATES[fcaf])/ALGEBRAIC[tfcaf];
RATES[ffp] = (ALGEBRAIC[fss] - STATES[ffp])/ALGEBRAIC[tffp];
RATES[fcafp] = (ALGEBRAIC[fcass] - STATES[fcafp])/ALGEBRAIC[tfcafp];
RATES[fcas] = (ALGEBRAIC[fcass] - STATES[fcas])/ALGEBRAIC[tfcas];
RATES[CaMKt] =  CONSTANTS[aCaMK]*ALGEBRAIC[CaMKb]*(ALGEBRAIC[CaMKb]+STATES[CaMKt]) -  CONSTANTS[bCaMK]*STATES[CaMKt];
RATES[nass] = ( - (ALGEBRAIC[ICaNa]+ 3.00000*ALGEBRAIC[INaCa_ss])*CONSTANTS[cm]*CONSTANTS[Acap])/( CONSTANTS[F]*CONSTANTS[vss]) - ALGEBRAIC[JdiffNa];
RATES[nai] = ( - (ALGEBRAIC[INa]+ALGEBRAIC[INaL]+ 3.00000*ALGEBRAIC[INaCa_i]+ 3.00000*ALGEBRAIC[INaK]+ALGEBRAIC[INab])*CONSTANTS[Acap]*CONSTANTS[cm])/( CONSTANTS[F]*CONSTANTS[vmyo])+( ALGEBRAIC[JdiffNa]*CONSTANTS[vss])/CONSTANTS[vmyo];
RATES[cansr] = ALGEBRAIC[Jup] - ( ALGEBRAIC[Jtr]*CONSTANTS[vjsr])/CONSTANTS[vnsr];
RATES[cajsr] =  ALGEBRAIC[Bcajsr]*(ALGEBRAIC[Jtr] - ALGEBRAIC[Jrel]);
RATES[cai] =  ALGEBRAIC[Bcai]*((( - ((ALGEBRAIC[IpCa]+ALGEBRAIC[ICab]) -  2.00000*ALGEBRAIC[INaCa_i])*CONSTANTS[cm]*CONSTANTS[Acap])/( 2.00000*CONSTANTS[F]*CONSTANTS[vmyo]) - ( ALGEBRAIC[Jup]*CONSTANTS[vnsr])/CONSTANTS[vmyo])+( ALGEBRAIC[Jdiff]*CONSTANTS[vss])/CONSTANTS[vmyo]);
RATES[Jrelnp] = (ALGEBRAIC[Jrel_inf] - STATES[Jrelnp])/ALGEBRAIC[tau_rel];
RATES[Jrelp] = (ALGEBRAIC[Jrel_infp] - STATES[Jrelp])/ALGEBRAIC[tau_relp];
#endif
#if defined K1_ORUDY2011 
RATES[xk1] = (ALGEBRAIC[xk1ss] - STATES[xk1])/ALGEBRAIC[txk1];
RATES[ki] = ( - ((ALGEBRAIC[Ito]+ALGEBRAIC[IKr]+ALGEBRAIC[IKs]+ALGEBRAIC[IK1]+ALGEBRAIC[IKb]+ALGEBRAIC[Istim]) -  2.00000*ALGEBRAIC[INaK])*CONSTANTS[cm]*CONSTANTS[Acap])/( CONSTANTS[F]*CONSTANTS[vmyo])+( ALGEBRAIC[JdiffK]*CONSTANTS[vss])/CONSTANTS[vmyo];
RATES[kss] = ( - ALGEBRAIC[ICaK]*CONSTANTS[cm]*CONSTANTS[Acap])/( CONSTANTS[F]*CONSTANTS[vss]) - ALGEBRAIC[JdiffK];
#endif
#if defined KB_ORUDY2011 
RATES[ki] = ( - ((ALGEBRAIC[Ito]+ALGEBRAIC[IKr]+ALGEBRAIC[IKs]+ALGEBRAIC[IK1]+ALGEBRAIC[IKb]+ALGEBRAIC[Istim]) -  2.00000*ALGEBRAIC[INaK])*CONSTANTS[cm]*CONSTANTS[Acap])/( CONSTANTS[F]*CONSTANTS[vmyo])+( ALGEBRAIC[JdiffK]*CONSTANTS[vss])/CONSTANTS[vmyo];
RATES[kss] = ( - ALGEBRAIC[ICaK]*CONSTANTS[cm]*CONSTANTS[Acap])/( CONSTANTS[F]*CONSTANTS[vss]) - ALGEBRAIC[JdiffK];
#endif
#if defined KR_ORUDY2011 
RATES[ki] = ( - ((ALGEBRAIC[Ito]+ALGEBRAIC[IKr]+ALGEBRAIC[IKs]+ALGEBRAIC[IK1]+ALGEBRAIC[IKb]+ALGEBRAIC[Istim]) -  2.00000*ALGEBRAIC[INaK])*CONSTANTS[cm]*CONSTANTS[Acap])/( CONSTANTS[F]*CONSTANTS[vmyo])+( ALGEBRAIC[JdiffK]*CONSTANTS[vss])/CONSTANTS[vmyo];
RATES[kss] = ( - ALGEBRAIC[ICaK]*CONSTANTS[cm]*CONSTANTS[Acap])/( CONSTANTS[F]*CONSTANTS[vss]) - ALGEBRAIC[JdiffK];
RATES[IC1] = (- ( CONSTANTS[A11]*exp( CONSTANTS[B11]*STATES[V])*STATES[IC1]*exp(( (CONSTANTS[Temp] - 20.0000)*log(CONSTANTS[q11]))/10.0000) -  CONSTANTS[A21]*exp( CONSTANTS[B21]*STATES[V])*STATES[IC2]*exp(( (CONSTANTS[Temp] - 20.0000)*log(CONSTANTS[q21]))/10.0000))+ CONSTANTS[A51]*exp( CONSTANTS[B51]*STATES[V])*STATES[C1]*exp(( (CONSTANTS[Temp] - 20.0000)*log(CONSTANTS[q51]))/10.0000)) -  CONSTANTS[A61]*exp( CONSTANTS[B61]*STATES[V])*STATES[IC1]*exp(( (CONSTANTS[Temp] - 20.0000)*log(CONSTANTS[q61]))/10.0000);
RATES[IC2] = ((( CONSTANTS[A11]*exp( CONSTANTS[B11]*STATES[V])*STATES[IC1]*exp(( (CONSTANTS[Temp] - 20.0000)*log(CONSTANTS[q11]))/10.0000) -  CONSTANTS[A21]*exp( CONSTANTS[B21]*STATES[V])*STATES[IC2]*exp(( (CONSTANTS[Temp] - 20.0000)*log(CONSTANTS[q21]))/10.0000)) - ( CONSTANTS[A3]*exp( CONSTANTS[B3]*STATES[V])*STATES[IC2]*exp(( (CONSTANTS[Temp] - 20.0000)*log(CONSTANTS[q3]))/10.0000) -  CONSTANTS[A4]*exp( CONSTANTS[B4]*STATES[V])*STATES[IO]*exp(( (CONSTANTS[Temp] - 20.0000)*log(CONSTANTS[q4]))/10.0000)))+ CONSTANTS[A52]*exp( CONSTANTS[B52]*STATES[V])*STATES[C2]*exp(( (CONSTANTS[Temp] - 20.0000)*log(CONSTANTS[q52]))/10.0000)) -  CONSTANTS[A62]*exp( CONSTANTS[B62]*STATES[V])*STATES[IC2]*exp(( (CONSTANTS[Temp] - 20.0000)*log(CONSTANTS[q62]))/10.0000);
RATES[C1] = - ( CONSTANTS[A1]*exp( CONSTANTS[B1]*STATES[V])*STATES[C1]*exp(( (CONSTANTS[Temp] - 20.0000)*log(CONSTANTS[q1]))/10.0000) -  CONSTANTS[A2]*exp( CONSTANTS[B2]*STATES[V])*STATES[C2]*exp(( (CONSTANTS[Temp] - 20.0000)*log(CONSTANTS[q2]))/10.0000)) - ( CONSTANTS[A51]*exp( CONSTANTS[B51]*STATES[V])*STATES[C1]*exp(( (CONSTANTS[Temp] - 20.0000)*log(CONSTANTS[q51]))/10.0000) -  CONSTANTS[A61]*exp( CONSTANTS[B61]*STATES[V])*STATES[IC1]*exp(( (CONSTANTS[Temp] - 20.0000)*log(CONSTANTS[q61]))/10.0000));
RATES[C2] = (( CONSTANTS[A1]*exp( CONSTANTS[B1]*STATES[V])*STATES[C1]*exp(( (CONSTANTS[Temp] - 20.0000)*log(CONSTANTS[q1]))/10.0000) -  CONSTANTS[A2]*exp( CONSTANTS[B2]*STATES[V])*STATES[C2]*exp(( (CONSTANTS[Temp] - 20.0000)*log(CONSTANTS[q2]))/10.0000)) - ( CONSTANTS[A31]*exp( CONSTANTS[B31]*STATES[V])*STATES[C2]*exp(( (CONSTANTS[Temp] - 20.0000)*log(CONSTANTS[q31]))/10.0000) -  CONSTANTS[A41]*exp( CONSTANTS[B41]*STATES[V])*STATES[O]*exp(( (CONSTANTS[Temp] - 20.0000)*log(CONSTANTS[q41]))/10.0000))) - ( CONSTANTS[A52]*exp( CONSTANTS[B52]*STATES[V])*STATES[C2]*exp(( (CONSTANTS[Temp] - 20.0000)*log(CONSTANTS[q52]))/10.0000) -  CONSTANTS[A62]*exp( CONSTANTS[B62]*STATES[V])*STATES[IC2]*exp(( (CONSTANTS[Temp] - 20.0000)*log(CONSTANTS[q62]))/10.0000));
RATES[O] = (( CONSTANTS[A31]*exp( CONSTANTS[B31]*STATES[V])*STATES[C2]*exp(( (CONSTANTS[Temp] - 20.0000)*log(CONSTANTS[q31]))/10.0000) -  CONSTANTS[A41]*exp( CONSTANTS[B41]*STATES[V])*STATES[O]*exp(( (CONSTANTS[Temp] - 20.0000)*log(CONSTANTS[q41]))/10.0000)) - ( CONSTANTS[A53]*exp( CONSTANTS[B53]*STATES[V])*STATES[O]*exp(( (CONSTANTS[Temp] - 20.0000)*log(CONSTANTS[q53]))/10.0000) -  CONSTANTS[A63]*exp( CONSTANTS[B63]*STATES[V])*STATES[IO]*exp(( (CONSTANTS[Temp] - 20.0000)*log(CONSTANTS[q63]))/10.0000))) - ( (( CONSTANTS[Kmax]*CONSTANTS[Ku]*exp( CONSTANTS[n]*log(STATES[D])))/(exp( CONSTANTS[n]*log(STATES[D]))+CONSTANTS[halfmax]))*STATES[O] -  CONSTANTS[Ku]*STATES[Obound]);
RATES[IO] = ((( CONSTANTS[A3]*exp( CONSTANTS[B3]*STATES[V])*STATES[IC2]*exp(( (CONSTANTS[Temp] - 20.0000)*log(CONSTANTS[q3]))/10.0000) -  CONSTANTS[A4]*exp( CONSTANTS[B4]*STATES[V])*STATES[IO]*exp(( (CONSTANTS[Temp] - 20.0000)*log(CONSTANTS[q4]))/10.0000))+ CONSTANTS[A53]*exp( CONSTANTS[B53]*STATES[V])*STATES[O]*exp(( (CONSTANTS[Temp] - 20.0000)*log(CONSTANTS[q53]))/10.0000)) -  CONSTANTS[A63]*exp( CONSTANTS[B63]*STATES[V])*STATES[IO]*exp(( (CONSTANTS[Temp] - 20.0000)*log(CONSTANTS[q63]))/10.0000)) - ( (( CONSTANTS[Kmax]*CONSTANTS[Ku]*exp( CONSTANTS[n]*log(STATES[D])))/(exp( CONSTANTS[n]*log(STATES[D]))+CONSTANTS[halfmax]))*STATES[IO] -  (( CONSTANTS[Ku]*CONSTANTS[A53]*exp( CONSTANTS[B53]*STATES[V])*exp(( (CONSTANTS[Temp] - 20.0000)*log(CONSTANTS[q53]))/10.0000))/( CONSTANTS[A63]*exp( CONSTANTS[B63]*STATES[V])*exp(( (CONSTANTS[Temp] - 20.0000)*log(CONSTANTS[q63]))/10.0000)))*STATES[IObound]);
RATES[D] = 0.;
RATES[IObound] = (( (( CONSTANTS[Kmax]*CONSTANTS[Ku]*exp( CONSTANTS[n]*log(STATES[D])))/(exp( CONSTANTS[n]*log(STATES[D]))+CONSTANTS[halfmax]))*STATES[IO] -  (( CONSTANTS[Ku]*CONSTANTS[A53]*exp( CONSTANTS[B53]*STATES[V])*exp(( (CONSTANTS[Temp] - 20.0000)*log(CONSTANTS[q53]))/10.0000))/( CONSTANTS[A63]*exp( CONSTANTS[B63]*STATES[V])*exp(( (CONSTANTS[Temp] - 20.0000)*log(CONSTANTS[q63]))/10.0000)))*STATES[IObound])+ (CONSTANTS[Kt]/(1.00000+exp(- (STATES[V] - CONSTANTS[Vhalf])/6.78900)))*STATES[Cbound]) -  CONSTANTS[Kt]*STATES[IObound];
RATES[Obound] = (( (( CONSTANTS[Kmax]*CONSTANTS[Ku]*exp( CONSTANTS[n]*log(STATES[D])))/(exp( CONSTANTS[n]*log(STATES[D]))+CONSTANTS[halfmax]))*STATES[O] -  CONSTANTS[Ku]*STATES[Obound])+ (CONSTANTS[Kt]/(1.00000+exp(- (STATES[V] - CONSTANTS[Vhalf])/6.78900)))*STATES[Cbound]) -  CONSTANTS[Kt]*STATES[Obound];
RATES[Cbound] = - ( (CONSTANTS[Kt]/(1.00000+exp(- (STATES[V] - CONSTANTS[Vhalf])/6.78900)))*STATES[Cbound] -  CONSTANTS[Kt]*STATES[Obound]) - ( (CONSTANTS[Kt]/(1.00000+exp(- (STATES[V] - CONSTANTS[Vhalf])/6.78900)))*STATES[Cbound] -  CONSTANTS[Kt]*STATES[IObound]);
#endif
#if defined KS_ORUDY2011 
RATES[xs1] = (ALGEBRAIC[xs1ss] - STATES[xs1])/ALGEBRAIC[txs1];
RATES[xs2] = (ALGEBRAIC[xs2ss] - STATES[xs2])/ALGEBRAIC[txs2];
RATES[ki] = ( - ((ALGEBRAIC[Ito]+ALGEBRAIC[IKr]+ALGEBRAIC[IKs]+ALGEBRAIC[IK1]+ALGEBRAIC[IKb]+ALGEBRAIC[Istim]) -  2.00000*ALGEBRAIC[INaK])*CONSTANTS[cm]*CONSTANTS[Acap])/( CONSTANTS[F]*CONSTANTS[vmyo])+( ALGEBRAIC[JdiffK]*CONSTANTS[vss])/CONSTANTS[vmyo];
RATES[kss] = ( - ALGEBRAIC[ICaK]*CONSTANTS[cm]*CONSTANTS[Acap])/( CONSTANTS[F]*CONSTANTS[vss]) - ALGEBRAIC[JdiffK];
RATES[CaMKt] =  CONSTANTS[aCaMK]*ALGEBRAIC[CaMKb]*(ALGEBRAIC[CaMKb]+STATES[CaMKt]) -  CONSTANTS[bCaMK]*STATES[CaMKt];
#endif
#if defined NA_ORUDY2011 
RATES[m] = (ALGEBRAIC[mss] - STATES[m])/ALGEBRAIC[tm];
RATES[j] = (ALGEBRAIC[jss] - STATES[j])/ALGEBRAIC[tj];
RATES[jp] = (ALGEBRAIC[jss] - STATES[jp])/ALGEBRAIC[tjp];
RATES[hf] = (ALGEBRAIC[hss] - STATES[hf])/ALGEBRAIC[thf];
RATES[hs] = (ALGEBRAIC[hss] - STATES[hs])/ALGEBRAIC[ths];
RATES[hsp] = (ALGEBRAIC[hssp] - STATES[hsp])/ALGEBRAIC[thsp];
RATES[nai] = ( - (ALGEBRAIC[INa]+ALGEBRAIC[INaL]+ 3.00000*ALGEBRAIC[INaCa_i]+ 3.00000*ALGEBRAIC[INaK]+ALGEBRAIC[INab])*CONSTANTS[Acap]*CONSTANTS[cm])/( CONSTANTS[F]*CONSTANTS[vmyo])+( ALGEBRAIC[JdiffNa]*CONSTANTS[vss])/CONSTANTS[vmyo];
RATES[nass] = ( - (ALGEBRAIC[ICaNa]+ 3.00000*ALGEBRAIC[INaCa_ss])*CONSTANTS[cm]*CONSTANTS[Acap])/( CONSTANTS[F]*CONSTANTS[vss]) - ALGEBRAIC[JdiffNa];
RATES[CaMKt] =  CONSTANTS[aCaMK]*ALGEBRAIC[CaMKb]*(ALGEBRAIC[CaMKb]+STATES[CaMKt]) -  CONSTANTS[bCaMK]*STATES[CaMKt];
#endif
#if defined NAB_ORUDY2011 
RATES[nass] = ( - (ALGEBRAIC[ICaNa]+ 3.00000*ALGEBRAIC[INaCa_ss])*CONSTANTS[cm]*CONSTANTS[Acap])/( CONSTANTS[F]*CONSTANTS[vss]) - ALGEBRAIC[JdiffNa];
RATES[nai] = ( - (ALGEBRAIC[INa]+ALGEBRAIC[INaL]+ 3.00000*ALGEBRAIC[INaCa_i]+ 3.00000*ALGEBRAIC[INaK]+ALGEBRAIC[INab])*CONSTANTS[Acap]*CONSTANTS[cm])/( CONSTANTS[F]*CONSTANTS[vmyo])+( ALGEBRAIC[JdiffNa]*CONSTANTS[vss])/CONSTANTS[vmyo];
#endif
#if defined NACA_I_ORUDY2011 
RATES[cai] =  ALGEBRAIC[Bcai]*((( - ((ALGEBRAIC[IpCa]+ALGEBRAIC[ICab]) -  2.00000*ALGEBRAIC[INaCa_i])*CONSTANTS[cm]*CONSTANTS[Acap])/( 2.00000*CONSTANTS[F]*CONSTANTS[vmyo]) - ( ALGEBRAIC[Jup]*CONSTANTS[vnsr])/CONSTANTS[vmyo])+( ALGEBRAIC[Jdiff]*CONSTANTS[vss])/CONSTANTS[vmyo]);
RATES[nai] = ( - (ALGEBRAIC[INa]+ALGEBRAIC[INaL]+ 3.00000*ALGEBRAIC[INaCa_i]+ 3.00000*ALGEBRAIC[INaK]+ALGEBRAIC[INab])*CONSTANTS[Acap]*CONSTANTS[cm])/( CONSTANTS[F]*CONSTANTS[vmyo])+( ALGEBRAIC[JdiffNa]*CONSTANTS[vss])/CONSTANTS[vmyo];
RATES[cansr] = ALGEBRAIC[Jup] - ( ALGEBRAIC[Jtr]*CONSTANTS[vjsr])/CONSTANTS[vnsr];
RATES[cass] =  ALGEBRAIC[Bcass]*((( - (ALGEBRAIC[ICaL] -  2.00000*ALGEBRAIC[INaCa_ss])*CONSTANTS[cm]*CONSTANTS[Acap])/( 2.00000*CONSTANTS[F]*CONSTANTS[vss])+( ALGEBRAIC[Jrel]*CONSTANTS[vjsr])/CONSTANTS[vss]) - ALGEBRAIC[Jdiff]);
RATES[nass] = ( - (ALGEBRAIC[ICaNa]+ 3.00000*ALGEBRAIC[INaCa_ss])*CONSTANTS[cm]*CONSTANTS[Acap])/( CONSTANTS[F]*CONSTANTS[vss]) - ALGEBRAIC[JdiffNa];
RATES[CaMKt] =  CONSTANTS[aCaMK]*ALGEBRAIC[CaMKb]*(ALGEBRAIC[CaMKb]+STATES[CaMKt]) -  CONSTANTS[bCaMK]*STATES[CaMKt];
RATES[cajsr] =  ALGEBRAIC[Bcajsr]*(ALGEBRAIC[Jtr] - ALGEBRAIC[Jrel]);
RATES[Jrelnp] = (ALGEBRAIC[Jrel_inf] - STATES[Jrelnp])/ALGEBRAIC[tau_rel];
RATES[Jrelp] = (ALGEBRAIC[Jrel_infp] - STATES[Jrelp])/ALGEBRAIC[tau_relp];
#endif
#if defined NACA_SS_ORUDY2011 
RATES[cai] =  ALGEBRAIC[Bcai]*((( - ((ALGEBRAIC[IpCa]+ALGEBRAIC[ICab]) -  2.00000*ALGEBRAIC[INaCa_i])*CONSTANTS[cm]*CONSTANTS[Acap])/( 2.00000*CONSTANTS[F]*CONSTANTS[vmyo]) - ( ALGEBRAIC[Jup]*CONSTANTS[vnsr])/CONSTANTS[vmyo])+( ALGEBRAIC[Jdiff]*CONSTANTS[vss])/CONSTANTS[vmyo]);
RATES[nai] = ( - (ALGEBRAIC[INa]+ALGEBRAIC[INaL]+ 3.00000*ALGEBRAIC[INaCa_i]+ 3.00000*ALGEBRAIC[INaK]+ALGEBRAIC[INab])*CONSTANTS[Acap]*CONSTANTS[cm])/( CONSTANTS[F]*CONSTANTS[vmyo])+( ALGEBRAIC[JdiffNa]*CONSTANTS[vss])/CONSTANTS[vmyo];
RATES[cansr] = ALGEBRAIC[Jup] - ( ALGEBRAIC[Jtr]*CONSTANTS[vjsr])/CONSTANTS[vnsr];
RATES[cass] =  ALGEBRAIC[Bcass]*((( - (ALGEBRAIC[ICaL] -  2.00000*ALGEBRAIC[INaCa_ss])*CONSTANTS[cm]*CONSTANTS[Acap])/( 2.00000*CONSTANTS[F]*CONSTANTS[vss])+( ALGEBRAIC[Jrel]*CONSTANTS[vjsr])/CONSTANTS[vss]) - ALGEBRAIC[Jdiff]);
RATES[nass] = ( - (ALGEBRAIC[ICaNa]+ 3.00000*ALGEBRAIC[INaCa_ss])*CONSTANTS[cm]*CONSTANTS[Acap])/( CONSTANTS[F]*CONSTANTS[vss]) - ALGEBRAIC[JdiffNa];
RATES[CaMKt] =  CONSTANTS[aCaMK]*ALGEBRAIC[CaMKb]*(ALGEBRAIC[CaMKb]+STATES[CaMKt]) -  CONSTANTS[bCaMK]*STATES[CaMKt];
RATES[cajsr] =  ALGEBRAIC[Bcajsr]*(ALGEBRAIC[Jtr] - ALGEBRAIC[Jrel]);
RATES[Jrelnp] = (ALGEBRAIC[Jrel_inf] - STATES[Jrelnp])/ALGEBRAIC[tau_rel];
RATES[Jrelp] = (ALGEBRAIC[Jrel_infp] - STATES[Jrelp])/ALGEBRAIC[tau_relp];
#endif
#if defined NAK_ORUDY2011 
RATES[nai] = ( - (ALGEBRAIC[INa]+ALGEBRAIC[INaL]+ 3.00000*ALGEBRAIC[INaCa_i]+ 3.00000*ALGEBRAIC[INaK]+ALGEBRAIC[INab])*CONSTANTS[Acap]*CONSTANTS[cm])/( CONSTANTS[F]*CONSTANTS[vmyo])+( ALGEBRAIC[JdiffNa]*CONSTANTS[vss])/CONSTANTS[vmyo];
RATES[ki] = ( - ((ALGEBRAIC[Ito]+ALGEBRAIC[IKr]+ALGEBRAIC[IKs]+ALGEBRAIC[IK1]+ALGEBRAIC[IKb]+ALGEBRAIC[Istim]) -  2.00000*ALGEBRAIC[INaK])*CONSTANTS[cm]*CONSTANTS[Acap])/( CONSTANTS[F]*CONSTANTS[vmyo])+( ALGEBRAIC[JdiffK]*CONSTANTS[vss])/CONSTANTS[vmyo];
RATES[kss] = ( - ALGEBRAIC[ICaK]*CONSTANTS[cm]*CONSTANTS[Acap])/( CONSTANTS[F]*CONSTANTS[vss]) - ALGEBRAIC[JdiffK];
RATES[nass] = ( - (ALGEBRAIC[ICaNa]+ 3.00000*ALGEBRAIC[INaCa_ss])*CONSTANTS[cm]*CONSTANTS[Acap])/( CONSTANTS[F]*CONSTANTS[vss]) - ALGEBRAIC[JdiffNa];
#endif
#if defined NAL_ORUDY2011 
RATES[mL] = (ALGEBRAIC[mLss] - STATES[mL])/ALGEBRAIC[tmL];
RATES[hL] = (ALGEBRAIC[hLss] - STATES[hL])/CONSTANTS[thL];
RATES[hLp] = (ALGEBRAIC[hLssp] - STATES[hLp])/CONSTANTS[thLp];
RATES[CaMKt] =  CONSTANTS[aCaMK]*ALGEBRAIC[CaMKb]*(ALGEBRAIC[CaMKb]+STATES[CaMKt]) -  CONSTANTS[bCaMK]*STATES[CaMKt];
RATES[nai] = ( - (ALGEBRAIC[INa]+ALGEBRAIC[INaL]+ 3.00000*ALGEBRAIC[INaCa_i]+ 3.00000*ALGEBRAIC[INaK]+ALGEBRAIC[INab])*CONSTANTS[Acap]*CONSTANTS[cm])/( CONSTANTS[F]*CONSTANTS[vmyo])+( ALGEBRAIC[JdiffNa]*CONSTANTS[vss])/CONSTANTS[vmyo];
RATES[nass] = ( - (ALGEBRAIC[ICaNa]+ 3.00000*ALGEBRAIC[INaCa_ss])*CONSTANTS[cm]*CONSTANTS[Acap])/( CONSTANTS[F]*CONSTANTS[vss]) - ALGEBRAIC[JdiffNa];
#endif
#if defined PCA_ORUDY2011 
RATES[cai] =  ALGEBRAIC[Bcai]*((( - ((ALGEBRAIC[IpCa]+ALGEBRAIC[ICab]) -  2.00000*ALGEBRAIC[INaCa_i])*CONSTANTS[cm]*CONSTANTS[Acap])/( 2.00000*CONSTANTS[F]*CONSTANTS[vmyo]) - ( ALGEBRAIC[Jup]*CONSTANTS[vnsr])/CONSTANTS[vmyo])+( ALGEBRAIC[Jdiff]*CONSTANTS[vss])/CONSTANTS[vmyo]);
RATES[cansr] = ALGEBRAIC[Jup] - ( ALGEBRAIC[Jtr]*CONSTANTS[vjsr])/CONSTANTS[vnsr];
RATES[CaMKt] =  CONSTANTS[aCaMK]*ALGEBRAIC[CaMKb]*(ALGEBRAIC[CaMKb]+STATES[CaMKt]) -  CONSTANTS[bCaMK]*STATES[CaMKt];
RATES[cajsr] =  ALGEBRAIC[Bcajsr]*(ALGEBRAIC[Jtr] - ALGEBRAIC[Jrel]);
RATES[Jrelnp] = (ALGEBRAIC[Jrel_inf] - STATES[Jrelnp])/ALGEBRAIC[tau_rel];
RATES[Jrelp] = (ALGEBRAIC[Jrel_infp] - STATES[Jrelp])/ALGEBRAIC[tau_relp];
RATES[cass] =  ALGEBRAIC[Bcass]*((( - (ALGEBRAIC[ICaL] -  2.00000*ALGEBRAIC[INaCa_ss])*CONSTANTS[cm]*CONSTANTS[Acap])/( 2.00000*CONSTANTS[F]*CONSTANTS[vss])+( ALGEBRAIC[Jrel]*CONSTANTS[vjsr])/CONSTANTS[vss]) - ALGEBRAIC[Jdiff]);
#endif
#if defined TO_ORUDY2011 
RATES[a] = (ALGEBRAIC[ass] - STATES[a])/ALGEBRAIC[ta];
RATES[ap] = (ALGEBRAIC[assp] - STATES[ap])/ALGEBRAIC[ta];
RATES[ki] = ( - ((ALGEBRAIC[Ito]+ALGEBRAIC[IKr]+ALGEBRAIC[IKs]+ALGEBRAIC[IK1]+ALGEBRAIC[IKb]+ALGEBRAIC[Istim]) -  2.00000*ALGEBRAIC[INaK])*CONSTANTS[cm]*CONSTANTS[Acap])/( CONSTANTS[F]*CONSTANTS[vmyo])+( ALGEBRAIC[JdiffK]*CONSTANTS[vss])/CONSTANTS[vmyo];
RATES[iF] = (ALGEBRAIC[iss] - STATES[iF])/ALGEBRAIC[tiF];
RATES[iS] = (ALGEBRAIC[iss] - STATES[iS])/ALGEBRAIC[tiS];
RATES[kss] = ( - ALGEBRAIC[ICaK]*CONSTANTS[cm]*CONSTANTS[Acap])/( CONSTANTS[F]*CONSTANTS[vss]) - ALGEBRAIC[JdiffK];
RATES[iFp] = (ALGEBRAIC[iss] - STATES[iFp])/ALGEBRAIC[tiFp];
RATES[iSp] = (ALGEBRAIC[iss] - STATES[iSp])/ALGEBRAIC[tiSp];
RATES[CaMKt] =  CONSTANTS[aCaMK]*ALGEBRAIC[CaMKb]*(ALGEBRAIC[CaMKb]+STATES[CaMKt]) -  CONSTANTS[bCaMK]*STATES[CaMKt];
#endif
*/

RATES[43] = CONSTANTS[205];
RATES[34] = (- ( CONSTANTS[66]*exp( CONSTANTS[67]*STATES[0])*STATES[34]*exp(( (CONSTANTS[102] - 20.0000)*log(CONSTANTS[68]))/10.0000) -  CONSTANTS[69]*exp( CONSTANTS[70]*STATES[0])*STATES[35]*exp(( (CONSTANTS[102] - 20.0000)*log(CONSTANTS[71]))/10.0000))+ CONSTANTS[78]*exp( CONSTANTS[79]*STATES[0])*STATES[36]*exp(( (CONSTANTS[102] - 20.0000)*log(CONSTANTS[80]))/10.0000)) -  CONSTANTS[87]*exp( CONSTANTS[88]*STATES[0])*STATES[34]*exp(( (CONSTANTS[102] - 20.0000)*log(CONSTANTS[89]))/10.0000);
RATES[35] = ((( CONSTANTS[66]*exp( CONSTANTS[67]*STATES[0])*STATES[34]*exp(( (CONSTANTS[102] - 20.0000)*log(CONSTANTS[68]))/10.0000) -  CONSTANTS[69]*exp( CONSTANTS[70]*STATES[0])*STATES[35]*exp(( (CONSTANTS[102] - 20.0000)*log(CONSTANTS[71]))/10.0000)) - ( CONSTANTS[60]*exp( CONSTANTS[61]*STATES[0])*STATES[35]*exp(( (CONSTANTS[102] - 20.0000)*log(CONSTANTS[62]))/10.0000) -  CONSTANTS[63]*exp( CONSTANTS[64]*STATES[0])*STATES[39]*exp(( (CONSTANTS[102] - 20.0000)*log(CONSTANTS[65]))/10.0000)))+ CONSTANTS[81]*exp( CONSTANTS[82]*STATES[0])*STATES[37]*exp(( (CONSTANTS[102] - 20.0000)*log(CONSTANTS[83]))/10.0000)) -  CONSTANTS[90]*exp( CONSTANTS[91]*STATES[0])*STATES[35]*exp(( (CONSTANTS[102] - 20.0000)*log(CONSTANTS[92]))/10.0000);
RATES[36] = - ( CONSTANTS[54]*exp( CONSTANTS[55]*STATES[0])*STATES[36]*exp(( (CONSTANTS[102] - 20.0000)*log(CONSTANTS[56]))/10.0000) -  CONSTANTS[57]*exp( CONSTANTS[58]*STATES[0])*STATES[37]*exp(( (CONSTANTS[102] - 20.0000)*log(CONSTANTS[59]))/10.0000)) - ( CONSTANTS[78]*exp( CONSTANTS[79]*STATES[0])*STATES[36]*exp(( (CONSTANTS[102] - 20.0000)*log(CONSTANTS[80]))/10.0000) -  CONSTANTS[87]*exp( CONSTANTS[88]*STATES[0])*STATES[34]*exp(( (CONSTANTS[102] - 20.0000)*log(CONSTANTS[89]))/10.0000));
RATES[37] = (( CONSTANTS[54]*exp( CONSTANTS[55]*STATES[0])*STATES[36]*exp(( (CONSTANTS[102] - 20.0000)*log(CONSTANTS[56]))/10.0000) -  CONSTANTS[57]*exp( CONSTANTS[58]*STATES[0])*STATES[37]*exp(( (CONSTANTS[102] - 20.0000)*log(CONSTANTS[59]))/10.0000)) - ( CONSTANTS[72]*exp( CONSTANTS[73]*STATES[0])*STATES[37]*exp(( (CONSTANTS[102] - 20.0000)*log(CONSTANTS[74]))/10.0000) -  CONSTANTS[75]*exp( CONSTANTS[76]*STATES[0])*STATES[38]*exp(( (CONSTANTS[102] - 20.0000)*log(CONSTANTS[77]))/10.0000))) - ( CONSTANTS[81]*exp( CONSTANTS[82]*STATES[0])*STATES[37]*exp(( (CONSTANTS[102] - 20.0000)*log(CONSTANTS[83]))/10.0000) -  CONSTANTS[90]*exp( CONSTANTS[91]*STATES[0])*STATES[35]*exp(( (CONSTANTS[102] - 20.0000)*log(CONSTANTS[92]))/10.0000));
RATES[38] = (( CONSTANTS[72]*exp( CONSTANTS[73]*STATES[0])*STATES[37]*exp(( (CONSTANTS[102] - 20.0000)*log(CONSTANTS[74]))/10.0000) -  CONSTANTS[75]*exp( CONSTANTS[76]*STATES[0])*STATES[38]*exp(( (CONSTANTS[102] - 20.0000)*log(CONSTANTS[77]))/10.0000)) - ( CONSTANTS[84]*exp( CONSTANTS[85]*STATES[0])*STATES[38]*exp(( (CONSTANTS[102] - 20.0000)*log(CONSTANTS[86]))/10.0000) -  CONSTANTS[93]*exp( CONSTANTS[94]*STATES[0])*STATES[39]*exp(( (CONSTANTS[102] - 20.0000)*log(CONSTANTS[95]))/10.0000))) - ( (( CONSTANTS[96]*CONSTANTS[97]*exp( CONSTANTS[98]*log(STATES[43])))/(exp( CONSTANTS[98]*log(STATES[43]))+CONSTANTS[99]))*STATES[38] -  CONSTANTS[97]*STATES[41]);
RATES[39] = ((( CONSTANTS[60]*exp( CONSTANTS[61]*STATES[0])*STATES[35]*exp(( (CONSTANTS[102] - 20.0000)*log(CONSTANTS[62]))/10.0000) -  CONSTANTS[63]*exp( CONSTANTS[64]*STATES[0])*STATES[39]*exp(( (CONSTANTS[102] - 20.0000)*log(CONSTANTS[65]))/10.0000))+ CONSTANTS[84]*exp( CONSTANTS[85]*STATES[0])*STATES[38]*exp(( (CONSTANTS[102] - 20.0000)*log(CONSTANTS[86]))/10.0000)) -  CONSTANTS[93]*exp( CONSTANTS[94]*STATES[0])*STATES[39]*exp(( (CONSTANTS[102] - 20.0000)*log(CONSTANTS[95]))/10.0000)) - ( (( CONSTANTS[96]*CONSTANTS[97]*exp( CONSTANTS[98]*log(STATES[43])))/(exp( CONSTANTS[98]*log(STATES[43]))+CONSTANTS[99]))*STATES[39] -  (( CONSTANTS[97]*CONSTANTS[84]*exp( CONSTANTS[85]*STATES[0])*exp(( (CONSTANTS[102] - 20.0000)*log(CONSTANTS[86]))/10.0000))/( CONSTANTS[93]*exp( CONSTANTS[94]*STATES[0])*exp(( (CONSTANTS[102] - 20.0000)*log(CONSTANTS[95]))/10.0000)))*STATES[40]);
RATES[40] = (( (( CONSTANTS[96]*CONSTANTS[97]*exp( CONSTANTS[98]*log(STATES[43])))/(exp( CONSTANTS[98]*log(STATES[43]))+CONSTANTS[99]))*STATES[39] -  (( CONSTANTS[97]*CONSTANTS[84]*exp( CONSTANTS[85]*STATES[0])*exp(( (CONSTANTS[102] - 20.0000)*log(CONSTANTS[86]))/10.0000))/( CONSTANTS[93]*exp( CONSTANTS[94]*STATES[0])*exp(( (CONSTANTS[102] - 20.0000)*log(CONSTANTS[95]))/10.0000)))*STATES[40])+ (CONSTANTS[100]/(1.00000+exp(- (STATES[0] - CONSTANTS[101])/6.78900)))*STATES[42]) -  CONSTANTS[100]*STATES[40];
RATES[41] = (( (( CONSTANTS[96]*CONSTANTS[97]*exp( CONSTANTS[98]*log(STATES[43])))/(exp( CONSTANTS[98]*log(STATES[43]))+CONSTANTS[99]))*STATES[38] -  CONSTANTS[97]*STATES[41])+ (CONSTANTS[100]/(1.00000+exp(- (STATES[0] - CONSTANTS[101])/6.78900)))*STATES[42]) -  CONSTANTS[100]*STATES[41];
RATES[42] = - ( (CONSTANTS[100]/(1.00000+exp(- (STATES[0] - CONSTANTS[101])/6.78900)))*STATES[42] -  CONSTANTS[100]*STATES[41]) - ( (CONSTANTS[100]/(1.00000+exp(- (STATES[0] - CONSTANTS[101])/6.78900)))*STATES[42] -  CONSTANTS[100]*STATES[40]);
ALGEBRAIC[3] = 1.00000/(1.00000+exp((STATES[0]+87.6100)/7.48800));
RATES[17] = (ALGEBRAIC[3] - STATES[17])/CONSTANTS[47];
ALGEBRAIC[4] = 1.00000/(1.00000+exp((STATES[0]+93.8100)/7.48800));
RATES[18] = (ALGEBRAIC[4] - STATES[18])/CONSTANTS[152];
ALGEBRAIC[1] = 1.00000/(1.00000+exp(- (STATES[0]+CONSTANTS[34])/CONSTANTS[35]));
ALGEBRAIC[13] = 1.00000/( CONSTANTS[38]*exp((STATES[0]+CONSTANTS[36])/CONSTANTS[37])+ CONSTANTS[39]*exp(- (STATES[0]+CONSTANTS[40])/CONSTANTS[41]));
RATES[10] = (ALGEBRAIC[1] - STATES[10])/ALGEBRAIC[13];
ALGEBRAIC[2] = 1.00000/(1.00000+exp(((STATES[0]+CONSTANTS[42]) - CONSTANTS[46])/CONSTANTS[43]));
ALGEBRAIC[14] = 1.00000/( 1.43200e-05*exp(- ((STATES[0]+1.19600) - CONSTANTS[46])/6.28500)+ 6.14900*exp(((STATES[0]+0.509600) - CONSTANTS[46])/20.2700));
RATES[11] = (ALGEBRAIC[2] - STATES[11])/ALGEBRAIC[14];
ALGEBRAIC[15] = 1.00000/( 0.00979400*exp(- ((STATES[0]+17.9500) - CONSTANTS[46])/28.0500)+ 0.334300*exp(((STATES[0]+5.73000) - CONSTANTS[46])/56.6600));
RATES[12] = (ALGEBRAIC[2] - STATES[12])/ALGEBRAIC[15];
ALGEBRAIC[5] = 1.00000/(1.00000+exp(- (STATES[0] - 14.3400)/14.8200));
ALGEBRAIC[17] = 1.05150/(1.00000/( 1.20890*(1.00000+exp(- (STATES[0] - 18.4099)/29.3814)))+3.50000/(1.00000+exp((STATES[0]+100.000)/29.3814)));
RATES[19] = (ALGEBRAIC[5] - STATES[19])/ALGEBRAIC[17];
ALGEBRAIC[7] = 1.00000/(1.00000+exp(- (STATES[0]+3.94000)/4.23000));
ALGEBRAIC[21] = 0.600000+1.00000/(exp( - 0.0500000*(STATES[0]+6.00000))+exp( 0.0900000*(STATES[0]+14.0000)));
RATES[25] = (ALGEBRAIC[7] - STATES[25])/ALGEBRAIC[21];
ALGEBRAIC[8] = 1.00000/(1.00000+exp((STATES[0]+19.5800)/3.69600));
ALGEBRAIC[22] = 7.00000+1.00000/( 0.00450000*exp(- (STATES[0]+20.0000)/10.0000)+ 0.00450000*exp((STATES[0]+20.0000)/10.0000));
RATES[26] = (ALGEBRAIC[8] - STATES[26])/ALGEBRAIC[22];
ALGEBRAIC[23] = 1000.00+1.00000/( 3.50000e-05*exp(- (STATES[0]+5.00000)/4.00000)+ 3.50000e-05*exp((STATES[0]+5.00000)/6.00000));
RATES[27] = (ALGEBRAIC[8] - STATES[27])/ALGEBRAIC[23];
ALGEBRAIC[19] = ALGEBRAIC[8];
RATES[30] = (ALGEBRAIC[19] - STATES[30])/CONSTANTS[157];
ALGEBRAIC[9] =  STATES[30]*1.00000;
ALGEBRAIC[20] = 1.00000/(CONSTANTS[51]/ALGEBRAIC[9]+pow(1.00000+CONSTANTS[50]/STATES[2], 4.00000));
RATES[33] =  ALGEBRAIC[20]*CONSTANTS[51] -  STATES[33]*ALGEBRAIC[9];
ALGEBRAIC[10] = 1.00000/(1.00000+exp(- (STATES[0]+11.6000)/8.93200));
ALGEBRAIC[25] = CONSTANTS[104]+1.00000/( 0.000232600*exp((STATES[0]+48.2800)/17.8000)+ 0.00129200*exp(- (STATES[0]+210.000)/230.000));
RATES[44] = (ALGEBRAIC[10] - STATES[44])/ALGEBRAIC[25];
ALGEBRAIC[11] = 1.00000/(1.00000+exp(- (STATES[0]+ 2.55380*CONSTANTS[3]+144.590)/( 1.56920*CONSTANTS[3]+3.81150)));
ALGEBRAIC[26] = 122.200/(exp(- (STATES[0]+127.200)/20.3600)+exp((STATES[0]+236.800)/69.3300));
RATES[46] = (ALGEBRAIC[11] - STATES[46])/ALGEBRAIC[26];
ALGEBRAIC[36] = ( CONSTANTS[20]*(1.00000 - STATES[1]))/(1.00000+CONSTANTS[21]/STATES[2]);
RATES[1] =  CONSTANTS[18]*ALGEBRAIC[36]*(ALGEBRAIC[36]+STATES[1]) -  CONSTANTS[19]*STATES[1];
ALGEBRAIC[16] = ALGEBRAIC[2];
ALGEBRAIC[27] = 2.03800+1.00000/( 0.0213600*exp(- ((STATES[0]+100.600) - CONSTANTS[46])/8.28100)+ 0.305200*exp(((STATES[0]+0.994100) - CONSTANTS[46])/38.4500));
RATES[13] = (ALGEBRAIC[16] - STATES[13])/ALGEBRAIC[27];
ALGEBRAIC[31] = 1.00000/(1.00000+exp(- (STATES[0] - 24.3400)/14.8200));
RATES[22] = (ALGEBRAIC[31] - STATES[22])/ALGEBRAIC[17];
ALGEBRAIC[32] = 7.00000+1.00000/( 0.0400000*exp(- (STATES[0] - 4.00000)/7.00000)+ 0.0400000*exp((STATES[0] - 4.00000)/7.00000));
RATES[28] = (ALGEBRAIC[19] - STATES[28])/ALGEBRAIC[32];
ALGEBRAIC[33] = 100.000+1.00000/( 0.000120000*exp(- STATES[0]/3.00000)+ 0.000120000*exp(STATES[0]/7.00000));
RATES[29] = (ALGEBRAIC[19] - STATES[29])/ALGEBRAIC[33];
ALGEBRAIC[34] =  2.50000*ALGEBRAIC[22];
RATES[31] = (ALGEBRAIC[8] - STATES[31])/ALGEBRAIC[34];
ALGEBRAIC[24] = ALGEBRAIC[10];
ALGEBRAIC[35] = 1.00000/( 0.0100000*exp((STATES[0] - 50.0000)/20.0000)+ 0.0193000*exp(- (STATES[0]+66.5400)/31.0000));
RATES[45] = (ALGEBRAIC[24] - STATES[45])/ALGEBRAIC[35];
ALGEBRAIC[28] = 1.00000/(1.00000+exp(((STATES[0]+89.1000) - CONSTANTS[46])/6.08600));
ALGEBRAIC[37] =  3.00000*ALGEBRAIC[15];
RATES[14] = (ALGEBRAIC[28] - STATES[14])/ALGEBRAIC[37];
ALGEBRAIC[38] =  1.46000*ALGEBRAIC[27];
RATES[15] = (ALGEBRAIC[16] - STATES[15])/ALGEBRAIC[38];
ALGEBRAIC[29] = 1.00000/(1.00000+exp(- (STATES[0]+42.8500)/5.26400));
ALGEBRAIC[39] = ALGEBRAIC[13];
RATES[16] = (ALGEBRAIC[29] - STATES[16])/ALGEBRAIC[39];
ALGEBRAIC[41] =  2.50000*ALGEBRAIC[32];
RATES[32] = (ALGEBRAIC[19] - STATES[32])/ALGEBRAIC[41];
ALGEBRAIC[6] = 1.00000/(1.00000+exp((STATES[0]+43.9400)/5.71100));
ALGEBRAIC[18] = (CONSTANTS[0]==1.00000 ? 1.00000 - 0.950000/(1.00000+exp((STATES[0]+70.0000)/5.00000)) : 1.00000);
ALGEBRAIC[30] = 4.56200+1.00000/( 0.393300*exp(- (STATES[0]+100.000)/100.000)+ 0.0800400*exp((STATES[0]+50.0000)/16.5900));
ALGEBRAIC[43] =  ALGEBRAIC[30]*ALGEBRAIC[18];
RATES[20] = (ALGEBRAIC[6] - STATES[20])/ALGEBRAIC[43];
ALGEBRAIC[40] = 23.6200+1.00000/( 0.00141600*exp(- (STATES[0]+96.5200)/59.0500)+ 1.78000e-08*exp((STATES[0]+114.100)/8.07900));
ALGEBRAIC[45] =  ALGEBRAIC[40]*ALGEBRAIC[18];
RATES[21] = (ALGEBRAIC[6] - STATES[21])/ALGEBRAIC[45];
ALGEBRAIC[47] = 1.35400+0.000100000/(exp((STATES[0] - 167.400)/15.8900)+exp(- (STATES[0] - 12.2300)/0.215400));
ALGEBRAIC[49] = 1.00000 - 0.500000/(1.00000+exp((STATES[0]+70.0000)/20.0000));
ALGEBRAIC[51] =  ALGEBRAIC[47]*ALGEBRAIC[49]*ALGEBRAIC[43];
RATES[23] = (ALGEBRAIC[6] - STATES[23])/ALGEBRAIC[51];
ALGEBRAIC[52] =  ALGEBRAIC[47]*ALGEBRAIC[49]*ALGEBRAIC[45];
RATES[24] = (ALGEBRAIC[6] - STATES[24])/ALGEBRAIC[52];
ALGEBRAIC[67] =  CONSTANTS[155]*STATES[26]+ CONSTANTS[170]*STATES[27];
ALGEBRAIC[68] = 0.300000+0.600000/(1.00000+exp((STATES[0] - 10.0000)/10.0000));
ALGEBRAIC[69] = 1.00000 - ALGEBRAIC[68];
ALGEBRAIC[70] =  ALGEBRAIC[68]*STATES[28]+ ALGEBRAIC[69]*STATES[29];
ALGEBRAIC[71] =  CONSTANTS[155]*STATES[31]+ CONSTANTS[170]*STATES[27];
ALGEBRAIC[72] =  ALGEBRAIC[68]*STATES[32]+ ALGEBRAIC[69]*STATES[29];
ALGEBRAIC[12] =  STATES[0]*CONSTANTS[149];
ALGEBRAIC[73] = ( 4.00000*CONSTANTS[169]*( STATES[2]*exp( 2.00000*ALGEBRAIC[12]) -  0.341000*CONSTANTS[2]))/CONSTANTS[174];
ALGEBRAIC[74] =  CONSTANTS[174]*(STATES[0] - CONSTANTS[158]);
ALGEBRAIC[75] = (- 1.00000e-07<=ALGEBRAIC[74]&&ALGEBRAIC[74]<=1.00000e-07 ?  ALGEBRAIC[73]*(1.00000 -  0.500000*ALGEBRAIC[74]) : ( ALGEBRAIC[73]*ALGEBRAIC[74])/(exp(ALGEBRAIC[74]) - 1.00000));
ALGEBRAIC[42] = ALGEBRAIC[36]+STATES[1];
ALGEBRAIC[82] = 1.00000/(1.00000+CONSTANTS[17]/ALGEBRAIC[42]);
ALGEBRAIC[83] =  (1.00000 - ALGEBRAIC[82])*CONSTANTS[156]*ALGEBRAIC[75]*STATES[25]*( ALGEBRAIC[67]*(1.00000 - STATES[33])+ STATES[30]*ALGEBRAIC[70]*STATES[33])+ ALGEBRAIC[82]*CONSTANTS[171]*ALGEBRAIC[75]*STATES[25]*( ALGEBRAIC[71]*(1.00000 - STATES[33])+ STATES[30]*ALGEBRAIC[72]*STATES[33]);
ALGEBRAIC[85] = ( CONSTANTS[166]*- ALGEBRAIC[83])/(1.00000+ 1.00000*pow(1.50000/STATES[8], 8.00000));
ALGEBRAIC[88] = (CONSTANTS[0]==2.00000 ?  ALGEBRAIC[85]*1.70000 : ALGEBRAIC[85]);
ALGEBRAIC[91] = CONSTANTS[146]/(1.00000+0.0123000/STATES[8]);
ALGEBRAIC[94] = (ALGEBRAIC[91]<0.00100000 ? 0.00100000 : ALGEBRAIC[91]);
RATES[47] = (ALGEBRAIC[88] - STATES[47])/ALGEBRAIC[94];
ALGEBRAIC[86] = ( CONSTANTS[180]*- ALGEBRAIC[83])/(1.00000+pow(1.50000/STATES[8], 8.00000));
ALGEBRAIC[89] = (CONSTANTS[0]==2.00000 ?  ALGEBRAIC[86]*1.70000 : ALGEBRAIC[86]);
ALGEBRAIC[92] = CONSTANTS[167]/(1.00000+0.0123000/STATES[8]);
ALGEBRAIC[95] = (ALGEBRAIC[92]<0.00100000 ? 0.00100000 : ALGEBRAIC[92]);
RATES[48] = (ALGEBRAIC[89] - STATES[48])/ALGEBRAIC[95];
ALGEBRAIC[53] =  (( CONSTANTS[4]*CONSTANTS[5])/CONSTANTS[6])*log(CONSTANTS[3]/STATES[5]);
ALGEBRAIC[61] = 1.00000/(1.00000+exp((STATES[0] - 213.600)/151.200));
ALGEBRAIC[62] = 1.00000 - ALGEBRAIC[61];
ALGEBRAIC[63] =  ALGEBRAIC[61]*STATES[20]+ ALGEBRAIC[62]*STATES[21];
ALGEBRAIC[64] =  ALGEBRAIC[61]*STATES[23]+ ALGEBRAIC[62]*STATES[24];
ALGEBRAIC[65] = 1.00000/(1.00000+CONSTANTS[17]/ALGEBRAIC[42]);
ALGEBRAIC[66] =  CONSTANTS[154]*(STATES[0] - ALGEBRAIC[53])*( (1.00000 - ALGEBRAIC[65])*STATES[19]*ALGEBRAIC[63]+ ALGEBRAIC[65]*STATES[22]*ALGEBRAIC[64]);
ALGEBRAIC[90] =  CONSTANTS[159]* pow((CONSTANTS[3]/5.40000), 1.0 / 2)*STATES[38]*(STATES[0] - ALGEBRAIC[53]);
ALGEBRAIC[54] =  (( CONSTANTS[4]*CONSTANTS[5])/CONSTANTS[6])*log((CONSTANTS[3]+ CONSTANTS[33]*CONSTANTS[1])/(STATES[5]+ CONSTANTS[33]*STATES[3]));
ALGEBRAIC[93] = 1.00000+0.600000/(1.00000+pow(3.80000e-05/STATES[9], 1.40000));
ALGEBRAIC[96] =  CONSTANTS[160]*ALGEBRAIC[93]*STATES[44]*STATES[45]*(STATES[0] - ALGEBRAIC[54]);
ALGEBRAIC[97] = 1.00000/(1.00000+exp(((STATES[0]+105.800) -  2.60000*CONSTANTS[3])/9.49300));
ALGEBRAIC[98] =  CONSTANTS[161]* pow(CONSTANTS[3], 1.0 / 2)*ALGEBRAIC[97]*STATES[46]*(STATES[0] - ALGEBRAIC[53]);
ALGEBRAIC[162] =  CONSTANTS[128]*exp(( (1.00000 - CONSTANTS[129])*STATES[0]*CONSTANTS[6])/( 3.00000*CONSTANTS[4]*CONSTANTS[5]));
ALGEBRAIC[166] = ( CONSTANTS[123]*pow(CONSTANTS[3]/CONSTANTS[131], 2.00000))/((pow(1.00000+CONSTANTS[1]/ALGEBRAIC[162], 3.00000)+pow(1.00000+CONSTANTS[3]/CONSTANTS[131], 2.00000)) - 1.00000);
ALGEBRAIC[163] = CONSTANTS[136]/(1.00000+CONSTANTS[135]/CONSTANTS[137]+STATES[3]/CONSTANTS[138]+STATES[5]/CONSTANTS[139]);
ALGEBRAIC[167] = ( CONSTANTS[124]*ALGEBRAIC[163]*CONSTANTS[135])/(1.00000+CONSTANTS[133]/CONSTANTS[134]);
ALGEBRAIC[161] =  CONSTANTS[127]*exp(( CONSTANTS[129]*STATES[0]*CONSTANTS[6])/( 3.00000*CONSTANTS[4]*CONSTANTS[5]));
ALGEBRAIC[164] = ( CONSTANTS[119]*pow(STATES[3]/ALGEBRAIC[161], 3.00000))/((pow(1.00000+STATES[3]/ALGEBRAIC[161], 3.00000)+pow(1.00000+STATES[5]/CONSTANTS[130], 2.00000)) - 1.00000);
ALGEBRAIC[165] = ( CONSTANTS[122]*pow(CONSTANTS[1]/ALGEBRAIC[162], 3.00000))/((pow(1.00000+CONSTANTS[1]/ALGEBRAIC[162], 3.00000)+pow(1.00000+CONSTANTS[3]/CONSTANTS[131], 2.00000)) - 1.00000);
ALGEBRAIC[168] = ( CONSTANTS[126]*pow(STATES[5]/CONSTANTS[130], 2.00000))/((pow(1.00000+STATES[3]/ALGEBRAIC[161], 3.00000)+pow(1.00000+STATES[5]/CONSTANTS[130], 2.00000)) - 1.00000);
ALGEBRAIC[169] =  CONSTANTS[203]*ALGEBRAIC[164]*CONSTANTS[202]+ ALGEBRAIC[165]*ALGEBRAIC[168]*ALGEBRAIC[167]+ CONSTANTS[202]*ALGEBRAIC[168]*ALGEBRAIC[167]+ ALGEBRAIC[167]*ALGEBRAIC[164]*CONSTANTS[202];
ALGEBRAIC[170] =  ALGEBRAIC[165]*CONSTANTS[201]*ALGEBRAIC[168]+ ALGEBRAIC[164]*CONSTANTS[202]*ALGEBRAIC[166]+ ALGEBRAIC[166]*CONSTANTS[201]*ALGEBRAIC[168]+ CONSTANTS[202]*ALGEBRAIC[166]*ALGEBRAIC[168];
ALGEBRAIC[171] =  CONSTANTS[202]*ALGEBRAIC[166]*CONSTANTS[203]+ ALGEBRAIC[167]*ALGEBRAIC[165]*CONSTANTS[201]+ ALGEBRAIC[165]*CONSTANTS[201]*CONSTANTS[203]+ ALGEBRAIC[166]*CONSTANTS[203]*CONSTANTS[201];
ALGEBRAIC[172] =  ALGEBRAIC[168]*ALGEBRAIC[167]*ALGEBRAIC[165]+ ALGEBRAIC[166]*CONSTANTS[203]*ALGEBRAIC[164]+ ALGEBRAIC[165]*CONSTANTS[203]*ALGEBRAIC[164]+ ALGEBRAIC[167]*ALGEBRAIC[165]*ALGEBRAIC[164];
ALGEBRAIC[173] = ALGEBRAIC[169]/(ALGEBRAIC[169]+ALGEBRAIC[170]+ALGEBRAIC[171]+ALGEBRAIC[172]);
ALGEBRAIC[174] = ALGEBRAIC[170]/(ALGEBRAIC[169]+ALGEBRAIC[170]+ALGEBRAIC[171]+ALGEBRAIC[172]);
ALGEBRAIC[177] =  3.00000*( ALGEBRAIC[173]*ALGEBRAIC[166] -  ALGEBRAIC[174]*ALGEBRAIC[167]);
ALGEBRAIC[175] = ALGEBRAIC[171]/(ALGEBRAIC[169]+ALGEBRAIC[170]+ALGEBRAIC[171]+ALGEBRAIC[172]);
ALGEBRAIC[176] = ALGEBRAIC[172]/(ALGEBRAIC[169]+ALGEBRAIC[170]+ALGEBRAIC[171]+ALGEBRAIC[172]);
ALGEBRAIC[178] =  2.00000*( ALGEBRAIC[176]*CONSTANTS[201] -  ALGEBRAIC[175]*ALGEBRAIC[164]);
ALGEBRAIC[179] =  CONSTANTS[204]*( CONSTANTS[7]*ALGEBRAIC[177]+ CONSTANTS[9]*ALGEBRAIC[178]);
ALGEBRAIC[180] = 1.00000/(1.00000+exp(- (STATES[0] - 14.4800)/18.3400));
ALGEBRAIC[181] =  CONSTANTS[163]*ALGEBRAIC[180]*(STATES[0] - ALGEBRAIC[53]);
ALGEBRAIC[0] = (TIME>=CONSTANTS[12]&&TIME<=CONSTANTS[13]&&(TIME - CONSTANTS[12]) -  floor((TIME - CONSTANTS[12])/CONSTANTS[15])*CONSTANTS[15]<=CONSTANTS[16] ? CONSTANTS[14] : 0.000000);
ALGEBRAIC[183] = (STATES[6] - STATES[5])/2.00000;
RATES[5] = ( - ((ALGEBRAIC[66]+ALGEBRAIC[90]+ALGEBRAIC[96]+ALGEBRAIC[98]+ALGEBRAIC[181]+ALGEBRAIC[0]) -  2.00000*ALGEBRAIC[179])*CONSTANTS[32]*CONSTANTS[183])/( CONSTANTS[6]*CONSTANTS[184])+( ALGEBRAIC[183]*CONSTANTS[187])/CONSTANTS[184];
ALGEBRAIC[79] = ( 0.750000*CONSTANTS[169]*( STATES[6]*exp(ALGEBRAIC[12]) - CONSTANTS[3]))/CONSTANTS[176];
ALGEBRAIC[80] =  CONSTANTS[176]*(STATES[0] - CONSTANTS[158]);
ALGEBRAIC[81] = (- 1.00000e-07<=ALGEBRAIC[80]&&ALGEBRAIC[80]<=1.00000e-07 ?  ALGEBRAIC[79]*(1.00000 -  0.500000*ALGEBRAIC[80]) : ( ALGEBRAIC[79]*ALGEBRAIC[80])/(exp(ALGEBRAIC[80]) - 1.00000));
ALGEBRAIC[87] =  (1.00000 - ALGEBRAIC[82])*CONSTANTS[173]*ALGEBRAIC[81]*STATES[25]*( ALGEBRAIC[67]*(1.00000 - STATES[33])+ STATES[30]*ALGEBRAIC[70]*STATES[33])+ ALGEBRAIC[82]*CONSTANTS[182]*ALGEBRAIC[81]*STATES[25]*( ALGEBRAIC[71]*(1.00000 - STATES[33])+ STATES[30]*ALGEBRAIC[72]*STATES[33]);
RATES[6] = ( - ALGEBRAIC[87]*CONSTANTS[32]*CONSTANTS[183])/( CONSTANTS[6]*CONSTANTS[187]) - ALGEBRAIC[183];
ALGEBRAIC[50] =  (( CONSTANTS[4]*CONSTANTS[5])/CONSTANTS[6])*log(CONSTANTS[1]/STATES[3]);
ALGEBRAIC[55] =  CONSTANTS[44]*STATES[11]+ CONSTANTS[151]*STATES[12];
ALGEBRAIC[56] =  CONSTANTS[44]*STATES[11]+ CONSTANTS[151]*STATES[14];
ALGEBRAIC[57] = 1.00000/(1.00000+CONSTANTS[17]/ALGEBRAIC[42]);
ALGEBRAIC[58] =  CONSTANTS[45]*(STATES[0] - ALGEBRAIC[50])*pow(STATES[10], 3.00000)*( (1.00000 - ALGEBRAIC[57])*ALGEBRAIC[55]*STATES[13]+ ALGEBRAIC[57]*ALGEBRAIC[56]*STATES[15]);
ALGEBRAIC[59] = 1.00000/(1.00000+CONSTANTS[17]/ALGEBRAIC[42]);
ALGEBRAIC[60] =  CONSTANTS[153]*(STATES[0] - ALGEBRAIC[50])*STATES[16]*( (1.00000 - ALGEBRAIC[59])*STATES[17]+ ALGEBRAIC[59]*STATES[18]);
ALGEBRAIC[127] = 1.00000/(1.00000+pow(CONSTANTS[117]/STATES[9], 2.00000));
ALGEBRAIC[100] = exp(( CONSTANTS[115]*STATES[0]*CONSTANTS[6])/( CONSTANTS[4]*CONSTANTS[5]));
ALGEBRAIC[107] = 1.00000+ (CONSTANTS[1]/CONSTANTS[108])*(1.00000+1.00000/ALGEBRAIC[100]);
ALGEBRAIC[108] = CONSTANTS[1]/( CONSTANTS[108]*ALGEBRAIC[100]*ALGEBRAIC[107]);
ALGEBRAIC[111] =  ALGEBRAIC[108]*CONSTANTS[112];
ALGEBRAIC[101] = 1.00000+ (STATES[3]/CONSTANTS[108])*(1.00000+ALGEBRAIC[100]);
ALGEBRAIC[102] = ( STATES[3]*ALGEBRAIC[100])/( CONSTANTS[108]*ALGEBRAIC[101]);
ALGEBRAIC[114] =  ALGEBRAIC[102]*CONSTANTS[112];
ALGEBRAIC[104] = 1.00000+ (STATES[3]/CONSTANTS[106])*(1.00000+STATES[3]/CONSTANTS[107]);
ALGEBRAIC[105] = ( STATES[3]*STATES[3])/( ALGEBRAIC[104]*CONSTANTS[106]*CONSTANTS[107]);
ALGEBRAIC[117] =  ALGEBRAIC[105]*ALGEBRAIC[102]*CONSTANTS[110];
ALGEBRAIC[118] =  ALGEBRAIC[108]*CONSTANTS[189]*CONSTANTS[110];
ALGEBRAIC[109] = 1.00000/ALGEBRAIC[107];
ALGEBRAIC[110] =  ALGEBRAIC[109]*CONSTANTS[111];
ALGEBRAIC[112] = ALGEBRAIC[110]+ALGEBRAIC[111];
ALGEBRAIC[99] = exp(( CONSTANTS[116]*STATES[0]*CONSTANTS[6])/( CONSTANTS[4]*CONSTANTS[5]));
ALGEBRAIC[103] = 1.00000/ALGEBRAIC[101];
ALGEBRAIC[113] = ( ALGEBRAIC[103]*CONSTANTS[111])/ALGEBRAIC[99];
ALGEBRAIC[115] = ALGEBRAIC[113]+ALGEBRAIC[114];
ALGEBRAIC[106] = 1.00000/ALGEBRAIC[104];
ALGEBRAIC[116] =  ALGEBRAIC[106]*STATES[9]*CONSTANTS[113];
ALGEBRAIC[119] =  CONSTANTS[192]*ALGEBRAIC[115]*(ALGEBRAIC[117]+ALGEBRAIC[116])+ CONSTANTS[193]*ALGEBRAIC[117]*(CONSTANTS[192]+ALGEBRAIC[112]);
ALGEBRAIC[120] =  CONSTANTS[191]*ALGEBRAIC[117]*(ALGEBRAIC[115]+CONSTANTS[193])+ ALGEBRAIC[115]*ALGEBRAIC[116]*(CONSTANTS[191]+ALGEBRAIC[118]);
ALGEBRAIC[121] =  CONSTANTS[191]*ALGEBRAIC[112]*(ALGEBRAIC[117]+ALGEBRAIC[116])+ ALGEBRAIC[118]*ALGEBRAIC[116]*(CONSTANTS[192]+ALGEBRAIC[112]);
ALGEBRAIC[122] =  CONSTANTS[192]*ALGEBRAIC[118]*(ALGEBRAIC[115]+CONSTANTS[193])+ ALGEBRAIC[112]*CONSTANTS[193]*(CONSTANTS[191]+ALGEBRAIC[118]);
ALGEBRAIC[123] = ALGEBRAIC[119]/(ALGEBRAIC[119]+ALGEBRAIC[120]+ALGEBRAIC[121]+ALGEBRAIC[122]);
ALGEBRAIC[124] = ALGEBRAIC[120]/(ALGEBRAIC[119]+ALGEBRAIC[120]+ALGEBRAIC[121]+ALGEBRAIC[122]);
ALGEBRAIC[125] = ALGEBRAIC[121]/(ALGEBRAIC[119]+ALGEBRAIC[120]+ALGEBRAIC[121]+ALGEBRAIC[122]);
ALGEBRAIC[126] = ALGEBRAIC[122]/(ALGEBRAIC[119]+ALGEBRAIC[120]+ALGEBRAIC[121]+ALGEBRAIC[122]);
ALGEBRAIC[128] = ( 3.00000*( ALGEBRAIC[126]*ALGEBRAIC[117] -  ALGEBRAIC[123]*ALGEBRAIC[118])+ ALGEBRAIC[125]*ALGEBRAIC[114]) -  ALGEBRAIC[124]*ALGEBRAIC[111];
ALGEBRAIC[129] =  ALGEBRAIC[124]*CONSTANTS[192] -  ALGEBRAIC[123]*CONSTANTS[191];
ALGEBRAIC[130] =  0.800000*CONSTANTS[194]*ALGEBRAIC[127]*( CONSTANTS[7]*ALGEBRAIC[128]+ CONSTANTS[8]*ALGEBRAIC[129]);
ALGEBRAIC[182] = ( CONSTANTS[142]*CONSTANTS[169]*( STATES[3]*exp(ALGEBRAIC[12]) - CONSTANTS[1]))/CONSTANTS[178];
ALGEBRAIC[184] =  CONSTANTS[178]*(STATES[0] - CONSTANTS[164]);
ALGEBRAIC[185] = (- 1.00000e-07<=ALGEBRAIC[184]&&ALGEBRAIC[184]<=1.00000e-07 ?  ALGEBRAIC[182]*(1.00000 -  0.500000*ALGEBRAIC[184]) : ( ALGEBRAIC[182]*ALGEBRAIC[184])/(exp(ALGEBRAIC[184]) - 1.00000));
ALGEBRAIC[187] = (STATES[4] - STATES[3])/2.00000;
RATES[3] = ( - (ALGEBRAIC[58]+ALGEBRAIC[60]+ 3.00000*ALGEBRAIC[130]+ 3.00000*ALGEBRAIC[179]+ALGEBRAIC[185])*CONSTANTS[183]*CONSTANTS[32])/( CONSTANTS[6]*CONSTANTS[184])+( ALGEBRAIC[187]*CONSTANTS[187])/CONSTANTS[184];
ALGEBRAIC[76] = ( 0.750000*CONSTANTS[169]*( STATES[4]*exp(ALGEBRAIC[12]) - CONSTANTS[1]))/CONSTANTS[175];
ALGEBRAIC[77] =  CONSTANTS[175]*(STATES[0] - CONSTANTS[158]);
ALGEBRAIC[78] = (- 1.00000e-07<=ALGEBRAIC[77]&&ALGEBRAIC[77]<=1.00000e-07 ?  ALGEBRAIC[76]*(1.00000 -  0.500000*ALGEBRAIC[77]) : ( ALGEBRAIC[76]*ALGEBRAIC[77])/(exp(ALGEBRAIC[77]) - 1.00000));
ALGEBRAIC[84] =  (1.00000 - ALGEBRAIC[82])*CONSTANTS[172]*ALGEBRAIC[78]*STATES[25]*( ALGEBRAIC[67]*(1.00000 - STATES[33])+ STATES[30]*ALGEBRAIC[70]*STATES[33])+ ALGEBRAIC[82]*CONSTANTS[181]*ALGEBRAIC[78]*STATES[25]*( ALGEBRAIC[71]*(1.00000 - STATES[33])+ STATES[30]*ALGEBRAIC[72]*STATES[33]);
ALGEBRAIC[157] = 1.00000/(1.00000+pow(CONSTANTS[117]/STATES[2], 2.00000));
ALGEBRAIC[137] = 1.00000+ (CONSTANTS[1]/CONSTANTS[108])*(1.00000+1.00000/ALGEBRAIC[100]);
ALGEBRAIC[138] = CONSTANTS[1]/( CONSTANTS[108]*ALGEBRAIC[100]*ALGEBRAIC[137]);
ALGEBRAIC[141] =  ALGEBRAIC[138]*CONSTANTS[112];
ALGEBRAIC[131] = 1.00000+ (STATES[4]/CONSTANTS[108])*(1.00000+ALGEBRAIC[100]);
ALGEBRAIC[132] = ( STATES[4]*ALGEBRAIC[100])/( CONSTANTS[108]*ALGEBRAIC[131]);
ALGEBRAIC[144] =  ALGEBRAIC[132]*CONSTANTS[112];
ALGEBRAIC[134] = 1.00000+ (STATES[4]/CONSTANTS[106])*(1.00000+STATES[4]/CONSTANTS[107]);
ALGEBRAIC[135] = ( STATES[4]*STATES[4])/( ALGEBRAIC[134]*CONSTANTS[106]*CONSTANTS[107]);
ALGEBRAIC[147] =  ALGEBRAIC[135]*ALGEBRAIC[132]*CONSTANTS[110];
ALGEBRAIC[148] =  ALGEBRAIC[138]*CONSTANTS[196]*CONSTANTS[110];
ALGEBRAIC[139] = 1.00000/ALGEBRAIC[137];
ALGEBRAIC[140] =  ALGEBRAIC[139]*CONSTANTS[111];
ALGEBRAIC[142] = ALGEBRAIC[140]+ALGEBRAIC[141];
ALGEBRAIC[133] = 1.00000/ALGEBRAIC[131];
ALGEBRAIC[143] = ( ALGEBRAIC[133]*CONSTANTS[111])/ALGEBRAIC[99];
ALGEBRAIC[145] = ALGEBRAIC[143]+ALGEBRAIC[144];
ALGEBRAIC[136] = 1.00000/ALGEBRAIC[134];
ALGEBRAIC[146] =  ALGEBRAIC[136]*STATES[2]*CONSTANTS[113];
ALGEBRAIC[149] =  CONSTANTS[199]*ALGEBRAIC[145]*(ALGEBRAIC[147]+ALGEBRAIC[146])+ CONSTANTS[200]*ALGEBRAIC[147]*(CONSTANTS[199]+ALGEBRAIC[142]);
ALGEBRAIC[150] =  CONSTANTS[198]*ALGEBRAIC[147]*(ALGEBRAIC[145]+CONSTANTS[200])+ ALGEBRAIC[145]*ALGEBRAIC[146]*(CONSTANTS[198]+ALGEBRAIC[148]);
ALGEBRAIC[151] =  CONSTANTS[198]*ALGEBRAIC[142]*(ALGEBRAIC[147]+ALGEBRAIC[146])+ ALGEBRAIC[148]*ALGEBRAIC[146]*(CONSTANTS[199]+ALGEBRAIC[142]);
ALGEBRAIC[152] =  CONSTANTS[199]*ALGEBRAIC[148]*(ALGEBRAIC[145]+CONSTANTS[200])+ ALGEBRAIC[142]*CONSTANTS[200]*(CONSTANTS[198]+ALGEBRAIC[148]);
ALGEBRAIC[153] = ALGEBRAIC[149]/(ALGEBRAIC[149]+ALGEBRAIC[150]+ALGEBRAIC[151]+ALGEBRAIC[152]);
ALGEBRAIC[154] = ALGEBRAIC[150]/(ALGEBRAIC[149]+ALGEBRAIC[150]+ALGEBRAIC[151]+ALGEBRAIC[152]);
ALGEBRAIC[155] = ALGEBRAIC[151]/(ALGEBRAIC[149]+ALGEBRAIC[150]+ALGEBRAIC[151]+ALGEBRAIC[152]);
ALGEBRAIC[156] = ALGEBRAIC[152]/(ALGEBRAIC[149]+ALGEBRAIC[150]+ALGEBRAIC[151]+ALGEBRAIC[152]);
ALGEBRAIC[158] = ( 3.00000*( ALGEBRAIC[156]*ALGEBRAIC[147] -  ALGEBRAIC[153]*ALGEBRAIC[148])+ ALGEBRAIC[155]*ALGEBRAIC[144]) -  ALGEBRAIC[154]*ALGEBRAIC[141];
ALGEBRAIC[159] =  ALGEBRAIC[154]*CONSTANTS[199] -  ALGEBRAIC[153]*CONSTANTS[198];
ALGEBRAIC[160] =  0.200000*CONSTANTS[194]*ALGEBRAIC[157]*( CONSTANTS[7]*ALGEBRAIC[158]+ CONSTANTS[8]*ALGEBRAIC[159]);
RATES[4] = ( - (ALGEBRAIC[84]+ 3.00000*ALGEBRAIC[160])*CONSTANTS[32]*CONSTANTS[183])/( CONSTANTS[6]*CONSTANTS[187]) - ALGEBRAIC[187];
ALGEBRAIC[190] = ( CONSTANTS[144]*STATES[9])/(CONSTANTS[145]+STATES[9]);
ALGEBRAIC[186] = ( CONSTANTS[143]*4.00000*CONSTANTS[169]*( STATES[9]*exp( 2.00000*ALGEBRAIC[12]) -  0.341000*CONSTANTS[2]))/CONSTANTS[179];
ALGEBRAIC[188] =  CONSTANTS[179]*(STATES[0] - CONSTANTS[165]);
ALGEBRAIC[189] = (- 1.00000e-07<=ALGEBRAIC[188]&&ALGEBRAIC[188]<=1.00000e-07 ?  ALGEBRAIC[186]*(1.00000 -  0.500000*ALGEBRAIC[188]) : ( ALGEBRAIC[186]*ALGEBRAIC[188])/(exp(ALGEBRAIC[188]) - 1.00000));
RATES[0] = - (ALGEBRAIC[58]+ALGEBRAIC[60]+ALGEBRAIC[66]+ALGEBRAIC[83]+ALGEBRAIC[84]+ALGEBRAIC[87]+ALGEBRAIC[90]+ALGEBRAIC[96]+ALGEBRAIC[98]+ALGEBRAIC[130]+ALGEBRAIC[160]+ALGEBRAIC[179]+ALGEBRAIC[185]+ALGEBRAIC[181]+ALGEBRAIC[190]+ALGEBRAIC[189]+ALGEBRAIC[0]);
ALGEBRAIC[191] = (STATES[2] - STATES[9])/0.200000;
ALGEBRAIC[192] = 1.00000/(1.00000+CONSTANTS[17]/ALGEBRAIC[42]);
ALGEBRAIC[193] =  CONSTANTS[147]*( (1.00000 - ALGEBRAIC[192])*STATES[47]+ ALGEBRAIC[192]*STATES[48]);
ALGEBRAIC[46] = 1.00000/(1.00000+( CONSTANTS[26]*CONSTANTS[27])/pow(CONSTANTS[27]+STATES[2], 2.00000)+( CONSTANTS[28]*CONSTANTS[29])/pow(CONSTANTS[29]+STATES[2], 2.00000));
RATES[2] =  ALGEBRAIC[46]*((( - (ALGEBRAIC[83] -  2.00000*ALGEBRAIC[160])*CONSTANTS[32]*CONSTANTS[183])/( 2.00000*CONSTANTS[6]*CONSTANTS[187])+( ALGEBRAIC[193]*CONSTANTS[186])/CONSTANTS[187]) - ALGEBRAIC[191]);
ALGEBRAIC[194] = ( CONSTANTS[168]*0.00437500*STATES[9])/(STATES[9]+0.000920000);
ALGEBRAIC[195] = ( CONSTANTS[168]*2.75000*0.00437500*STATES[9])/((STATES[9]+0.000920000) - 0.000170000);
ALGEBRAIC[196] = 1.00000/(1.00000+CONSTANTS[17]/ALGEBRAIC[42]);
ALGEBRAIC[197] = ( 0.00393750*STATES[7])/15.0000;
ALGEBRAIC[198] =  CONSTANTS[148]*(( (1.00000 - ALGEBRAIC[196])*ALGEBRAIC[194]+ ALGEBRAIC[196]*ALGEBRAIC[195]) - ALGEBRAIC[197]);
ALGEBRAIC[44] = 1.00000/(1.00000+( CONSTANTS[150]*CONSTANTS[23])/pow(CONSTANTS[23]+STATES[9], 2.00000)+( CONSTANTS[24]*CONSTANTS[25])/pow(CONSTANTS[25]+STATES[9], 2.00000));
RATES[9] =  ALGEBRAIC[44]*((( - ((ALGEBRAIC[190]+ALGEBRAIC[189]) -  2.00000*ALGEBRAIC[130])*CONSTANTS[32]*CONSTANTS[183])/( 2.00000*CONSTANTS[6]*CONSTANTS[184]) - ( ALGEBRAIC[198]*CONSTANTS[185])/CONSTANTS[184])+( ALGEBRAIC[191]*CONSTANTS[187])/CONSTANTS[184]);
ALGEBRAIC[199] = (STATES[7] - STATES[8])/100.000;
RATES[7] = ALGEBRAIC[198] - ( ALGEBRAIC[199]*CONSTANTS[186])/CONSTANTS[185];
ALGEBRAIC[48] = 1.00000/(1.00000+( CONSTANTS[30]*CONSTANTS[31])/pow(CONSTANTS[31]+STATES[8], 2.00000));
RATES[8] =  ALGEBRAIC[48]*(ALGEBRAIC[199] - ALGEBRAIC[193]);
}
