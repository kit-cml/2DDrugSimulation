/*
   There are a total of 198 entries in the algebraic variable array.
   There are a total of 41+3 entries in each of the rate and state variable arrays.
   There are a total of 139+3 entries in the constant variable array.
 */

#include "Ohara_Rudy_2011.hpp"
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
 * ALGEBRAIC[vffrt] is vffrt in component membrane (coulomb_per_mole).
 * ALGEBRAIC[vfrt] is vfrt in component membrane (dimensionless).
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
 * CONSTANTS[amp] is amp in component membrane (microA_per_microF).
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
 * CONSTANTS[GKr_b] is GKr_b in component IKr (milliS_per_microF).
 * CONSTANTS[GKr] is GKr in component IKr (milliS_per_microF).
 * ALGEBRAIC[xrss] is xrss in component IKr (dimensionless).
 * ALGEBRAIC[txrf] is txrf in component IKr (millisecond).
 * ALGEBRAIC[txrs] is txrs in component IKr (millisecond).
 * ALGEBRAIC[Axrf] is Axrf in component IKr (dimensionless).
 * ALGEBRAIC[Axrs] is Axrs in component IKr (dimensionless).
 * STATES[xrf] is xrf in component IKr (dimensionless).
 * STATES[xrs] is xrs in component IKr (dimensionless).
 * ALGEBRAIC[xr] is xr in component IKr (dimensionless).
 * ALGEBRAIC[rkr] is rkr in component IKr (dimensionless).
 * CONSTANTS[GKs_b] is GKs_b in component IKs (milliS_per_microF).
 * CONSTANTS[GKs] is GKs in component IKs (milliS_per_microF).
 * ALGEBRAIC[xs1ss] is xs1ss in component IKs (dimensionless).
 * ALGEBRAIC[xs2ss] is xs2ss in component IKs (dimensionless).
 * ALGEBRAIC[txs1] is txs1 in component IKs (millisecond).
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
 * CONSTANTS[PCab] is PCab in component ICab (milliS_per_microF).
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
 * ALGEBRAIC[tau_rel_temp] is tau_rel_temp in component ryr (millisecond).
 * ALGEBRAIC[tau_relp_temp] is tau_relp_temp in component ryr (millisecond).
 * CONSTANTS[upScale] is upScale in component SERCA (dimensionless).
 * ALGEBRAIC[Jupnp] is Jupnp in component SERCA (millimolar_per_millisecond).
 * ALGEBRAIC[Jupp] is Jupp in component SERCA (millimolar_per_millisecond).
 * ALGEBRAIC[fJupp] is fJupp in component SERCA (dimensionless).
 * ALGEBRAIC[Jleak] is Jleak in component SERCA (millimolar_per_millisecond).
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
 * RATES[xrf] is d/dt xrf in component IKr (dimensionless).
 * RATES[xrs] is d/dt xrs in component IKr (dimensionless).
 * RATES[xs1] is d/dt xs1 in component IKs (dimensionless).
 * RATES[xs2] is d/dt xs2 in component IKs (dimensionless).
 * RATES[xk1] is d/dt xk1 in component IK1 (dimensionless).
 * RATES[Jrelnp] is d/dt Jrelnp in component ryr (dimensionless).
 * RATES[Jrelp] is d/dt Jrelp in component ryr (dimensionless).
 */


Ohara_Rudy_2011::Ohara_Rudy_2011()
{
algebraic_size = 198;
constants_size = 139+3;
states_size = 41+3;
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
printf("DEALLOCATE P_CELL SUCCESS!!\n");
}

void Ohara_Rudy_2011::___initConsts()
{
STATES[V] = -87;
STATES[qnet] = 0.0;
STATES[INaL_AUC] = 0.0;
STATES[ICaL_AUC] = 0.0;
CONSTANTS[amp] = -80;
CONSTANTS[duration] = 0.5;
CONSTANTS[stim_start] = 10.0;
CONSTANTS[stim_period] = 1000.0;
CONSTANTS[R] = 8314;
CONSTANTS[T] = 310;
CONSTANTS[F] = 96485;
CONSTANTS[cm] = 1;
CONSTANTS[rad] = 0.0011;
CONSTANTS[L] = 0.01;
CONSTANTS[vcell] =  1000.00*3.14000*CONSTANTS[rad]*CONSTANTS[rad]*CONSTANTS[L];
#if defined CAB_ORUDY2011 
CONSTANTS[PCab] = 2.5e-8;
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
CONSTANTS[PCa_b] = 0.0001;
CONSTANTS[PCa] = (CONSTANTS[celltype]==1.00000 ?  CONSTANTS[PCa_b]*1.20000 : CONSTANTS[celltype]==2.00000 ?  CONSTANTS[PCa_b]*2.50000 : CONSTANTS[PCa_b]);
CONSTANTS[PCap] =  1.10000*CONSTANTS[PCa];
CONSTANTS[PCaK] =  0.000357400*CONSTANTS[PCa];
CONSTANTS[PCaKp] =  0.000357400*CONSTANTS[PCap];
CONSTANTS[ko] = 5.4;
#endif
#if defined CAL_ORUDY2011 
CONSTANTS[bt] = 4.75;
CONSTANTS[btp] =  1.25000*CONSTANTS[bt];
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
CONSTANTS[k2n] = 1000;
CONSTANTS[tjca] = 75.0000;
CONSTANTS[Aff] = 0.600000;
CONSTANTS[Afs] = 1.00000 - CONSTANTS[Aff];
CONSTANTS[cao] = 1.8;
CONSTANTS[KmCaMK] = 0.15;
CONSTANTS[PCa_b] = 0.0001;
CONSTANTS[PCa] = (CONSTANTS[celltype]==1.00000 ?  CONSTANTS[PCa_b]*1.20000 : CONSTANTS[celltype]==2.00000 ?  CONSTANTS[PCa_b]*2.50000 : CONSTANTS[PCa_b]);
CONSTANTS[PCap] =  1.10000*CONSTANTS[PCa];
#endif
#if defined CANA_ORUDY2011 
CONSTANTS[bt] = 4.75;
CONSTANTS[btp] =  1.25000*CONSTANTS[bt];
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
CONSTANTS[k2n] = 1000;
CONSTANTS[tjca] = 75.0000;
CONSTANTS[Aff] = 0.600000;
CONSTANTS[Afs] = 1.00000 - CONSTANTS[Aff];
CONSTANTS[cao] = 1.8;
CONSTANTS[KmCaMK] = 0.15;
CONSTANTS[PCa_b] = 0.0001;
CONSTANTS[PCa] = (CONSTANTS[celltype]==1.00000 ?  CONSTANTS[PCa_b]*1.20000 : CONSTANTS[celltype]==2.00000 ?  CONSTANTS[PCa_b]*2.50000 : CONSTANTS[PCa_b]);
CONSTANTS[PCap] =  1.10000*CONSTANTS[PCa];
CONSTANTS[PCaNa] =  0.00125000*CONSTANTS[PCa];
CONSTANTS[PCaNap] =  0.00125000*CONSTANTS[PCap];
CONSTANTS[nao] = 140;
#endif
#if defined K1_ORUDY2011 
CONSTANTS[GK1_b] = 0.1908;
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
CONSTANTS[GKr_b] = 0.046;
CONSTANTS[GKr] = (CONSTANTS[celltype]==1.00000 ?  CONSTANTS[GKr_b]*1.30000 : CONSTANTS[celltype]==2.00000 ?  CONSTANTS[GKr_b]*0.800000 : CONSTANTS[GKr_b]);
CONSTANTS[ko] = 5.4;
CONSTANTS[Ageo] =  2.00000*3.14000*CONSTANTS[rad]*CONSTANTS[rad]+ 2.00000*3.14000*CONSTANTS[rad]*CONSTANTS[L];
CONSTANTS[Acap] =  2.00000*CONSTANTS[Ageo];
CONSTANTS[vmyo] =  0.680000*CONSTANTS[vcell];
CONSTANTS[vss] =  0.0200000*CONSTANTS[vcell];
#endif
#if defined KS_ORUDY2011 
CONSTANTS[GKs_b] = 0.0034;
CONSTANTS[GKs] = (CONSTANTS[celltype]==1.00000 ?  CONSTANTS[GKs_b]*1.40000 : CONSTANTS[GKs_b]);
CONSTANTS[ko] = 5.4;
CONSTANTS[Ageo] =  2.00000*3.14000*CONSTANTS[rad]*CONSTANTS[rad]+ 2.00000*3.14000*CONSTANTS[rad]*CONSTANTS[L];
CONSTANTS[Acap] =  2.00000*CONSTANTS[Ageo];
CONSTANTS[vmyo] =  0.680000*CONSTANTS[vcell];
CONSTANTS[vss] =  0.0200000*CONSTANTS[vcell];
#endif
#if defined NA_ORUDY2011 
CONSTANTS[GNa] = 75;
CONSTANTS[CaMKo] = 0.05;
CONSTANTS[KmCaM] = 0.0015;
CONSTANTS[KmCaMK] = 0.15;
CONSTANTS[nao] = 140;
CONSTANTS[Ahf] = 0.99;
CONSTANTS[Ahs] = 1.00000 - CONSTANTS[Ahf];
CONSTANTS[Ageo] =  2.00000*3.14000*CONSTANTS[rad]*CONSTANTS[rad]+ 2.00000*3.14000*CONSTANTS[rad]*CONSTANTS[L];
CONSTANTS[Acap] =  2.00000*CONSTANTS[Ageo];
CONSTANTS[vmyo] =  0.680000*CONSTANTS[vcell];
CONSTANTS[vss] =  0.0200000*CONSTANTS[vcell];
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
#endif
#if defined NAB_ORUDY2011 
CONSTANTS[PNab] = 3.75e-10;
CONSTANTS[nao] = 140;
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
CONSTANTS[GNaL_b] = 0.0075;
CONSTANTS[GNaL] = (CONSTANTS[celltype]==1.00000 ?  CONSTANTS[GNaL_b]*0.600000 : CONSTANTS[GNaL_b]);
CONSTANTS[CaMKo] = 0.05;
CONSTANTS[KmCaM] = 0.0015;
CONSTANTS[KmCaMK] = 0.15;
CONSTANTS[nao] = 140;
CONSTANTS[thL] = 200;
CONSTANTS[thLp] =  3.00000*CONSTANTS[thL];
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
#endif
#if defined PCA_ORUDY2011 
CONSTANTS[GpCa] = 0.0005;
CONSTANTS[KmCap] = 0.0005;
CONSTANTS[vss] =  0.0200000*CONSTANTS[vcell];
CONSTANTS[Ageo] =  2.00000*3.14000*CONSTANTS[rad]*CONSTANTS[rad]+ 2.00000*3.14000*CONSTANTS[rad]*CONSTANTS[L];
CONSTANTS[Acap] =  2.00000*CONSTANTS[Ageo];
CONSTANTS[vmyo] =  0.680000*CONSTANTS[vcell];
CONSTANTS[vnsr] =  0.0552000*CONSTANTS[vcell];
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
STATES[cai] = 1e-4;
STATES[nass] = 7;
STATES[nai] = 7;
STATES[cansr] = 1.2;
STATES[CaMKt] = 0;
STATES[cajsr] = 1.2;
STATES[Jrelnp] = 0;
STATES[Jrelp] = 0;
STATES[cass] = 1e-4;
#endif
#if defined CAK_ORUDY2011 
STATES[cass] = 1e-4;
STATES[d] = 0;
STATES[CaMKt] = 0;
STATES[ff] = 1;
STATES[fs] = 1;
STATES[fcaf] = 1;
STATES[nca] = 0;
STATES[jca] = 1;
STATES[fcas] = 1;
STATES[ffp] = 1;
STATES[fcafp] = 1;
STATES[Jrelnp] = 0;
STATES[Jrelp] = 0;
STATES[cai] = 1e-4;
STATES[cansr] = 1.2;
STATES[cajsr] = 1.2;
STATES[ki] = 145;
STATES[kss] = 145;
#endif
#if defined CAL_ORUDY2011 
STATES[cass] = 1e-4;
STATES[d] = 0;
STATES[CaMKt] = 0;
STATES[ff] = 1;
STATES[fs] = 1;
STATES[fcaf] = 1;
STATES[nca] = 0;
STATES[jca] = 1;
STATES[fcas] = 1;
STATES[ffp] = 1;
STATES[fcafp] = 1;
STATES[Jrelnp] = 0;
STATES[Jrelp] = 0;
STATES[cai] = 1e-4;
STATES[cansr] = 1.2;
STATES[cajsr] = 1.2;
#endif
#if defined CANA_ORUDY2011 
STATES[d] = 0;
STATES[cass] = 1e-4;
STATES[ff] = 1;
STATES[fs] = 1;
STATES[nca] = 0;
STATES[jca] = 1;
STATES[fcaf] = 1;
STATES[ffp] = 1;
STATES[fcafp] = 1;
STATES[fcas] = 1;
STATES[CaMKt] = 0;
STATES[nass] = 7;
STATES[nai] = 7;
STATES[cansr] = 1.2;
STATES[cajsr] = 1.2;
STATES[cai] = 1e-4;
STATES[Jrelnp] = 0;
STATES[Jrelp] = 0;
#endif
#if defined K1_ORUDY2011 
STATES[xk1] = 1;
STATES[ki] = 145;
STATES[kss] = 145;
#endif
#if defined KB_ORUDY2011 
STATES[ki] = 145;
STATES[kss] = 145;
#endif
#if defined KR_ORUDY2011 
STATES[ki] = 145;
STATES[xrf] = 0;
STATES[xrs] = 0;
STATES[kss] = 145;
#endif
#if defined KS_ORUDY2011 
STATES[xs1] = 0;
STATES[xs2] = 0;
STATES[ki] = 145;
STATES[kss] = 145;
#endif
#if defined NA_ORUDY2011 
STATES[m] = 0;
STATES[j] = 1;
STATES[jp] = 1;
STATES[nai] = 7;
STATES[hf] = 1;
STATES[hs] = 1;
STATES[hsp] = 1;
STATES[CaMKt] = 0;
#endif
#if defined NAB_ORUDY2011 
STATES[nass] = 7;
STATES[nai] = 7;
#endif
#if defined NACA_I_ORUDY2011 
STATES[cai] = 1e-4;
STATES[nai] = 7;
STATES[cansr] = 1.2;
STATES[cass] = 1e-4;
STATES[nass] = 7;
STATES[CaMKt] = 0;
STATES[cajsr] = 1.2;
STATES[Jrelnp] = 0;
STATES[Jrelp] = 0;
#endif
#if defined NACA_SS_ORUDY2011 
STATES[cai] = 1e-4;
STATES[nai] = 7;
STATES[cansr] = 1.2;
STATES[cass] = 1e-4;
STATES[nass] = 7;
STATES[CaMKt] = 0;
STATES[cajsr] = 1.2;
STATES[Jrelnp] = 0;
STATES[Jrelp] = 0;
#endif
#if defined NAK_ORUDY2011 
STATES[nai] = 7;
STATES[ki] = 145;
STATES[kss] = 145;
STATES[nass] = 7;
#endif
#if defined NAL_ORUDY2011 
STATES[mL] = 0;
STATES[hL] = 1;
STATES[hLp] = 1;
STATES[CaMKt] = 0;
STATES[nai] = 7;
STATES[nass] = 7;
#endif
#if defined PCA_ORUDY2011 
STATES[cai] = 1e-4;
STATES[cansr] = 1.2;
STATES[CaMKt] = 0;
STATES[cajsr] = 1.2;
STATES[Jrelnp] = 0;
STATES[Jrelp] = 0;
STATES[cass] = 1e-4;
#endif
#if defined TO_ORUDY2011 
STATES[a] = 0;
STATES[ki] = 145;
STATES[iF] = 1;
STATES[iS] = 1;
STATES[kss] = 145;
STATES[CaMKt] = 0;
#endif
}

void Ohara_Rudy_2011::initConsts()
{
CONSTANTS[celltype] = 0;
___initConsts();
}

void Ohara_Rudy_2011::initConsts(double type)
{
CONSTANTS[celltype] = type;
___initConsts();
}

void Ohara_Rudy_2011::initConsts(double type, double D, double *hill, bool is_dutta)
{
initConsts(type);
// Apply Dutta's modification from the original ORd model
if(is_dutta == true){
  // new dutta
  CONSTANTS[GKs] *= 1.870;
  CONSTANTS[GKr] *= 1.013;
  CONSTANTS[GK1] *= 1.698;
  CONSTANTS[PCa] *= 1.007;
  CONSTANTS[GNaL] *= 2.661;
}
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
#if defined SINGLE_CELL
//ALGEBRAIC[Istim] = (TIME<=CONSTANTS[duration] ? CONSTANTS[amp] : 0.000000);
ALGEBRAIC[Istim] = (TIME>=CONSTANTS[stim_start] && (TIME - CONSTANTS[stim_start]) - floor((TIME - CONSTANTS[stim_start])/CONSTANTS[stim_period])*CONSTANTS[stim_period]<=CONSTANTS[duration] ? CONSTANTS[amp] : 0.000000);
#endif
#if defined CAB_ORUDY2011 
ALGEBRAIC[vffrt] = ( STATES[V]*CONSTANTS[F]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[T]);
ALGEBRAIC[vfrt] = ( STATES[V]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[T]);
ALGEBRAIC[ICab] = ( CONSTANTS[PCab]*4.00000*ALGEBRAIC[vffrt]*( STATES[cai]*exp( 2.00000*ALGEBRAIC[vfrt]) -  0.341000*CONSTANTS[cao]))/(exp( 2.00000*ALGEBRAIC[vfrt]) - 1.00000);
ALGEBRAIC[JdiffNa] = (STATES[nass] - STATES[nai])/2.00000;
ALGEBRAIC[Bcai] = 1.00000/(1.00000+( CONSTANTS[cmdnmax]*CONSTANTS[kmcmdn])/pow(CONSTANTS[kmcmdn]+STATES[cai], 2.00000)+( CONSTANTS[trpnmax]*CONSTANTS[kmtrpn])/pow(CONSTANTS[kmtrpn]+STATES[cai], 2.00000));
ALGEBRAIC[CaMKb] = ( CONSTANTS[CaMKo]*(1.00000 - STATES[CaMKt]))/(1.00000+CONSTANTS[KmCaM]/STATES[cass]);
ALGEBRAIC[CaMKa] = ALGEBRAIC[CaMKb]+STATES[CaMKt];
ALGEBRAIC[fJupp] = 1.00000/(1.00000+CONSTANTS[KmCaMK]/ALGEBRAIC[CaMKa]);
ALGEBRAIC[Jupnp] = ( CONSTANTS[upScale]*0.00437500*STATES[cai])/(STATES[cai]+0.000920000);
ALGEBRAIC[Jupp] = ( CONSTANTS[upScale]*2.75000*0.00437500*STATES[cai])/((STATES[cai]+0.000920000) - 0.000170000);
ALGEBRAIC[Jleak] = ( 0.00393750*STATES[cansr])/15.0000;
ALGEBRAIC[Jup] = ( (1.00000 - ALGEBRAIC[fJupp])*ALGEBRAIC[Jupnp]+ ALGEBRAIC[fJupp]*ALGEBRAIC[Jupp]) - ALGEBRAIC[Jleak];
ALGEBRAIC[Jtr] = (STATES[cansr] - STATES[cajsr])/100.000;
ALGEBRAIC[Bcajsr] = 1.00000/(1.00000+( CONSTANTS[csqnmax]*CONSTANTS[kmcsqn])/pow(CONSTANTS[kmcsqn]+STATES[cajsr], 2.00000));
ALGEBRAIC[Bcass] = 1.00000/(1.00000+( CONSTANTS[BSRmax]*CONSTANTS[KmBSR])/pow(CONSTANTS[KmBSR]+STATES[cass], 2.00000)+( CONSTANTS[BSLmax]*CONSTANTS[KmBSL])/pow(CONSTANTS[KmBSL]+STATES[cass], 2.00000));
ALGEBRAIC[Jdiff] = (STATES[cass] - STATES[cai])/0.200000;
ALGEBRAIC[JdiffNa] = (STATES[nass] - STATES[nai])/2.00000;
ALGEBRAIC[fJrelp] = 1.00000/(1.00000+CONSTANTS[KmCaMK]/ALGEBRAIC[CaMKa]);
ALGEBRAIC[Jrel] =  (1.00000 - ALGEBRAIC[fJrelp])*STATES[Jrelnp]+ ALGEBRAIC[fJrelp]*STATES[Jrelp];
ALGEBRAIC[Jrel_inf_temp] = ( CONSTANTS[a_rel]*- ALGEBRAIC[ICaL])/(1.00000+ 1.00000*pow(1.50000/STATES[cajsr], 8.00000));
ALGEBRAIC[tau_relp_temp] = CONSTANTS[btp]/(1.00000+0.0123000/STATES[cajsr]);
ALGEBRAIC[Jrel_temp] = ( CONSTANTS[a_relp]*- ALGEBRAIC[ICaL])/(1.00000+pow(1.50000/STATES[cajsr], 8.00000));
ALGEBRAIC[tau_rel_temp] = CONSTANTS[bt]/(1.00000+0.0123000/STATES[cajsr]);
ALGEBRAIC[tau_rel] = (ALGEBRAIC[tau_rel_temp]<0.00100000 ? 0.00100000 : ALGEBRAIC[tau_rel_temp]);
ALGEBRAIC[Jrel_inf] = (CONSTANTS[celltype]==2.00000 ?  ALGEBRAIC[Jrel_inf_temp]*1.70000 : ALGEBRAIC[Jrel_inf_temp]);
ALGEBRAIC[Jrel_infp] = (CONSTANTS[celltype]==2.00000 ?  ALGEBRAIC[Jrel_temp]*1.70000 : ALGEBRAIC[Jrel_temp]);
ALGEBRAIC[tau_relp] = (ALGEBRAIC[tau_relp_temp]<0.00100000 ? 0.00100000 : ALGEBRAIC[tau_relp_temp]);
#endif
#if defined CAK_ORUDY2011 
ALGEBRAIC[CaMKb] = ( CONSTANTS[CaMKo]*(1.00000 - STATES[CaMKt]))/(1.00000+CONSTANTS[KmCaM]/STATES[cass]);
ALGEBRAIC[CaMKa] = ALGEBRAIC[CaMKb]+STATES[CaMKt];
ALGEBRAIC[vffrt] = ( STATES[V]*CONSTANTS[F]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[T]);
ALGEBRAIC[vfrt] = ( STATES[V]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[T]);
ALGEBRAIC[fICaLp] = 1.00000/(1.00000+CONSTANTS[KmCaMK]/ALGEBRAIC[CaMKa]);
ALGEBRAIC[PhiCaK] = ( 1.00000*ALGEBRAIC[vffrt]*( 0.750000*STATES[kss]*exp( 1.00000*ALGEBRAIC[vfrt]) -  0.750000*CONSTANTS[ko]))/(exp( 1.00000*ALGEBRAIC[vfrt]) - 1.00000);
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
ALGEBRAIC[ths] = 1.00000/( 0.00979400*exp(- (STATES[V]+17.9500)/28.0500)+ 0.334300*exp((STATES[V]+5.73000)/56.6600));
ALGEBRAIC[tffp] =  2.50000*ALGEBRAIC[tff];
ALGEBRAIC[fcass] = ALGEBRAIC[fss];
ALGEBRAIC[tfcas] = 100.000+1.00000/( 0.000120000*exp(- STATES[V]/3.00000)+ 0.000120000*exp(STATES[V]/7.00000));
ALGEBRAIC[Bcass] = 1.00000/(1.00000+( CONSTANTS[BSRmax]*CONSTANTS[KmBSR])/pow(CONSTANTS[KmBSR]+STATES[cass], 2.00000)+( CONSTANTS[BSLmax]*CONSTANTS[KmBSL])/pow(CONSTANTS[KmBSL]+STATES[cass], 2.00000));
ALGEBRAIC[fJrelp] = 1.00000/(1.00000+CONSTANTS[KmCaMK]/ALGEBRAIC[CaMKa]);
ALGEBRAIC[Jrel] =  (1.00000 - ALGEBRAIC[fJrelp])*STATES[Jrelnp]+ ALGEBRAIC[fJrelp]*STATES[Jrelp];
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
ALGEBRAIC[Jup] = ( (1.00000 - ALGEBRAIC[fJupp])*ALGEBRAIC[Jupnp]+ ALGEBRAIC[fJupp]*ALGEBRAIC[Jupp]) - ALGEBRAIC[Jleak];
ALGEBRAIC[Bcajsr] = 1.00000/(1.00000+( CONSTANTS[csqnmax]*CONSTANTS[kmcsqn])/pow(CONSTANTS[kmcsqn]+STATES[cajsr], 2.00000));
ALGEBRAIC[Jtr] = (STATES[cansr] - STATES[cajsr])/100.000;
#endif
#if defined CAL_ORUDY2011 
ALGEBRAIC[CaMKb] = ( CONSTANTS[CaMKo]*(1.00000 - STATES[CaMKt]))/(1.00000+CONSTANTS[KmCaM]/STATES[cass]);
ALGEBRAIC[CaMKa] = ALGEBRAIC[CaMKb]+STATES[CaMKt];
ALGEBRAIC[fICaLp] = 1.00000/(1.00000+CONSTANTS[KmCaMK]/ALGEBRAIC[CaMKa]);
ALGEBRAIC[vffrt] = ( STATES[V]*CONSTANTS[F]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[T]);
ALGEBRAIC[vfrt] = ( STATES[V]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[T]);
ALGEBRAIC[PhiCaL] = ( 4.00000*ALGEBRAIC[vffrt]*( STATES[cass]*exp( 2.00000*ALGEBRAIC[vfrt]) -  0.341000*CONSTANTS[cao]))/(exp( 2.00000*ALGEBRAIC[vfrt]) - 1.00000);
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
ALGEBRAIC[Jup] = ( (1.00000 - ALGEBRAIC[fJupp])*ALGEBRAIC[Jupnp]+ ALGEBRAIC[fJupp]*ALGEBRAIC[Jupp]) - ALGEBRAIC[Jleak];
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
ALGEBRAIC[Jrel] =  (1.00000 - ALGEBRAIC[fJrelp])*STATES[Jrelnp]+ ALGEBRAIC[fJrelp]*STATES[Jrelp];
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
ALGEBRAIC[vffrt] = ( STATES[V]*CONSTANTS[F]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[T]);
ALGEBRAIC[vfrt] = ( STATES[V]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[T]);
ALGEBRAIC[PhiCaNa] = ( 1.00000*ALGEBRAIC[vffrt]*( 0.750000*STATES[nass]*exp( 1.00000*ALGEBRAIC[vfrt]) -  0.750000*CONSTANTS[nao]))/(exp( 1.00000*ALGEBRAIC[vfrt]) - 1.00000);
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
ALGEBRAIC[Jup] = ( (1.00000 - ALGEBRAIC[fJupp])*ALGEBRAIC[Jupnp]+ ALGEBRAIC[fJupp]*ALGEBRAIC[Jupp]) - ALGEBRAIC[Jleak];
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
ALGEBRAIC[Jrel] =  (1.00000 - ALGEBRAIC[fJrelp])*STATES[Jrelnp]+ ALGEBRAIC[fJrelp]*STATES[Jrelp];
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
ALGEBRAIC[rkr] = ( (1.00000/(1.00000+exp((STATES[V]+55.0000)/75.0000)))*1.00000)/(1.00000+exp((STATES[V] - 10.0000)/30.0000));
ALGEBRAIC[xrss] = 1.00000/(1.00000+exp(- (STATES[V]+8.33700)/6.78900));
ALGEBRAIC[txrf] = 12.9800+1.00000/( 0.365200*exp((STATES[V] - 31.6600)/3.86900)+ 4.12300e-05*exp(- (STATES[V] - 47.7800)/20.3800));
ALGEBRAIC[txrs] = 1.86500+1.00000/( 0.0662900*exp((STATES[V] - 34.7000)/7.35500)+ 1.12800e-05*exp(- (STATES[V] - 29.7400)/25.9400));
ALGEBRAIC[Axrf] = 1.00000/(1.00000+exp((STATES[V]+54.8100)/38.2100));
ALGEBRAIC[Axrs] = 1.00000 - ALGEBRAIC[Axrf];
ALGEBRAIC[xr] =  ALGEBRAIC[Axrf]*STATES[xrf]+ ALGEBRAIC[Axrs]*STATES[xrs];
ALGEBRAIC[IKr] =  CONSTANTS[GKr]* pow((CONSTANTS[ko]/5.40000), 1.0 / 2)*ALGEBRAIC[xr]*ALGEBRAIC[rkr]*(STATES[V] - ALGEBRAIC[EK]);
ALGEBRAIC[JdiffK] = (STATES[kss] - STATES[ki])/2.00000;
#endif
#if defined KS_ORUDY2011 
ALGEBRAIC[EKs] =  (( CONSTANTS[R]*CONSTANTS[T])/CONSTANTS[F])*log((CONSTANTS[ko]+ CONSTANTS[PKNa]*CONSTANTS[nao])/(STATES[ki]+ CONSTANTS[PKNa]*STATES[nai]));
ALGEBRAIC[KsCa] = 1.00000+0.600000/(1.00000+pow(3.80000e-05/STATES[cai], 1.40000));
ALGEBRAIC[IKs] =  CONSTANTS[GKs]*ALGEBRAIC[KsCa]*STATES[xs1]*STATES[xs2]*(STATES[V] - ALGEBRAIC[EKs]);
ALGEBRAIC[xs1ss] = 1.00000/(1.00000+exp(- (STATES[V]+11.6000)/8.93200));
ALGEBRAIC[xs2ss] = ALGEBRAIC[xs1ss];
ALGEBRAIC[txs1] = 817.300+1.00000/( 0.000232600*exp((STATES[V]+48.2800)/17.8000)+ 0.00129200*exp(- (STATES[V]+210.000)/230.000));
ALGEBRAIC[txs2] = 1.00000/( 0.0100000*exp((STATES[V] - 50.0000)/20.0000)+ 0.0193000*exp(- (STATES[V]+66.5400)/31.0000));
ALGEBRAIC[JdiffK] = (STATES[kss] - STATES[ki])/2.00000;
#endif
#if defined NA_ORUDY2011 
ALGEBRAIC[ENa] =  (( CONSTANTS[R]*CONSTANTS[T])/CONSTANTS[F])*log(CONSTANTS[nao]/STATES[nai]);
ALGEBRAIC[CaMKb] = ( CONSTANTS[CaMKo]*(1.00000 - STATES[CaMKt]))/(1.00000+CONSTANTS[KmCaM]/STATES[cass]);
ALGEBRAIC[CaMKa] = ALGEBRAIC[CaMKb]+STATES[CaMKt];
ALGEBRAIC[fINap] = 1.00000/(1.00000+CONSTANTS[KmCaMK]/ALGEBRAIC[CaMKa]);
ALGEBRAIC[h] =  CONSTANTS[Ahf]*STATES[hf]+ CONSTANTS[Ahs]*STATES[hs];
ALGEBRAIC[hp] =  CONSTANTS[Ahf]*STATES[hf]+ CONSTANTS[Ahs]*STATES[hsp];
ALGEBRAIC[INa] =  CONSTANTS[GNa]*(STATES[V] - ALGEBRAIC[ENa])*pow(STATES[m], 3.00000)*( (1.00000 - ALGEBRAIC[fINap])*ALGEBRAIC[h]*STATES[j]+ ALGEBRAIC[fINap]*ALGEBRAIC[hp]*STATES[jp]);
ALGEBRAIC[mss] = 1.00000/(1.00000+exp(- (STATES[V]+CONSTANTS[mssV1])/CONSTANTS[mssV2]));
ALGEBRAIC[tm] = 1.00000/( CONSTANTS[mtD1]*exp((STATES[V]+CONSTANTS[mtV1])/CONSTANTS[mtV2])+ CONSTANTS[mtD2]*exp(- (STATES[V]+CONSTANTS[mtV3])/CONSTANTS[mtV4]));
ALGEBRAIC[hss] = 1.00000/(1.00000+exp((STATES[V]+CONSTANTS[hssV1])/CONSTANTS[hssV2]));
ALGEBRAIC[ths] = 1.00000/( 0.00979400*exp(- (STATES[V]+17.9500)/28.0500)+ 0.334300*exp((STATES[V]+5.73000)/56.6600));
ALGEBRAIC[thf] = 1.00000/( 1.43200e-05*exp(- (STATES[V]+1.19600)/6.28500)+ 6.14900*exp((STATES[V]+0.509600)/20.2700));
ALGEBRAIC[tj] = 2.03800+1.00000/( 0.0213600*exp(- (STATES[V]+100.600)/8.28100)+ 0.305200*exp((STATES[V]+0.994100)/38.4500));
ALGEBRAIC[jss] = ALGEBRAIC[hss];
ALGEBRAIC[JdiffNa] = (STATES[nass] - STATES[nai])/2.00000;
ALGEBRAIC[tjp] =  1.46000*ALGEBRAIC[tj];
ALGEBRAIC[hssp] = 1.00000/(1.00000+exp((STATES[V]+89.1000)/6.08600));
ALGEBRAIC[thsp] =  3.00000*ALGEBRAIC[ths];
#endif
#if defined NAB_ORUDY2011 
ALGEBRAIC[vffrt] = ( STATES[V]*CONSTANTS[F]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[T]);
ALGEBRAIC[vfrt] = ( STATES[V]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[T]);
ALGEBRAIC[INab] = ( CONSTANTS[PNab]*ALGEBRAIC[vffrt]*( STATES[nai]*exp(ALGEBRAIC[vfrt]) - CONSTANTS[nao]))/(exp(ALGEBRAIC[vfrt]) - 1.00000);
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
ALGEBRAIC[Jup] = ( (1.00000 - ALGEBRAIC[fJupp])*ALGEBRAIC[Jupnp]+ ALGEBRAIC[fJupp]*ALGEBRAIC[Jupp]) - ALGEBRAIC[Jleak];
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
ALGEBRAIC[Jrel] =  (1.00000 - ALGEBRAIC[fJrelp])*STATES[Jrelnp]+ ALGEBRAIC[fJrelp]*STATES[Jrelp];
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
ALGEBRAIC[Jup] = ( (1.00000 - ALGEBRAIC[fJupp])*ALGEBRAIC[Jupnp]+ ALGEBRAIC[fJupp]*ALGEBRAIC[Jupp]) - ALGEBRAIC[Jleak];
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
ALGEBRAIC[Jrel] =  (1.00000 - ALGEBRAIC[fJrelp])*STATES[Jrelnp]+ ALGEBRAIC[fJrelp]*STATES[Jrelp];
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
ALGEBRAIC[ENa] =  (( CONSTANTS[R]*CONSTANTS[T])/CONSTANTS[F])*log(CONSTANTS[nao]/STATES[nai]);
ALGEBRAIC[CaMKb] = ( CONSTANTS[CaMKo]*(1.00000 - STATES[CaMKt]))/(1.00000+CONSTANTS[KmCaM]/STATES[cass]);
ALGEBRAIC[CaMKa] = ALGEBRAIC[CaMKb]+STATES[CaMKt];
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
ALGEBRAIC[Jup] = ( (1.00000 - ALGEBRAIC[fJupp])*ALGEBRAIC[Jupnp]+ ALGEBRAIC[fJupp]*ALGEBRAIC[Jupp]) - ALGEBRAIC[Jleak];
ALGEBRAIC[Jtr] = (STATES[cansr] - STATES[cajsr])/100.000;
ALGEBRAIC[Bcajsr] = 1.00000/(1.00000+( CONSTANTS[csqnmax]*CONSTANTS[kmcsqn])/pow(CONSTANTS[kmcsqn]+STATES[cajsr], 2.00000));
ALGEBRAIC[fJrelp] = 1.00000/(1.00000+CONSTANTS[KmCaMK]/ALGEBRAIC[CaMKa]);
ALGEBRAIC[Jrel] =  (1.00000 - ALGEBRAIC[fJrelp])*STATES[Jrelnp]+ ALGEBRAIC[fJrelp]*STATES[Jrelp];
ALGEBRAIC[Jrel_inf_temp] = ( CONSTANTS[a_rel]*- ALGEBRAIC[ICaL])/(1.00000+ 1.00000*pow(1.50000/STATES[cajsr], 8.00000));
ALGEBRAIC[tau_rel_temp] = CONSTANTS[bt]/(1.00000+0.0123000/STATES[cajsr]);
ALGEBRAIC[Jrel_temp] = ( CONSTANTS[a_relp]*- ALGEBRAIC[ICaL])/(1.00000+pow(1.50000/STATES[cajsr], 8.00000));
ALGEBRAIC[tau_rel] = (ALGEBRAIC[tau_rel_temp]<0.00100000 ? 0.00100000 : ALGEBRAIC[tau_rel_temp]);
ALGEBRAIC[Jrel_inf] = (CONSTANTS[celltype]==2.00000 ?  ALGEBRAIC[Jrel_inf_temp]*1.70000 : ALGEBRAIC[Jrel_inf_temp]);
ALGEBRAIC[Jrel_infp] = (CONSTANTS[celltype]==2.00000 ?  ALGEBRAIC[Jrel_temp]*1.70000 : ALGEBRAIC[Jrel_temp]);
ALGEBRAIC[tau_relp] = (ALGEBRAIC[tau_relp_temp]<0.00100000 ? 0.00100000 : ALGEBRAIC[tau_relp_temp]);
#endif
#if defined TO_ORUDY2011 
ALGEBRAIC[EK] =  (( CONSTANTS[R]*CONSTANTS[T])/CONSTANTS[F])*log(CONSTANTS[ko]/STATES[ki]);
ALGEBRAIC[CaMKb] = ( CONSTANTS[CaMKo]*(1.00000 - STATES[CaMKt]))/(1.00000+CONSTANTS[KmCaM]/STATES[cass]);
ALGEBRAIC[CaMKa] = ALGEBRAIC[CaMKb]+STATES[CaMKt];
ALGEBRAIC[fItop] = 1.00000/(1.00000+CONSTANTS[KmCaMK]/ALGEBRAIC[CaMKa]);
ALGEBRAIC[AiF] = 1.00000/(1.00000+exp((STATES[V] - 213.600)/151.200));
ALGEBRAIC[AiS] = 1.00000 - ALGEBRAIC[AiF];
ALGEBRAIC[i] =  ALGEBRAIC[AiF]*STATES[iF]+ ALGEBRAIC[AiS]*STATES[iS];
ALGEBRAIC[Ito] =  CONSTANTS[Gto]*(STATES[V] - ALGEBRAIC[EK])*( (1.00000 - ALGEBRAIC[fItop])*STATES[a]*ALGEBRAIC[i]+ ALGEBRAIC[fItop]*STATES[ap]*ALGEBRAIC[ip]);
ALGEBRAIC[ass] = 1.00000/(1.00000+exp(- (STATES[V] - 14.3400)/14.8200));
ALGEBRAIC[ta] = 1.05150/(1.00000/( 1.20890*(1.00000+exp(- (STATES[V] - 18.4099)/29.3814)))+3.50000/(1.00000+exp((STATES[V]+100.000)/29.3814)));
ALGEBRAIC[JdiffK] = (STATES[kss] - STATES[ki])/2.00000;
ALGEBRAIC[iss] = 1.00000/(1.00000+exp((STATES[V]+43.9400)/5.71100));
ALGEBRAIC[tiF_b] = 4.56200+1.00000/( 0.393300*exp(- (STATES[V]+100.000)/100.000)+ 0.0800400*exp((STATES[V]+50.0000)/16.5900));
ALGEBRAIC[delta_epi] = (CONSTANTS[celltype]==1.00000 ? 1.00000 - 0.950000/(1.00000+exp((STATES[V]+70.0000)/5.00000)) : 1.00000);
ALGEBRAIC[tiS_b] = 23.6200+1.00000/( 0.00141600*exp(- (STATES[V]+96.5200)/59.0500)+ 1.78000e-08*exp((STATES[V]+114.100)/8.07900));
ALGEBRAIC[tiF] =  ALGEBRAIC[tiF_b]*ALGEBRAIC[delta_epi];
ALGEBRAIC[tiS] =  ALGEBRAIC[tiS_b]*ALGEBRAIC[delta_epi];
#endif


#if defined SINGLE_CELL
RATES[V] = - (ALGEBRAIC[INa]+ALGEBRAIC[INaL]+ALGEBRAIC[Ito]+ALGEBRAIC[ICaL]+ALGEBRAIC[ICaNa]+ALGEBRAIC[ICaK]+ALGEBRAIC[IKr]+ALGEBRAIC[IKs]+ALGEBRAIC[IK1]+ALGEBRAIC[INaCa_i]+ALGEBRAIC[INaCa_ss]+ALGEBRAIC[INaK]+ALGEBRAIC[INab]+ALGEBRAIC[IKb]+ALGEBRAIC[IpCa]+ALGEBRAIC[ICab]+ALGEBRAIC[Istim]);
RATES[qnet] = (ALGEBRAIC[INaL]+ALGEBRAIC[ICaL]+ALGEBRAIC[Ito]+ALGEBRAIC[IKr]+ALGEBRAIC[IKs]+ALGEBRAIC[IK1]);
RATES[INaL_AUC] = ALGEBRAIC[INaL];
RATES[ICaL_AUC] = ALGEBRAIC[ICaL];
#endif
#if defined CAB_ORUDY2011 
RATES[cai] =  ALGEBRAIC[Bcai]*((( - ((ALGEBRAIC[IpCa]+ALGEBRAIC[ICab]) -  2.00000*ALGEBRAIC[INaCa_i])*CONSTANTS[cm]*CONSTANTS[Acap])/( 2.00000*CONSTANTS[F]*CONSTANTS[vmyo]) - ( ALGEBRAIC[Jup]*CONSTANTS[vnsr])/CONSTANTS[vmyo])+( ALGEBRAIC[Jdiff]*CONSTANTS[vss])/CONSTANTS[vmyo]);
RATES[nass] = ( - (ALGEBRAIC[ICaNa]+ 3.00000*ALGEBRAIC[INaCa_ss])*CONSTANTS[cm]*CONSTANTS[Acap])/( CONSTANTS[F]*CONSTANTS[vss]) - ALGEBRAIC[JdiffNa];
RATES[nai] = ( - (ALGEBRAIC[INa]+ALGEBRAIC[INaL]+ 3.00000*ALGEBRAIC[INaCa_i]+ 3.00000*ALGEBRAIC[INaK]+ALGEBRAIC[INab])*CONSTANTS[Acap]*CONSTANTS[cm])/( CONSTANTS[F]*CONSTANTS[vmyo])+( ALGEBRAIC[JdiffNa]*CONSTANTS[vss])/CONSTANTS[vmyo];
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
RATES[xrf] = (ALGEBRAIC[xrss] - STATES[xrf])/ALGEBRAIC[txrf];
RATES[xrs] = (ALGEBRAIC[xrss] - STATES[xrs])/ALGEBRAIC[txrs];
RATES[kss] = ( - ALGEBRAIC[ICaK]*CONSTANTS[cm]*CONSTANTS[Acap])/( CONSTANTS[F]*CONSTANTS[vss]) - ALGEBRAIC[JdiffK];
#endif
#if defined KS_ORUDY2011 
RATES[xs1] = (ALGEBRAIC[xs1ss] - STATES[xs1])/ALGEBRAIC[txs1];
RATES[xs2] = (ALGEBRAIC[xs2ss] - STATES[xs2])/ALGEBRAIC[txs2];
RATES[ki] = ( - ((ALGEBRAIC[Ito]+ALGEBRAIC[IKr]+ALGEBRAIC[IKs]+ALGEBRAIC[IK1]+ALGEBRAIC[IKb]+ALGEBRAIC[Istim]) -  2.00000*ALGEBRAIC[INaK])*CONSTANTS[cm]*CONSTANTS[Acap])/( CONSTANTS[F]*CONSTANTS[vmyo])+( ALGEBRAIC[JdiffK]*CONSTANTS[vss])/CONSTANTS[vmyo];
RATES[kss] = ( - ALGEBRAIC[ICaK]*CONSTANTS[cm]*CONSTANTS[Acap])/( CONSTANTS[F]*CONSTANTS[vss]) - ALGEBRAIC[JdiffK];
#endif
#if defined NA_ORUDY2011 
RATES[m] = (ALGEBRAIC[mss] - STATES[m])/ALGEBRAIC[tm];
RATES[j] = (ALGEBRAIC[jss] - STATES[j])/ALGEBRAIC[tj];
RATES[jp] = (ALGEBRAIC[jss] - STATES[jp])/ALGEBRAIC[tjp];
RATES[nai] = ( - (ALGEBRAIC[INa]+ALGEBRAIC[INaL]+ 3.00000*ALGEBRAIC[INaCa_i]+ 3.00000*ALGEBRAIC[INaK]+ALGEBRAIC[INab])*CONSTANTS[Acap]*CONSTANTS[cm])/( CONSTANTS[F]*CONSTANTS[vmyo])+( ALGEBRAIC[JdiffNa]*CONSTANTS[vss])/CONSTANTS[vmyo];
RATES[hf] = (ALGEBRAIC[hss] - STATES[hf])/ALGEBRAIC[thf];
RATES[hs] = (ALGEBRAIC[hss] - STATES[hs])/ALGEBRAIC[ths];
RATES[hsp] = (ALGEBRAIC[hssp] - STATES[hsp])/ALGEBRAIC[thsp];
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
RATES[ki] = ( - ((ALGEBRAIC[Ito]+ALGEBRAIC[IKr]+ALGEBRAIC[IKs]+ALGEBRAIC[IK1]+ALGEBRAIC[IKb]+ALGEBRAIC[Istim]) -  2.00000*ALGEBRAIC[INaK])*CONSTANTS[cm]*CONSTANTS[Acap])/( CONSTANTS[F]*CONSTANTS[vmyo])+( ALGEBRAIC[JdiffK]*CONSTANTS[vss])/CONSTANTS[vmyo];
RATES[iF] = (ALGEBRAIC[iss] - STATES[iF])/ALGEBRAIC[tiF];
RATES[iS] = (ALGEBRAIC[iss] - STATES[iS])/ALGEBRAIC[tiS];
RATES[kss] = ( - ALGEBRAIC[ICaK]*CONSTANTS[cm]*CONSTANTS[Acap])/( CONSTANTS[F]*CONSTANTS[vss]) - ALGEBRAIC[JdiffK];
RATES[CaMKt] =  CONSTANTS[aCaMK]*ALGEBRAIC[CaMKb]*(ALGEBRAIC[CaMKb]+STATES[CaMKt]) -  CONSTANTS[bCaMK]*STATES[CaMKt];
#endif
}
