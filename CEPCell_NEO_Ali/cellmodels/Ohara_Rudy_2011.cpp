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
