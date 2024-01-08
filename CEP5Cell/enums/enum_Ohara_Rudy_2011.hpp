#ifndef EN_OHARA_RUDY_2011_HPP
#define EN_OHARA_RUDY_2011_HPP

enum E_ALGEBRAIC_T{
  vffrt = 29,
  vfrt = 39,
  i_Na = 62,
  i_NaL = 64,
  i_to = 70,
  i_CaL = 81,
  i_CaNa = 82,
  i_CaK = 85,
  i_Kr = 96,
  i_Ks = 98,
  i_K1 = 100,
  i_NaCa_i = 132,
  i_NaCa_ss = 162,
  i_NaK = 181,
  i_Nab = 184,
  i_Kb = 183,
  i_pCa = 188,
  i_Cab = 186,
  i_stim = 12,
  CaMKb = 45,
  CaMKa = 47,
  JdiffNa = 187,
  Jdiff = 189,
  Jup = 196,
  JdiffK = 185,
  Jrel = 191,
  Jtr = 197,
  Bcai = 49,
  Bcajsr = 53,
  Bcass = 51,
  ENa = 56,
  EK = 57,
  EKs = 58,
  mss = 0,
  tm = 13,
  hss = 1,
  thf = 14,
  ths = 15,
  h = 59,
  jss = 16,
  tj = 30,
  hssp = 31,
  thsp = 40,
  hp = 60,
  tjp = 41,
  fINap = 61,
  mLss = 32,
  tmL = 42,
  hLss = 2,
  hLssp = 3,
  fINaLp = 63,
  ass = 4,
  ta = 17,
  iss = 5,
  delta_epi = 18,
  tiF_b = 33,
  tiS_b = 43,
  tiF = 46,
  tiS = 48,
  AiF = 65,
  AiS = 66,
  i = 67,
  assp = 34,
  dti_develop = 50,
  dti_recover = 52,
  tiFp = 54,
  tiSp = 55,
  ip = 68,
  fItop = 69,
  dss = 6,
  fss = 7,
  f = 71,
  fcass = 19,
  Afcaf = 72,
  Afcas = 73,
  fca = 74,
  fp = 75,
  fcap = 76,
  km2n = 8,
  anca = 20,
  PhiCaL = 77,
  PhiCaNa = 78,
  PhiCaK = 79,
  fICaLp = 80,
  td = 21,
  tff = 22,
  tfs = 23,
  tfcaf = 35,
  tfcas = 36,
  tffp = 37,
  tfcafp = 44,
  xrss = 9,
  txrf = 24,
  txrs = 25,
  Axrf = 88,
  Axrs = 91,
  xr = 94,
  rkr = 95,
  xs1ss = 10,
  xs2ss = 26,
  txs1 = 27,
  KsCa = 97,
  txs2 = 38,
  xk1ss = 11,
  txk1 = 28,
  rk1 = 99,
  hna = 102,
  hca = 101,
  h1_i = 103,
  h2_i = 104,
  h3_i = 105,
  h4_i = 106,
  h5_i = 107,
  h6_i = 108,
  h7_i = 109,
  h8_i = 110,
  h9_i = 111,
  k3p_i = 112,
  k3pp_i = 113,
  k3_i = 114,
  k4_i = 117,
  k4p_i = 115,
  k4pp_i = 116,
  k6_i = 118,
  k7_i = 119,
  k8_i = 120,
  x1_i = 121,
  x2_i = 122,
  x3_i = 123,
  x4_i = 124,
  E1_i = 125,
  E2_i = 126,
  E3_i = 127,
  E4_i = 128,
  allo_i = 129,
  JncxNa_i = 130,
  JncxCa_i = 131,
  h1_ss = 133,
  h2_ss = 134,
  h3_ss = 135,
  h4_ss = 136,
  h5_ss = 137,
  h6_ss = 138,
  h7_ss = 139,
  h8_ss = 140,
  h9_ss = 141,
  k3p_ss = 142,
  k3pp_ss = 143,
  k3_ss = 144,
  k4_ss = 147,
  k4p_ss = 145,
  k4pp_ss = 146,
  k6_ss = 148,
  k7_ss = 149,
  k8_ss = 150,
  x1_ss = 151,
  x2_ss = 152,
  x3_ss = 153,
  x4_ss = 154,
  E1_ss = 155,
  E2_ss = 156,
  E3_ss = 157,
  E4_ss = 158,
  allo_ss = 159,
  JncxNa_ss = 160,
  JncxCa_ss = 161,
  Knai = 163,
  Knao = 164,
  P = 165,
  a1 = 166,
  b2 = 167,
  a3 = 168,
  b3 = 169,
  b4 = 170,
  x1 = 171,
  x2 = 172,
  x3 = 173,
  x4 = 174,
  E1 = 175,
  E2 = 176,
  E3 = 177,
  E4 = 178,
  JnakNa = 179,
  JnakK = 180,
  xkb = 182,
  Jrel_inf = 86,
  tau_rel = 92,
  Jrel_infp = 87,
  Jrel_temp = 84,
  tau_relp = 93,
  Jrel_inf_temp = 83,
  fJrelp = 190,
  tau_rel_temp = 89,
  tau_relp_temp = 90,
  Jupnp = 192,
  Jupp = 193,
  fJupp = 194,
  Jleak = 195,
};

enum E_CONSTANTS_T{
  celltype = 0,
  nao = 1,
  cao = 2,
  ko = 3,
  R = 4,
  T = 5,
  F = 6,
  zna = 7,
  zca = 8,
  zk = 9,
  L = 10,
  rad = 11,
  vcell = 104,
  Ageo = 113,
  Acap = 117,
  vmyo = 118,
  vnsr = 119,
  vjsr = 120,
  vss = 121,
  stim_amplitude = 12,
  stim_duration = 13,
  KmCaMK = 14,
  aCaMK = 15,
  bCaMK = 16,
  CaMKo = 17,
  KmCaM = 18,
  cmdnmax_b = 19,
  cmdnmax = 93,
  kmcmdn = 20,
  trpnmax = 21,
  kmtrpn = 22,
  BSRmax = 23,
  KmBSR = 24,
  BSLmax = 25,
  KmBSL = 26,
  csqnmax = 27,
  kmcsqn = 28,
  cm = 29,
  PKNa = 30,
  mssV1 = 31,
  mssV2 = 32,
  mtV1 = 33,
  mtV2 = 34,
  mtD1 = 35,
  mtD2 = 36,
  mtV3 = 37,
  mtV4 = 38,
  hssV1 = 39,
  hssV2 = 40,
  Ahs = 94,
  Ahf = 41,
  GNa = 42,
  thL = 43,
  thLp = 95,
  GNaL_b = 44,
  GNaL = 96,
  Gto_b = 45,
  Gto = 97,
  Kmn = 46,
  k2n = 47,
  PCa_b = 48,
  Aff = 98,
  Afs = 109,
  PCa = 99,
  PCap = 110,
  PCaNa = 111,
  PCaK = 112,
  PCaNap = 115,
  PCaKp = 116,
  tjca = 100,
  GKr_b = 49,
  GKr = 101,
  GKs_b = 50,
  GKs = 102,
  GK1 = 103,
  GK1_b = 51,
  kna1 = 52,
  kna2 = 53,
  kna3 = 54,
  kasymm = 55,
  wna = 56,
  wca = 57,
  wnaca = 58,
  kcaon = 59,
  kcaoff = 60,
  qna = 61,
  qca = 62,
  KmCaAct = 63,
  Gncx_b = 64,
  Gncx = 128,
  h10_i = 122,
  h11_i = 123,
  h12_i = 124,
  k1_i = 125,
  k2_i = 126,
  k5_i = 127,
  h10_ss = 129,
  h11_ss = 130,
  h12_ss = 131,
  k1_ss = 132,
  k2_ss = 133,
  k5_ss = 134,
  k1p = 65,
  k1m = 66,
  k2p = 67,
  k2m = 68,
  k3p = 69,
  k3m = 70,
  k4p = 71,
  k4m = 72,
  Knai0 = 73,
  Knao0 = 74,
  delta = 75,
  Kki = 76,
  Kko = 77,
  MgADP = 78,
  MgATP = 79,
  Kmgatp = 80,
  H = 81,
  eP = 82,
  Khp = 83,
  Knap = 84,
  Kxkur = 85,
  Pnak_b = 86,
  Pnak = 138,
  b1 = 135,
  a2 = 136,
  a4 = 137,
  GKb_b = 87,
  GKb = 105,
  PNab = 88,
  PCab = 89,
  GpCa = 90,
  KmCap = 91,
  bt = 92,
  a_rel = 106,
  btp = 107,
  a_relp = 114,
  upScale = 108,
  stim_start = 139,
  stim_period = 140,
};

enum E_STATES_T{
  V = 0,
  CaMKt = 1,
  cass = 2,
  nai = 3,
  nass = 4,
  ki = 5,
  kss = 6,
  cansr = 7,
  cajsr = 8,
  Ca_i = 9,
  m = 10,
  hf = 11,
  hs = 12,
  j = 13,
  hsp = 14,
  jp = 15,
  mL = 16,
  hL = 17,
  hLp = 18,
  a = 19,
  iF = 20,
  iS = 21,
  ap = 22,
  iFp = 23,
  iSp = 24,
  d = 25,
  ff = 26,
  fs = 27,
  fcaf = 28,
  fcas = 29,
  jca = 30,
  ffp = 31,
  fcafp = 32,
  nca = 33,
  xrf = 34,
  xrs = 35,
  xs1 = 36,
  xs2 = 37,
  xk1 = 38,
  Jrelnp = 39,
  Jrelp = 40,
  qnet = 41,
  INaL_AUC = 42,
  ICaL_AUC = 43
};

#endif
