#include "A1656D.hpp"
#include "../enums/enum_tentusscher_noble_noble_panfilov_2004.hpp"

#include <cmath>
#include <cstring>
#include <cstdio>

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
/*
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
*/

  // Inje Code
  double AdjustF[5] = { 0.5, 0.4189418, 1.333233, 2.539508, 2.559347 };
  double wt[24] = { AdjustF[0] * 520,-33.4609808234912,-5.63753034812371,-81.8951228209157,5.62771308469908,
		0.663712605742043, 146.921178619256, 29.2400346266816, -52.0, 0.207360629474796,
		2.89437515408304, 242.28779205519, 431.478085811286, -50.0, 0.327111948539544,
		183.467968356691, 169.515305828317, 1339.70833771837, -77.0, 8.42604440633291,
		0.9817 };

  const double mt[24] = { AdjustF[1]*170,-44.6935990478942,-5.36860593965297,-69.4760433770628,
		4.89788835743973,5.07788041282928,196.176009392084,85.0862531683836,-68.7643477250666,0,13.5675103597582,
		-181.415695718764,2.77155685948369,-36.689899890811,0,1.38*325.875005324949,128.854673106975,165.944408699613,-3.42290346461849,17.5510189408494,0.8627533044};

  const double mex[24] = { AdjustF[2] * 57.8,
			-52.1119427677749,-3.73592397914885,
			-76.2669236828744,6.34359640549037,
			0.80565634931842,297.47496625459,-96.6719372417807,-50.3330164735987,0.00465853471402205,
			5.27272981752212,-407.717612458511,12.7902611681218,-55.918825676663,1.13835720080851,
			1.46788942367968,-11.5330246256226,672.047842462476,-48.139905930177,56.136,
			0.99 };

  const double fle[24] = { AdjustF[3] * 34,-54.6090051693093,-3.61452208306045,-73.5332567849385,8.2200766426964,1.34682372999204,122.631638305188,53.1271812972241,-61.5364333672071,0.394163124659821,11.0209831376502,168.052880317447,35.4129546482693,53.8381598146392,0,9.1979219482718,-15.9885736099013,362.57139759009,-10.1381427423855,17.083592607715,0.5 }; 

  const double ran[24] = { AdjustF[4] * 34, -50.2700247654263, -4.95018322071719, -70.5298452364176, 4.86595516223751, 3.39378511922962, 161.293500879792, 90.4267497230428, -63.2899528911703, 0.135348276255773, 19.6068812657396, 30.8819532857551, 84.4559330615111, 14.4826774271407, 1.74483131548136, 91.0693475950254, -103.071290691799, -84.1400593760986, 6.12647640696309, 5.69707375208435, 0.88044 };


  const double RT = (CONSTANTS[R]*CONSTANTS[T])/1000.0;
  static double mtemp=STATES[m], htemp=STATES[h], jtemp=STATES[j];


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

/*
  ALGEBRAIC[m_inf] = pow( 1.0 / (1.0 + exp((STATES[V]-wt[1]) / wt[2])), 1.0/3.0 );
  ALGEBRAIC[tau_m] = wt[4] / (exp(wt[5] * (STATES[V] - wt[3]) / RT) + exp(-wt[6] * (STATES[V] - wt[3]) / RT)) + wt[7];
  ALGEBRAIC[h_inf] = 1.0 / (1.0 + exp((STATES[V] - wt[8]) / wt[9]));
  ALGEBRAIC[tau_h] = wt[11] / (exp(wt[12] * (STATES[V] - wt[10]) / RT) + exp(-wt[13] * (STATES[V] - wt[10]) / RT)) + wt[14];
  ALGEBRAIC[j_inf] = ALGEBRAIC[h_inf];
  ALGEBRAIC[tau_j] = wt[16] / (exp(wt[17] * (STATES[V] - wt[15]) / RT) + exp(-wt[18] * (STATES[V] - wt[15]) / RT)) + wt[19];
*/

  // Inje Code
  ALGEBRAIC[tau_m] = wt[5] / (exp(wt[6] * (STATES[V] - wt[8]) / RT) + exp(-wt[7] * (STATES[V] - wt[8]) / RT)) + wt[9];
  ALGEBRAIC[m_inf] = pow(1.0 / (1.0 + exp((STATES[V] - wt[1]) / wt[2])), 1.0 / 3.0);
  ALGEBRAIC[tau_h] = wt[10] / (exp(wt[11] * (STATES[V] - wt[13]) / RT) + exp(-wt[12] * (STATES[V] - wt[13]) / RT)) + wt[14];
  ALGEBRAIC[h_inf] = 1.0 / (1.0 + exp((STATES[V] - wt[3]) / wt[4]));
  ALGEBRAIC[tau_j] = wt[15] / (exp(wt[16] * (STATES[V] - wt[18]) / RT) + exp(-wt[17] * (STATES[V] - wt[18]) / RT)) + wt[19];
  ALGEBRAIC[j_inf] = ALGEBRAIC[h_inf];

/*
  mtemp = ALGEBRAIC[m_inf]-( ( ALGEBRAIC[m_inf]-0. )*exp( -0.01/ALGEBRAIC[tau_m] ) ); 
  htemp = ALGEBRAIC[h_inf]-( ( ALGEBRAIC[h_inf]-0. )*exp( -0.01/ALGEBRAIC[tau_h] ) ); 
  jtemp = ALGEBRAIC[j_inf]-( ( ALGEBRAIC[j_inf]-0. )*exp( -0.01/ALGEBRAIC[tau_j] ) );
*/
/*
  mtemp = ( ( mtemp - ALGEBRAIC[m_inf] )*exp( -0.01/ALGEBRAIC[tau_m] ) ) + ALGEBRAIC[m_inf];
  htemp = ( ( htemp - ALGEBRAIC[h_inf] )*exp( -0.01/ALGEBRAIC[tau_h] ) ) + ALGEBRAIC[h_inf]; 
  jtemp = ( ( jtemp - ALGEBRAIC[j_inf] )*exp( -0.01/ALGEBRAIC[tau_j] ) ) + ALGEBRAIC[j_inf];


  ALGEBRAIC[i_Na] = wt[0]*pow(mtemp, 3.0)*((wt[20]*htemp) + ((1.0-wt[20])*jtemp))*(STATES[V] - ALGEBRAIC[E_Na]);
*/
}
