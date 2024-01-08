#ifndef TENTUSSCHER_NOBLE_NOBLE_PANFILOV_2004_B_HPP
#define TENTUSSCHER_NOBLE_NOBLE_PANFILOV_2004_B_HPP

#include "patch_clamp.hpp"
#include "../enums/enum_tentusscher_noble_noble_panfilov_2004.hpp"

#define NA_TN2004
//#define NA_A1656D_TN2004
#define BNA_TN2004
#define NAK_TN2004
#define NACA_TN2004
#define KR_TN2004
#define KS_TN2004
#define K1_TN2004
#define TO_TN2004
#define PK_TN2004
#define CAL_TN2004
#define BCA_TN2004
#define PCA_TN2004
#define CAUP_TN2004
#define CALEAK_TN2004
#define CAREL_TN2004

#define TISSUE

class tentusscher_noble_noble_panfilov_2004_b : public patch_clamp
{
public:
  tentusscher_noble_noble_panfilov_2004_b();
  ~tentusscher_noble_noble_panfilov_2004_b();
  void initConsts ();
  void computeRates( double TIME, double *CONSTANTS, double *RATES, double *STATES, double *ALGEBRAIC );
  void solveAnalytical(double dt);
};


#endif

