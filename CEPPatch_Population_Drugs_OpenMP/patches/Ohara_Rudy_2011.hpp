#ifndef OHARA_RUDY_2011_HPP
#define OHARA_RUDY_2011_HPP

#include "patch_clamp.hpp"
#include "../enums/enum_Ohara_Rudy_2011.hpp"

#define NA_ORUDY2011
#define NAL_ORUDY2011
#define TO_ORUDY2011
#define CAL_ORUDY2011
#define CANA_ORUDY2011
#define CAK_ORUDY2011
#define KR_ORUDY2011
#define KS_ORUDY2011
#define K1_ORUDY2011
#define NACA_I_ORUDY2011
#define NACA_SS_ORUDY2011
#define NAK_ORUDY2011
#define NAB_ORUDY2011
#define KB_ORUDY2011
#define PCA_ORUDY2011
#define CAB_ORUDY2011

class Ohara_Rudy_2011 : public patch_clamp
{
public:
  Ohara_Rudy_2011();
  ~Ohara_Rudy_2011();
  void initConsts ();
  void initConsts (double type);
  void initConsts (double type, double D, double *hill, bool is_dutta);
  void computeRates( double TIME, double *CONSTANTS, double *RATES, double *STATES, double *ALGEBRAIC );
private:
  void ___initConsts();
};


#endif

