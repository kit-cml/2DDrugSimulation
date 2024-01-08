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

#if defined(NA_ORUDY2011) && defined(NAL_ORUDY2011) && defined(TO_ORUDY2011) && defined(CAL_ORUDY2011) && defined(CANA_ORUDY2011) && defined(CAK_ORUDY2011) && defined(KR_ORUDY2011) && defined(KS_ORUDY2011) && defined(K1_ORUDY2011) && defined(NACA_I_ORUDY2011) && defined(NACA_SS_ORUDY2011) && defined(NAK_ORUDY2011) && defined(NAB_ORUDY2011) && defined(KB_ORUDY2011) && defined(PCA_ORUDY2011) && defined(CAB_ORUDY2011)
  #define SINGLE_CELL
#endif

class Ohara_Rudy_2011 : public patch_clamp
{
public:
  Ohara_Rudy_2011();
  ~Ohara_Rudy_2011();
  void ___initConsts();
  void initConsts ();
  void initConsts (double type);
  void initConsts (double type, double D, double *hill, bool is_dutta);
  void computeRates( double TIME, double *CONSTANTS, double *RATES, double *STATES, double *ALGEBRAIC );
private:
};


#endif

