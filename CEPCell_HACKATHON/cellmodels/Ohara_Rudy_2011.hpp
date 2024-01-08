#ifndef OHARA_RUDY_2011_HPP
#define OHARA_RUDY_2011_HPP

#include "patch_clamp.hpp"
#include "../enums/enum_Ohara_Rudy_2011.hpp"

class Ohara_Rudy_2011 : public patch_clamp
{
public:
  Ohara_Rudy_2011();
  ~Ohara_Rudy_2011();
  void ___applyDutta();
  void ___initConsts();
  void initConsts ();
  void initConsts (double type);
  void initConsts (double type, bool is_dutta);
  void initConsts (double type, double D, double *hill, bool is_dutta);
  void computeRates( double TIME, double *CONSTANTS, double *RATES, double *STATES, double *ALGEBRAIC );
  void solveAnalytical(double dt);
private:
};


#endif

