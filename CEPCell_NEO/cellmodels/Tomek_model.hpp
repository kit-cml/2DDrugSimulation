#ifndef TOMEK_MODEL_HPP
#define TOMEK_MODEL_HPP

#include "patch_clamp.hpp"
#include "../enums/enum_Tomek_model.hpp"


class Tomek_model : public patch_clamp
{
private:
  void ___initEndo();
  void ___initEpi();
  void ___initMid();
public:
  Tomek_model();
  ~Tomek_model();
  void initConsts ();
  void initConsts (double celltype);
  void initConsts (double celltype, double D, double *hill);
  void computeRates( double TIME, double *CONSTANTS, double *RATES, double *STATES, double *ALGEBRAIC );
  void solveAnalytical( double dt );
};


#endif

