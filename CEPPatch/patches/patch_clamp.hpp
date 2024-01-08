#ifndef PATCH_CLAMP_HPP
#define PATCH_CLAMP_HPP

class patch_clamp
{
protected:
  patch_clamp(){}
public:
  unsigned int algebraic_size;
  unsigned int constants_size;
  unsigned int states_size;
  double *ALGEBRAIC;
  double *CONSTANTS;
  double *RATES;
  double *STATES;
  virtual void initConsts() = 0;
  virtual void initConsts(double type){}
  virtual void initConsts (double type, double conc, double *hill){}
  virtual void computeRates(double TIME, double *CONSTANTS, double *RATES, double *STATES, double *ALGEBRAIC) = 0;
};

#endif
