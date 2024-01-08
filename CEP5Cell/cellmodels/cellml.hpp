#ifndef CELLML_HPP
#define CELLML_HPP

#include "../mutations/mutation.hpp"


class CellML 
{
protected:
  CellML() {}
public:
  int isMutated;
  int isEctopic;
  int isS1;
  int is_hill;
  unsigned int algebraic_size;
  unsigned int constants_size;
  unsigned int states_size;
  double *ALGEBRAIC;
  double *CONSTANTS;
  double *RATES;
  double *STATES;
  double *STATES_INIT;
  Mutation *mutation;
  virtual ~CellML() {};
  virtual void initConsts(){}
  virtual void initConsts(double type){}
  virtual void initConsts (double type, double conc, double *hill, const char *drug_name){}
  virtual void initConsts(double type, double conc, double *hill, bool is_dutta, const char *drug_name){}
  virtual void computeRates(double TIME, double *CONSTANTS, double *RATES, double *STATES, double *ALGEBRAIC) = 0;
  virtual void solveAnalytical(double dt) = 0;
};

#endif
