#ifndef OHARA_RUDY_CIPA_V1_2017_HPP
#define OHARA_RUDY_CIPA_V1_2017_HPP

#include "cellml.hpp"
#include "../enums/enum_ohara_rudy_cipa_v1_2017.hpp"

/**
 *  The CellML child class that contains the information from ohara_rudy_cipa_v1_2017.cellml file.
 *
 *  @see CellML
 */
class ohara_rudy_cipa_v1_2017 : public CellML
{
public:
  ohara_rudy_cipa_v1_2017();
  ~ohara_rudy_cipa_v1_2017();
  void initConsts ();
  void initConsts (double type);
  void initConsts (double type, double conc, double *herg, const char *drug_name);
  void computeRates( double TIME, double *CONSTANTS, double *RATES, double *STATES, double *ALGEBRAIC );
  void solveAnalytical( double dt );
private:
  void ___initConsts ();
  void ___initStates_bepridil ();
};


#endif

