#ifndef OHARA_RUDY_2011_HPP
#define OHARA_RUDY_2011_HPP

#include "cellml.hpp"
#include "../enums/enum_Ohara_Rudy_2011.hpp"

/**
 *  The CellML child class that contains the information from Ohara_Rudy_2011.cellml file.
 *
 *  @see CellML
 */
class Ohara_Rudy_2011 : public CellML
{
public:
  Ohara_Rudy_2011();
  ~Ohara_Rudy_2011();
  void initConsts ();
  void initConsts (double type);
  void initConsts (double type, double D, double *hill, bool is_dutta, const char *drug_name);
  void computeRates( double TIME, double *CONSTANTS, double *RATES, double *STATES, double *ALGEBRAIC );
  void solveAnalytical( double dt );
private:
  void ___initConsts ();
  void ___initStates_bepridil (); 
};


#endif
