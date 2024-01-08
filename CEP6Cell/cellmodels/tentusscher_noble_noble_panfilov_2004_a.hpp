#ifndef TENTUSSCHER_NOBLE_NOBLE_PANFILOV_2004_A_HPP
#define TENTUSSCHER_NOBLE_NOBLE_PANFILOV_2004_A_HPP

#include "cellml.hpp"
#include "../enums/enum_tentusscher_noble_noble_panfilov_2004.hpp"

/**
 *  The CellML child class that contains the information from tentusscher_noble_noble_panfilov_2004_a.cellml file.
 *
 *  @see CellML
 */
class tentusscher_noble_noble_panfilov_2004_a : public CellML
{
public:
  tentusscher_noble_noble_panfilov_2004_a();
  ~tentusscher_noble_noble_panfilov_2004_a();
  void initConsts ();
  void computeRates( double TIME, double *CONSTANTS, double *RATES, double *STATES, double *ALGEBRAIC );
  void solveAnalytical( double dt );
};


#endif

