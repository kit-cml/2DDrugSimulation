#ifndef COURTEMANCHE_RAMIREZ_NATTEL_1998_HPP
#define COURTEMANCHE_RAMIREZ_NATTEL_1998_HPP

#include "cellml.hpp"
#include "../enums/enum_courtemanche_ramirez_nattel_1998.hpp"

/**
 *  The CellML child class that contains the information from courtemanche_ramirez_nattel_1998.cellml file.
 *
 *  @see CellML
 */
class courtemanche_ramirez_nattel_1998 : public CellML
{
public:
  courtemanche_ramirez_nattel_1998();
  ~courtemanche_ramirez_nattel_1998();
  void initConsts ();
  void computeRates( double TIME, double *CONSTANTS, double *RATES, double *STATES, double *ALGEBRAIC );
  void solveAnalytical( double dt );
};


#endif
