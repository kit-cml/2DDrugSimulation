#ifndef TEN_TUSSCHER_MODEL_2006_IK1KO_EPI_UNITS_HPP
#define TEN_TUSSCHER_MODEL_2006_IK1KO_EPI_UNITS_HPP

#include "patch_clamp.hpp"
#include "../enums/enum_ten_tusscher_model_2006.hpp"

/**
 *  The CellML child class that contains the information from ten_tusscher_model_2006_IK1Ko_epi_units.cellml file.
 *
 *  @see CellML
 */
class ten_tusscher_model_2006_IK1Ko_epi_units : public patch_clamp
{
public:
  ten_tusscher_model_2006_IK1Ko_epi_units();
  ~ten_tusscher_model_2006_IK1Ko_epi_units();
  void initConsts ();
  void computeRates( double TIME, double *CONSTANTS, double *RATES, double *STATES, double *ALGEBRAIC );
  void solveAnalytical( double dt );
};


#endif

