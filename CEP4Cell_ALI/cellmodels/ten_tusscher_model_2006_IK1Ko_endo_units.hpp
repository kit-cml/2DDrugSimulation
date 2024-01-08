#ifndef TEN_TUSSCHER_MODEL_2006_IK1KO_ENDO_UNITS_HPP
#define TEN_TUSSCHER_MODEL_2006_IK1KO_ENDO_UNITS_HPP

#include "cellml.hpp"
#include "../enums/enum_ten_tusscher_model_2006_IK1Ko.hpp"

/**
 *  The CellML child class that contains the information from ten_tusscher_model_2006_IK1Ko_endo_units.cellml file.
 *
 *  @see CellML
 */
class ten_tusscher_model_2006_IK1Ko_endo_units : public CellML
{
  public:
    ten_tusscher_model_2006_IK1Ko_endo_units();
    ~ten_tusscher_model_2006_IK1Ko_endo_units();
    void initConsts ();
    void computeRates( double TIME, double *CONSTANTS, double *RATES, double *STATES, double *ALGEBRAIC );
    void solveAnalytical( double dt );
};


#endif

