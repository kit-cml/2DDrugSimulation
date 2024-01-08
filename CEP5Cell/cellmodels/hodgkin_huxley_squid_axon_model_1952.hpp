#ifndef HODGKIN_HUXLEY_SQUID_AXON_MODEL_1952_HPP
#define HODGKIN_HUXLEY_SQUID_AXON_MODEL_1952_HPP

#include "cellml.hpp"
#include "../enums/enum_hodgkin_huxley_squid_axon_model_1952.hpp"

/**
 *  The CellML child class that contains the information from hodgkin_huxley_squid_axon_model_1952.cellml file.
 *
 *  @see CellML
 */
class hodgkin_huxley_squid_axon_model_1952 : public CellML
{
public:
  hodgkin_huxley_squid_axon_model_1952();
  ~hodgkin_huxley_squid_axon_model_1952();
  void initConsts ();
  void computeRates( double TIME, double *CONSTANTS, double *RATES, double *STATES, double *ALGEBRAIC );
  void solveAnalytical( double dt );
};


#endif

