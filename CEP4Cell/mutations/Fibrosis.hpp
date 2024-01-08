#ifndef FIBROSIS_HPP
#define FIBROSIS_HPP

#include "mutation.hpp"

class Fibrosis : public Mutation{

  public:
    Fibrosis();
    void mutate(double* ALGEBRAIC, double* STATES, double* CONSTANTS);
};


#endif
