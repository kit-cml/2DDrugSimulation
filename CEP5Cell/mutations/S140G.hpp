#ifndef S140G_HPP
#define S140G_HPP

#include "mutation.hpp"

class S140G : public Mutation{

  private:
    double phi;
  public:
    S140G( double phi );
    void mutate(double* ALGEBRAIC, double* STATES, double* CONSTANTS);
};



#endif
