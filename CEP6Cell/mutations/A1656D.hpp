#ifndef A1656D_HPP
#define A1656D_HPP

#include "mutation.hpp"

class A1656D : public Mutation{

  private:
    int type;
  public:
    A1656D( const char* str );
    void mutate(double* ALGEBRAIC, double* STATES, double* CONSTANTS);
};



#endif
