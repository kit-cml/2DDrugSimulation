#ifndef AT_FIB_HPP
#define AT_FIB_HPP

#include "mutation.hpp"

class At_Fib : public Mutation{

  public:
    At_Fib();
    void mutate( double *ALGEBRAIC, double *STATES, double *CONSTANTS );
};



#endif
