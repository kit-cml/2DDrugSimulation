#ifndef MUTATION_HPP
#define MUTATION_HPP

#include "../commons/usermacro.hpp"

class Mutation 
{
protected:
  Mutation(){};
public:
  virtual void mutate(double* ALGEBRAIC, double* STATES, double* CONSTANTS) = 0;
};

#endif
