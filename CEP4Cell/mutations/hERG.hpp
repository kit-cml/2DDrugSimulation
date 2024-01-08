#ifndef L532P_HPP
#define L532P_HPP

#include "mutation.hpp"

class hERG : public Mutation{

  private:
    int type;
  public:
    hERG(const char* str);
    void mutate(double* ALGEBRAIC, double* STATES, double* CONSTANTS);
};



#endif
