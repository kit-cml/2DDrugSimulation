#ifndef G229D_HPP
#define G229D_HPP

#include "mutation.hpp"

class G229D : public Mutation{

  private:
    bool is_WT;
  public:
    G229D(const char *str);
    void mutate(double* ALGEBRAIC, double* STATES, double* CONSTANTS);
};



#endif
