#ifndef BSPM_H
#define BSPM_H
#include "fem.h"

class Cbspm : public Cfem
{
private:
  static const int DIM = 3;
  double gpt[4][3];
  FILE* fpecg, *fbspm;
  int iter;
  char fname[100];
  FILE* fpvtk;

public:
  ECGLEADS ecgl;
  Cbspm(Cmesh* pmeshtorso);
  //Cbspm(Cmesh* p_meshsource, Cmesh* p_mesh);
  double volt[10000];
  int processing(int argc,char **args, param_t pm); // override
  void writeVTK(); // override
  void writeoutput(); // override

  bool writeECG(double time);
  bool writeBSPM(double time);
};

#endif
