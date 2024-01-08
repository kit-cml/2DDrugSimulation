#ifndef LINE_H
#define LINE_H

#include "fem.h"
#include "bspm.h"

class Cfiber : public Cfem
{

private:
  static const int REPOL = -89; // number of maximal elements
  Cbspm* pbspm;
  int* ISTIMNODE; // [NUM_NODE];
  int* ifERP; //[NUM_NODE];   // pj Vm
  int* ifastele; //[NUM_NODE];
  double speeds;
  int nde, nel, refract;
  int irate;
  double* perwave;
public:
  Cfiber(Cmesh* pmeshfiber, param_t parameter);
  Cfiber(Cmesh* pmeshfiber, Cbspm* pbspm, param_t parameter);
  FILE* fpvm;
  double* Volt; //[NUM_NODE];   // pj Vm
  double* Volt_t; //[NUM_NODE];   // i
  static const int DEPOL = 30; // number of maximal elements
  // overrider
  int processing(int argc, char **args, param_t parameter);
  void writeVTK();
  void writeoutput();

  void AVpacing();
  int PropagationPF();
  void writePF(double time);
  bool setMeshConditions();
  void openVTK(FILE* fp);

  void computeECG(); 
  void computeBSPM(); 
  double Radiate(double* elecrodei);
};

#endif
