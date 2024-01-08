#ifndef TET_DRUG_H
#define TET_DRUG_H

#include "fem.h" 
#include "bspm.h" 
#include "fiber.h"
#include "DrugSimulation/modules/glob_type.hpp"
#include "DrugSimulation/modules/cipa_t.hpp"

#ifdef TOMEK_2019
#include "DrugSimulation/cellmodels/Tomek_model.hpp"
#else
#include "DrugSimulation/cellmodels/Ohara_Rudy_2011.hpp"
#endif


#define PETSC 
//#define BSPMONLY
#define SOLVER

#include <vector>


using std::vector;

class Cheart : public Cfem
{
private :

  Cfiber* pfiber;  // purkinje model object
  Cbspm* pbspm;  // bspm model object

  double* vm_glob;
  double* initStates;  // intitial states of cell
  double** Glob_cell;
  FILE *fpvm, *fpcai, *fpapd, *fpvmpeak;
  int* mutatenode;
  bool* s1nodes; //[NUM_NODE];
  bool* s2nodes; // [all_nodes]
  bool* ectopic_nodes; // [all_nodes]
  int* s1_num;  // [all_nodes]
  int* s2_num;  // [all_nodes]
  int num_s1nodes;
  int num_s2nodes;
  int ilocal;
  int ectopic_start_time;

  int iter_count;
  int iter_bcl;
  int pace_count;
  bool idid, idid2, idid3; 

  double dt,dt_set;

  
  double qnet;
  FILE *fpqnet;
  drug_t ic50;
  double conc;
  int sample_id;
  double concs[4];
  char buffer[255];
  char drug_name[50];
  cipa_t temp_result;

  double inal_auc, ical_auc;
  double vm_repol30, vm_repol50, vm_repol90;
  double t_depol;
  double t_ca_peak, ca_amp50, ca_amp90;
  double cad50_prev, cad50_curr, cad90_prev, cad90_curr;


  bool FindSynapses();
  void GetTargetnode();
  void StimTissue();
public :
  Cheart(Cmesh* pmesh);
  Cheart(Cmesh* pmesh, Cfiber* pfiber);
  Cheart(Cmesh* pmesh, Cbspm* pbspm);
  Cheart(Cmesh* pmesh, Cfiber* pfiber, Cbspm* pbspm);
  // overrider
  int processing(int argc,char **args, param_t parameter); // override
  void preprocess();
  void cellular_reaction();
  void applyS1S2(int iter_bcl);
  void apply_ectopic();
  void apply_s2();
  void write_result(int iter_dtw_vm, int iter_dtw_output, int iter_dtw_vtk);
  void write_vtk(int iter_dtw_vm, int iter_dtw_output, int iter_dtw_vtk);
  void initialize_Glob_cell(int iopt);
  void assemble_Voltsc();

  void set_gpt_Phi_dPhiS();
  void writeVTK(); // override
  void writeVTKFeature();
  void writeVTKqnet();
  void writeScalar();
  void writeoutput(); // override
  void memAllocationHeart(); // override
  bool loadECsurf();
  bool loadPacingsites(char* pace1, char* pace2);

  void apply_drug(param_t pm);

  double Radiate(double* electrode);
  void solveBSPM(double time);
  int ifInList(int array[], int size, int n);
  drug_t get_IC50_data_from_file(const char* file_name);

};
#endif
