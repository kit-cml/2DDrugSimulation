#ifndef FEM_H
#define FEM_H 

#include <vector>
#include <queue>

#include "CEPCell_NEO/modules/param.hpp"
#include "mesh.h"

using std::vector;

static char help[] = "Solves a linear system in parallel with KSP.  Also\n\
illustrates setting a user-defined shell preconditioner and using the\n\
macro __FUNCT__ to define routine names for use in error handling.\n\
Input parameters include:\n\
  -user_defined_pc : Activate a user-defined preconditioner\n\n";

#define NVAR 5  // number of variables for results writing
#define NUM_CELLPARAM 60
class Cfem
{
protected:
  param_t pm;
  double time;
  double t_out;
  double** vd; //[NODE][1]  
  double** vm_ele; //[NODE][1]  
  double** vel; //[NODE][1]  
  double** gpt; //[NODE][DIM]  gauss points
  double** Phi; //[NODE][NODE]  shape function
  double*** dPhiX; //[NODE][NODE][DIM] gradient of shape function in real coordinate
  double*** dPhiXt; // [NODE][NODE][DIM] Transpose of dPhiX
  double*** dPhiX2; // [NODE][NODE][DIM]
  double*** Jacob; //[NODE][DIM][DIM] jacobian
  double*** xJ; //[NODE][DIM][DIM] inverse of jacobian
  double*** xJt; //[NODE][DIM][DIM] Transpose of xJ
  double*** dPhiS; //[NODE][NODE][DIM] gradient of shape function in local coordinate
  double*** dPhiSt; //[NODE][NODE][DIM] gradient of shape function in local coordinate
  double** xyele;  // [NODE][DIM] local coordinate
  double* detJ; // [NODE] determinant of jacobian
  double** amat; // local stiffness matrix
  double** bmat; // local lumped matrix

  Mat Amat;  // global stiffness matrix
  Mat Bmat;  // global lumped matrix
  Vec bb;
  Vec Volt[NVAR]; // local vectors for saving state variables (Vm, Cai, etc)
  Vec Volt1[NVAR]; // local vectors for saving algebraic variables (Ina, IKs, etc)
  PetscScalar VolumeConductor; // variable for saving cardiac potential from all nodes. 
  PetscScalar* Voltsc[NVAR]; // global vectors for saving state variables
  PetscScalar* Voltsc1[NVAR]; // global vectors for saving state variables

  vector<string> vtkout_states_vec; // string vector for states name
  vector<string> vtkout_algebraic_vec; // string vector for algebraic name

  Vec V_SEQ[NVAR], V_SEQ1[NVAR];
  VecScatter ctx[NVAR], ctx1[NVAR];

  KSP ksp;
  int* element_proc; // [all_elements]
  int* node_proc; // [all_nodes]
  int* local_nodes; // [all_nodes]
  bool* isFibroticnode; 
  //vector<bool> isFibroticVec;

public :
  Cmesh* pmesh, *pmesh_surf, *pmesh_endosurf;
  double** DN; //sigma_in; // [all_nodes][2] difussion tensor
  double* DE; //sigma_ele; // [all_elements] difussion tensor
  double*** Dtensor; // [dimension][dimension] difussion tensor
  const static double DC = 0.00154; // difussion coefficient
//  const static double DC = 0.00154*4; // difussion coefficient


  virtual int processing(int argc,char **args, param_t pm)=0;
  virtual void writeVTK()=0;
  virtual void writeoutput()=0;
  virtual void openVTK(FILE* fp);

  Cfem(Cmesh* pmesh);
  void MemAllocation();
  void prepare_Amat_Bmat_bb_Voltsc(PetscScalar state0);
  void assemble_Amat_Bmat(double dt);

  void set_gpt();
  void set_gpt_tet();
  void set_gpt_hex();
  void set_gpt_square();
  void set_gpt_line();
  void set_gpt_prism();
  void set_Phi_dPhiS();

  void inverseMat(double** Ja, double** xJ, double* detJ);
  bool MultMatVec(double** D, double* vd, double* pvel);
  void MultMatMat(int row, int colrow, int col, double** amat, double** bmat, double** cmat);
  void set_xyele(int element, double** xyele);

  void set_Phi_dPhiS_tet();
  void set_Phi_dPhiS_line();
  void set_Phi_dPhiS_square();
  void set_Phi_dPhiS_hex();
  void set_Phi_dPhiS_prism();

//  void set_Jacob(double** xyele, double*** dPhiS, double*** Jacob);
  void set_amat_bmat(int element, double dt, double** amat, double** bmat);
  void solve_Bmat(double dt);

  bool vec_cross(double* pval, double* a, double* b);
  bool vec_normalize(double* pval, double* a);
  bool vec_sub(double* pval, double* p1, double* p2);
  void cross_product3D(double** Ja, double** xJ);
  void cross_product2D(double** Ja, double** xJ);
  double det3D(double** mat);
  double det2D(double** mat);
  void assign_conductivity();
  void assign_nodes_elements_proc();
  void make_diffusion_tensor();
  void printArray(double*** array, int isize, int jsize, int ksize, char* name);
  void printArray(double** array, int isize, int jsize, char* name);
  void Transition(int row, int col, double** amat, double** bmat);
};
#endif
