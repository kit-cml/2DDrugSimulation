#include "bspm.h"
#include <omp.h>
#include <pthread.h>

Cbspm::Cbspm(Cmesh* pmesh): Cfem(pmesh)
{
  fpecg = fopen("plotECG.plt", "wt");

  double Dpar = DC*pm.sf_diff;
  double Dnor = Dpar; //Dpar/7.0;
  double delD = Dnor-Dpar;
  for(int i=0; i<pmesh->dimension; i++){
    for(int j=0; j<pmesh->dimension; j++)
      DN[i][j] = delD;  //= delD*fibori[i]*fibori[j];
    DN[i][i] += Dpar;
  }

}
/*
Cbspm::Cbspm(Cmesh* p_meshsource, Cmesh* p_mesh): Cfem(p_meshsource, p_mesh)
{
  DTensor = new double* [3];
  for(int i=0; i<3; i++)
    DTensor[i] = new double [3];
  time =0.;
  fpecg = fopen("plotECG.plt", "wt");
}*/


bool Cbspm::writeBSPM(double time)
{
  //sprintf(fname,"./bspm/Torso%d.vtk",(int)floor(time+0.5));
  sprintf(fname,"./Torso%d.vtk",(int)time);
  fbspm=fopen(fname,"wt");

  fprintf(fbspm, "# vtk DataFile Version 3.0\n");
  fprintf(fbspm, "vtk output\n");
  fprintf(fbspm, "ASCII\n");
  fprintf(fbspm, "DATASET UNSTRUCTURED_GRID\n");

  fprintf(fbspm, "POINTS %d float\n", pmesh->all_nodes);
  for(int i=0; i<pmesh->all_nodes; i++)
    fprintf(fbspm, "%lf %lf %lf\n", pmesh->XY[i][x], pmesh->XY[i][y], pmesh->XY[i][z]);

  fprintf(fbspm, "CELLS %d %d\n", pmesh->all_elements, pmesh->all_elements*4);
  for(int i=0; i<pmesh->all_elements; i++)
    fprintf(fbspm, "3 %d %d %d\n", pmesh->EC[i][0], pmesh->EC[i][1], pmesh->EC[i][2]);

  fprintf(fbspm, "CELL_TYPES %d\n", pmesh->all_elements);
  for(int i=0; i<pmesh->all_elements; i++)
    fprintf(fbspm, "5\n");  // 5 : triangle

  fprintf(fbspm, "POINT_DATA %d\n", pmesh->all_nodes);
  fprintf(fbspm, "SCALARS Voltage float\n");
  fprintf(fbspm, "LOOKUP_TABLE default\n");
  for(int i=0; i<pmesh->all_nodes; i++)
    fprintf(fbspm, "%lf\n", volt[i]);  // Vm
  fclose(fbspm);
  return true;
}


bool Cbspm::writeECG(double time)
{
  PetscPrintf(PETSC_COMM_WORLD, "writeecg\n");
  double Ewct = (ecgl.vnode[0]+ecgl.vnode[1]+ecgl.vnode[2])/3.0;
  ecgl.xL[0] = ecgl.vnode[1]-ecgl.vnode[0];
  ecgl.xL[1] = ecgl.vnode[2]-ecgl.vnode[0];
  ecgl.xL[2] = ecgl.vnode[2]-ecgl.vnode[1];
  ecgl.av[0] = 1.5*(ecgl.vnode[0]-Ewct);
  ecgl.av[1] = 1.5*(ecgl.vnode[1]-Ewct);
  ecgl.av[2] = 1.5*(ecgl.vnode[2]-Ewct);
  ecgl.precordial[0] = ecgl.vnode[3]-Ewct;
  ecgl.precordial[1] = ecgl.vnode[4]-Ewct;
  ecgl.precordial[2] = ecgl.vnode[5]-Ewct;
  ecgl.precordial[3] = ecgl.vnode[6]-Ewct;
  ecgl.precordial[4] = ecgl.vnode[7]-Ewct;
  ecgl.precordial[5] = ecgl.vnode[8]-Ewct;

  fprintf(fpecg, "%lf ", time);
  for(int i=0; i<3; i++) fprintf(fpecg, "%lf ", ecgl.xL[i]);
  for(int i=0; i<3; i++) fprintf(fpecg, "%lf ", ecgl.av[i]);
  for(int i=0; i<6; i++) fprintf(fpecg, "%lf ", ecgl.precordial[i]);
  fprintf(fpecg, "\n");

  return true;
}
int Cbspm::processing(int argc,char **args, param_t pm) // override
{
  return 1;
}

void Cbspm::writeVTK() // override
{

}

void Cbspm::writeoutput() // override
{

}

