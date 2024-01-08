#include "fem.h"
#include "CEPCell_NEO/modules/globals.hpp"

Cfem::Cfem(Cmesh* pmesh)
{
  this->pmesh = pmesh;
  MemAllocation();
}

void Cfem::MemAllocation()
{
  int i, j, k;
  element_proc = new int [pmesh->all_elements];
  for(int i=0; i<pmesh->all_elements; i++)
    element_proc[i] = 0;

  isFibroticnode = new bool [pmesh->all_nodes];

  local_nodes = new int [pmesh->all_nodes];
  node_proc = new int [pmesh->all_nodes];

  for(int i=0; i<pmesh->all_nodes; i++){
    local_nodes[i] = 0;
    node_proc[i] = 0;
  }
  
  // Diffusion coefficient tensor (4*4)
  DE = new double [pmesh->all_elements]; // sigma_ele / element diffusion coeff
  DN = new double* [pmesh->all_nodes]; //sigma_in /node diffusion coeff
  for(i=0; i<pmesh->all_nodes; i++)
    DN[i] = new double [2];
  if(strcmp(pmesh->type, "fiber_real")){
    Dtensor = new double** [pmesh->all_elements];
    for(i=0; i<pmesh->all_elements; i++){
      Dtensor[i] = new double* [pmesh->dimension];
      for(j=0; j<pmesh->dimension; j++)
        Dtensor[i][j] = new double [pmesh->dimension];
    }
  }

  // Gaussian points (4*3)
  gpt = new double* [pmesh->node]; 
  for(i=0; i<pmesh->node; i++)  
    gpt[i] = new double [pmesh->dimension];

  // Shape function (4*4)
  Phi = new double* [pmesh->node]; 
  for(i=0; i<pmesh->node; i++)  
    Phi[i] = new double [pmesh->node];

  vd = new double* [pmesh->node];
  for(i=0; i<pmesh->node; i++)  
    vd[i] = new double [1];

  vm_ele = new double* [pmesh->node];
  for(i=0; i<pmesh->node; i++)
    vm_ele[i] = new double [1];

  vel = new double* [pmesh->node];
  for(i=0; i<pmesh->node; i++)
    vel[i] = new double [1];
  
  PetscPrintf(PETSC_COMM_WORLD, "before dPhiS %d\n", pmesh->dimension);

  // local coordiate Derivative of shape function (4*4*3)
  dPhiS= new double** [pmesh->node]; 
  for(i=0; i<pmesh->node; i++){  
    dPhiS[i] = new double* [pmesh->node];
    for(j=0; j<pmesh->node; j++){
      dPhiS[i][j] = new double [pmesh->dimension];
/*      for(k=0; k<pmesh->dimension; k++){
        dPhiS[i][j][k] = 0.;
      }
*/
    }

  }
  PetscPrintf(PETSC_COMM_WORLD, "after dPhiS\n");
  PetscPrintf(PETSC_COMM_WORLD, "before dPhiSt\n");
  dPhiSt= new double** [pmesh->node]; 
  for(i=0; i<pmesh->node; i++) {  
    dPhiSt[i] = new double* [pmesh->dimension];
    for(j=0; j<pmesh->dimension; j++) {
      dPhiSt[i][j] = new double [pmesh->node];
/*      for(k=0; k<pmesh->node; k++) {
        dPhiSt[i][j][k] = 0.;
      }
*/
    }
  }
  
  PetscPrintf(PETSC_COMM_WORLD, "before dPhiSt\n");

  // xJ: inverse matrix of Jacobian (4*3*3)
  // xJt: transpose of xJ (4*3*3)
  xJ = new double** [pmesh->node];
  xJt = new double** [pmesh->node];
  for(i=0; i<pmesh->node; i++){
    xJ[i] = new double* [pmesh->dimension];
    xJt[i] = new double* [pmesh->dimension];
    for(j=0; j<pmesh->dimension; j++){
      xJ[i][j] = new double [pmesh->dimension];
      xJt[i][j] = new double [pmesh->dimension];
    }
  }

  // determinant of Jacobian (4)
  detJ = new double [pmesh->node];
  // global coordinate derivate of shape function (4*4*3)
  dPhiX = new double** [pmesh->node];
  for(i=0; i<pmesh->node; i++){
    dPhiX[i] = new double* [pmesh->node];
    for(j=0; j<pmesh->node; j++)
      dPhiX[i][j] = new double [pmesh->dimension];
  }
  // transpose of dPhiX (4*3*4)
  dPhiXt = new double** [pmesh->node];
  for(i=0; i<pmesh->node; i++){
    dPhiXt[i] = new double* [pmesh->dimension];
    for(j=0; j<pmesh->dimension; j++){
      dPhiXt[i][j] = new double [pmesh->node];
    }
  }
  // (4*4*3) 
  dPhiX2 = new double** [pmesh->node];
  for(i=0; i<pmesh->node; i++){
    dPhiX2[i] = new double* [pmesh->node];
    for(j=0; j<pmesh->node; j++){
      dPhiX2[i][j] = new double [pmesh->dimension];
    }
  }

  xyele = new double* [pmesh->node]; //(4*3)
  for(i=0; i<pmesh->node; i++)
    xyele[i] = new double [pmesh->dimension];

  amat = new double* [pmesh->node]; // (4*4)
  bmat = new double* [pmesh->node]; // (4*4)
  for(i=0; i<pmesh->node; i++){
    amat[i] = new double [pmesh->node];
    bmat[i] = new double [pmesh->node];
  }

  Jacob = new double** [pmesh->node]; // (4*3*3)
  for(i=0; i<pmesh->node; i++){
    Jacob[i] = new double* [pmesh->dimension];
    for(j=0; j<pmesh->dimension; j++)
      Jacob[i][j] = new double [pmesh->dimension];
  }
  PetscPrintf(PETSC_COMM_WORLD, "Cfem::MemAllocation fem: %d\n", pmesh->all_nodes);
}

bool Cfem::MultMatVec(double** D, double* vd, double* pvel)
{
  pvel[x] = D[0][0]*vd[x] + D[0][1]*vd[y] + D[0][2]*vd[z];
  pvel[y] = D[1][0]*vd[x] + D[1][1]*vd[y] + D[1][2]*vd[z];
  pvel[z] = D[2][0]*vd[x] + D[2][1]*vd[y] + D[2][2]*vd[z];
  return true;
}


void Cfem::make_diffusion_tensor()
{
  for(int ele=0; ele<pmesh->all_elements; ele++){
    for(int i=0; i<pmesh->dimension; i++){   // intracellular diffusion tensor
      for(int j=0; j<pmesh->dimension; j++){
        if(i==j) Dtensor[ele][i][j]=DE[ele];
        else Dtensor[ele][i][j]=0.0;
      }
    }
  }

/*
  for(i=0; i<pmesh->dimension; i++){   // extracellular diffusion tensor
    for(j=0; j<pmesh->dimension; j++){
      if(i==j) Dextra[i][j] = Dintra[i][j]+Sige_par;
      else Dextra[i][j] = Dintra[i][j];
    }
  }*/
}

bool Cfem::vec_cross(double* pval, double* a, double* b)
{
  pval[x] = a[y]*b[z] - a[z]*b[y];
  pval[y] = a[z]*b[x] - a[x]*b[z];
  pval[z] = a[x]*b[y] - a[y]*b[x];
  return true;
}

bool Cfem::vec_normalize(double* pval, double* a)
{
  double amp = sqrt(a[x]*a[x] + a[y]*a[y] +a[z]*a[z]);

  pval[x] = a[x]/amp;
  pval[y] = a[y]/amp;
  pval[z] = a[z]/amp;
  return true;
}

bool Cfem::vec_sub(double* pval, double* p1, double* p2)
{
  pval[x] = p2[x] - p1[x];  // vector A
  pval[y] = p2[y] - p1[y];
  pval[z] = p2[z] - p1[z];
  return true;
}


void Cfem::set_gpt_tet()
{
  double aaa = (5.-sqrt(5.))/20.;     // 0.138197 : intergration point
  double bbb = (5.+3.*sqrt(5.))/20.;  // 0.58541 : integration point

  gpt[0][ss] = aaa;
  gpt[0][tt] = aaa;
  gpt[0][zz] = aaa;

  gpt[1][ss] = aaa;
  gpt[1][tt] = aaa;
  gpt[1][zz] = bbb;

  gpt[2][ss] = aaa;
  gpt[2][tt] = bbb;
  gpt[2][zz] = aaa;

  gpt[3][ss] = bbb;
  gpt[3][tt] = aaa;
  gpt[3][zz] = aaa;
}

void Cfem::set_gpt_hex()
{ 
  double aaa[2] = {-0.57735, 0.57735};
  int itx = 2;
  int ity = 2;
  int itz = 2;
  int node_number=0;
  for(int ix=0; ix<itx; ix++){
    for(int iy=0; iy<ity; iy++){
      for(int iz=0; iz<itz; iz++){
        gpt[node_number][ss] = aaa[ix];
        gpt[node_number][tt] = aaa[iy];
        gpt[node_number][zz] = aaa[iz];
        node_number++;
      }
    }
  }
}

void Cfem::set_gpt_line()
{
  gpt[0][0] = -0.57735;
  gpt[1][0] = 0.57735;
}

void Cfem::set_gpt_square()
{
  double aaa = -0.57735;
  double bbb = 0.57735;
  gpt[0][ss] = aaa;
  gpt[0][tt] = aaa;

  gpt[1][ss] = bbb;
  gpt[1][tt] = aaa;

  gpt[2][ss] = bbb;
  gpt[2][tt] = bbb;

  gpt[3][ss] = aaa;
  gpt[3][tt] = bbb;
}

void Cfem::set_gpt_prism()
{
  double a=0.445948490915965;
  double b=0.091576213509771;

  double sq3=sqrt(1./3.);

  gpt[0][0]=a;
  gpt[0][1]=a;
  gpt[0][2]=sq3;
  gpt[1][0]=1.-2.*a;
  gpt[1][1]=a;
  gpt[1][2]=sq3;
  gpt[2][0]=a;
  gpt[2][1]=1.-2.*a;
  gpt[2][2]=sq3;
  gpt[3][0]=b;
  gpt[3][1]=b;
  gpt[3][2]=sq3;
  gpt[4][0]=1.-2.*b;
  gpt[4][1]=b;
  gpt[4][2]=sq3;
  gpt[5][0]=b;
  gpt[5][1]=1.-2.*b;
  gpt[5][2]=sq3;
  gpt[6][0]=a;
  gpt[6][1]=a;
  gpt[6][2]=-sq3;
  gpt[7][0]=1.-2.*a;
  gpt[7][1]=a;
  gpt[7][2]=-sq3;
  gpt[8][0]=a;
  gpt[8][1]=1.-2.*a;
  gpt[8][2]=-sq3;
  gpt[9][0]=b;
  gpt[9][1]=b;
  gpt[9][2]=-sq3;
  gpt[10][0]=1.-2.*b;
  gpt[10][1]=b;
  gpt[10][2]=-sq3;
  gpt[11][0]=b;
  gpt[11][1]=1.-2.*b;
  gpt[11][2]=-sq3;

  gpt[0][3]=0.111690794839005;
  gpt[1][3]=0.111690794839005;
  gpt[2][3]=0.111690794839005;
  gpt[3][3]=0.054975871827661;
  gpt[4][3]=0.054975871827661;
  gpt[5][3]=0.054975871827661;
  gpt[6][3]=0.111690794839005;
  gpt[7][3]=0.111690794839005;
  gpt[8][3]=0.111690794839005;
  gpt[9][3]=0.054975871827661;
  gpt[10][3]=0.054975871827661;
  gpt[11][3]=0.054975871827661;
}

void Cfem::set_xyele(int element, double** xyele)
{
  int global_node;
  for(int node=0; node<pmesh->node; node++) {
    global_node = pmesh->EC[element][node];
    for(int dim=0; dim<pmesh->dimension; dim++)
      xyele[node][dim] = pmesh->XY[global_node][dim];
//    Diff_local[node][0] = DN[global_node][0];
  }
}

void Cfem::set_Phi_dPhiS_tet()
{
  int gauss;
  double s, t, z;
  for(gauss=0; gauss<pmesh->node; gauss++){
    s = gpt[gauss][ss];
    t = gpt[gauss][tt];
    z = gpt[gauss][zz];

    Phi[gauss][0] = 1.-s-t-z;
    Phi[gauss][1] = s;
    Phi[gauss][2] = t;
    Phi[gauss][3] = z;

    dPhiS[gauss][0][ss] = -1.;
    dPhiS[gauss][0][tt] = -1.;
    dPhiS[gauss][0][zz] = -1.;

    dPhiS[gauss][1][ss] = 1.;
    dPhiS[gauss][1][tt] = 0.;
    dPhiS[gauss][1][zz] = 0.;

    dPhiS[gauss][2][ss] = 0.;
    dPhiS[gauss][2][tt] = 1.;
    dPhiS[gauss][2][zz] = 0.;
 
    dPhiS[gauss][3][ss] = 0.;
    dPhiS[gauss][3][tt] = 0.;
    dPhiS[gauss][3][zz] = 1.;
    Transition(pmesh->node, pmesh->dimension, dPhiS[gauss], dPhiSt[gauss]);
  }
}


void Cfem::set_Phi_dPhiS_line()
{
  int gauss;
  double s;
  for(gauss=0; gauss<pmesh->node; gauss++){
    PetscPrintf(PETSC_COMM_WORLD, "line1\n");
    // gauss point at the node
    s = gpt[gauss][ss];
    // Two shape function values at the gauss point
    Phi[gauss][0] = (1-s)/2.;
    Phi[gauss][1] = (1+s)/2.;
    PetscPrintf(PETSC_COMM_WORLD, "line2\n");
    // Two shape function derivate values at the gauss point.
    dPhiS[gauss][0][ss] = -1./2.;
    dPhiS[gauss][1][ss] = 1./2.;
    Transition(pmesh->node, pmesh->dimension, dPhiS[gauss], dPhiSt[gauss]);
    PetscPrintf(PETSC_COMM_WORLD, "line3\n");
  }
  PetscPrintf(PETSC_COMM_WORLD, "line4\n");
}

// calculate Phi and PhiS values at every nodes in a local element.
void Cfem::set_Phi_dPhiS_square()
{
  int gauss;
  double s, t;
  for(gauss=0; gauss<pmesh->node; gauss++){
    // gauss point at the node
    s = gpt[gauss][ss];
    t = gpt[gauss][tt];
    // Four shape function values at the gauss point
    Phi[gauss][0] = (1.-s)*(1.-t)/4.;
    Phi[gauss][1] = (1.+s)*(1.-t)/4.;
    Phi[gauss][2] = (1.+s)*(1.+t)/4.;
    Phi[gauss][3] = (1.-s)*(1.+t)/4.;
    // Eight shape function derivate values at the gauss point.
    dPhiS[gauss][0][ss] = -(1.-t)/4.;
    dPhiS[gauss][0][tt] = -(1.-s)/4.;

    dPhiS[gauss][1][ss] = +(1.-t)/4.;
    dPhiS[gauss][1][tt] = -(1.+s)/4.;

    dPhiS[gauss][2][ss] = +(1.+t)/4.;
    dPhiS[gauss][2][tt] = +(1.+s)/4.;

    dPhiS[gauss][3][ss] = -(1.+t)/4.;
    dPhiS[gauss][3][tt] = +(1.-s)/4.;
    Transition(pmesh->node, pmesh->dimension, dPhiS[gauss], dPhiSt[gauss]);
  }
}

void Cfem::set_Phi_dPhiS_hex()
{
  int gauss;
  double s, t, z;
  for(gauss=0; gauss<pmesh->node; gauss++){
    s = gpt[gauss][ss];
    t = gpt[gauss][tt];
    z = gpt[gauss][zz];
    // Phi function
    Phi[gauss][0] = (1-s)*(1-t)*(1-z)/8.;
    Phi[gauss][1] = (1+s)*(1-t)*(1-z)/8.;
    Phi[gauss][2] = (1+s)*(1+t)*(1-z)/8.;
    Phi[gauss][3] = (1-s)*(1+t)*(1-z)/8.;
    Phi[gauss][4] = (1-s)*(1-t)*(1+z)/8.;
    Phi[gauss][5] = (1+s)*(1-t)*(1+z)/8.;
    Phi[gauss][6] = (1+s)*(1+t)*(1+z)/8.;
    Phi[gauss][7] = (1-s)*(1+t)*(1+z)/8.;
    // d shape function/ds, dt, dz
    dPhiS[gauss][0][ss] = -(1-t)*(1-z)/8.;
    dPhiS[gauss][0][tt] = -(1-s)*(1-z)/8.;
    dPhiS[gauss][0][zz] = -(1-s)*(1-t)/8.;

    dPhiS[gauss][1][ss] =  (1-t)*(1-z)/8.;
    dPhiS[gauss][1][tt] = -(1+s)*(1-z)/8.;
    dPhiS[gauss][1][zz] = -(1+s)*(1-t)/8.;

    dPhiS[gauss][2][ss] =  (1+t)*(1-z)/8.;
    dPhiS[gauss][2][tt] =  (1+s)*(1-z)/8.;
    dPhiS[gauss][2][zz] = -(1+s)*(1+t)/8.;

    dPhiS[gauss][3][ss] = -(1+t)*(1-z)/8.;
    dPhiS[gauss][3][tt] =  (1-s)*(1-z)/8.;
    dPhiS[gauss][3][zz] = -(1-s)*(1+t)/8.;

    dPhiS[gauss][4][ss] = -(1-t)*(1+z)/8.;
    dPhiS[gauss][4][tt] = -(1-s)*(1+z)/8.;
    dPhiS[gauss][4][zz] =  (1-s)*(1-t)/8.;
  
    dPhiS[gauss][5][ss] =  (1-t)*(1+z)/8.;
    dPhiS[gauss][5][tt] = -(1+s)*(1+z)/8.;
    dPhiS[gauss][5][zz] =  (1+s)*(1-t)/8.;

    dPhiS[gauss][6][ss] =  (1+t)*(1+z)/8.;
    dPhiS[gauss][6][tt] =  (1+s)*(1+z)/8.;
    dPhiS[gauss][6][zz] =  (1+s)*(1+t)/8.;

    dPhiS[gauss][7][ss] = -(1+t)*(1+z)/8.;
    dPhiS[gauss][7][tt] =  (1-s)*(1+z)/8.;
    dPhiS[gauss][7][zz] =  (1-s)*(1+t)/8.;
    Transition(pmesh->node, pmesh->dimension, dPhiS[gauss], dPhiSt[gauss]);
  }
}

void Cfem::set_Phi_dPhiS_prism()
{
  double a=0.445948490915965;
  double b=0.091576213509771;
  double r, s; 
  for(int ipoin=0;ipoin<12;ipoin++){
    r=gpt[ipoin][0];
    s=gpt[ipoin][1];
    //t=gpt[ipoin][2];
    //w=gpt[ipoin][3];

    Phi[ipoin][0]=(1.-r-s)*a;
    Phi[ipoin][1]=r*a;
    Phi[ipoin][2]=s*a;
    Phi[ipoin][3]=(1.-r-s)*b;
    Phi[ipoin][4]=r*b;
    Phi[ipoin][5]=s*b;

    dPhiS[ipoin][0][0]=-a;
    dPhiS[ipoin][0][1]=-a;
    dPhiS[ipoin][0][2]=-(1.-r-s)/2.;
    dPhiS[ipoin][1][0]=a;
    dPhiS[ipoin][1][1]=0.0;
    dPhiS[ipoin][1][2]=-0.5*r;
    dPhiS[ipoin][2][0]=0.0;
    dPhiS[ipoin][2][1]=a;
    dPhiS[ipoin][2][2]=-0.5*s;
    dPhiS[ipoin][3][0]=-b;
    dPhiS[ipoin][3][1]=-b;
    dPhiS[ipoin][3][2]=(1.-r-s)/2.;
    dPhiS[ipoin][4][0]= b;
    dPhiS[ipoin][4][1]=0.0;
    dPhiS[ipoin][4][2]=0.5*r;
    dPhiS[ipoin][5][0]=0.0;
    dPhiS[ipoin][5][1]=b;
    dPhiS[ipoin][5][2]=0.5*s;
    Transition(pmesh->node, pmesh->dimension, dPhiS[ipoin], dPhiSt[ipoin]);
  }
}

/*
void Cfem::set_Jacob(double** xyele, double*** dPhiS, double*** Jacob)
{
  int ix, jx, ia, gpt;
  for(gpt=0; gpt<pmesh->node; gpt++){
    for(ix=0; ix<pmesh->dimension; ix++){  // 2
      for(jx=0; jx<pmesh->dimension; jx++){  // 2
          Jacob[gpt][ix][jx] = 0.0;
      }
    }
  }
  for(gpt=0; gpt<pmesh->node; gpt++){

    for(ix=0; ix<pmesh->dimension; ix++){  // 2
      for(jx=0; jx<pmesh->dimension; jx++){  // 2
        for(ia=0; ia<pmesh->node; ia++){  // 4
          Jacob[gpt][ix][jx] += xyele[ia][jx]*dPhiS[gpt][ia][ix];
        }
      }
    }
  }
}
*/

double Cfem::det2D(double** mat)
{
  return mat[0][0]*mat[1][1] - mat[0][1]*mat[1][0];
}

double Cfem::det3D(double** mat)
{
  return  mat[0][0]*(mat[1][1]*mat[2][2]-mat[1][2]*mat[2][1])
         -mat[0][1]*(mat[1][0]*mat[2][2]-mat[1][2]*mat[2][0])
         +mat[0][2]*(mat[1][0]*mat[2][1]-mat[1][1]*mat[2][0]);
}


void Cfem::inverseMat(double** Ja, double** xJ, double* detJ)
{
  if(pmesh->dimension==3){
    for(int i=0; i<pmesh->node; i++){
      *detJ = det3D(Ja); 

      xJ[0][0]= (Ja[1][1]*Ja[2][2]-Ja[1][2]*Ja[2][1])/(*detJ);
      xJ[1][0]=-(Ja[1][0]*Ja[2][2]-Ja[1][2]*Ja[2][0])/(*detJ);
      xJ[2][0]= (Ja[1][0]*Ja[2][1]-Ja[1][1]*Ja[2][0])/(*detJ);

      xJ[0][1]=-(Ja[0][1]*Ja[2][2]-Ja[0][2]*Ja[2][1])/(*detJ);
      xJ[1][1]= (Ja[0][0]*Ja[2][2]-Ja[0][2]*Ja[2][0])/(*detJ);
      xJ[2][1]=-(Ja[0][0]*Ja[2][1]-Ja[0][1]*Ja[2][0])/(*detJ);

      xJ[0][2]= (Ja[0][1]*Ja[1][2]-Ja[0][2]*Ja[1][1])/(*detJ);
      xJ[1][2]=-(Ja[0][0]*Ja[1][2]-Ja[0][2]*Ja[1][0])/(*detJ);
      xJ[2][2]= (Ja[0][0]*Ja[1][1]-Ja[0][1]*Ja[1][0])/(*detJ);
    }
  }
  else if(pmesh->dimension==2){
    for(int i=0; i<pmesh->node; i++){
      *detJ = det2D(Ja);
      xJ[0][0]= Ja[1][1]/(*detJ);
      xJ[0][1]=-Ja[0][1]/(*detJ);
      xJ[1][0]=-Ja[1][0]/(*detJ);
      xJ[1][1]= Ja[0][0]/(*detJ);
    }
  }
  else{
    **xJ = 3;
  }
}
/*
// detJ, dPhiS[][], xJ[][] --> dPhiX[][]
void Cfem::set_dPhiX(double*** dPhiS, double*** xJ, double*** dPhiX)
{
  int ib, ia, i, j;
  for(ib=0; ib<pmesh->node; ib++){
    for(ia=0; ia<pmesh->node; ia++)
      for(i=0; i<pmesh->dimension; i++)
        dPhiX[ib][ia][i] = 0.0;
  
    for(ia=0; ia<pmesh->node; ia++)
      for(i=0; i<pmesh->dimension; i++)
        for(j=0; j<pmesh->dimension; j++) 
          dPhiX[ib][ia][i] += xJ[ib][i][j]*dPhiS[ib][ia][j];
  }
}*/

// row: row of amat, colrow: col of amat and row of bmat, col: col of bmat.
void Cfem::MultMatMat(int row, int colrow, int col, double** amat, double** bmat, double** cmat)
{
  for(int i=0; i<row; i++){
    for(int j=0; j<col; j++){
      cmat[i][j] = 0.;
      for(int k=0; k<colrow; k++){
        cmat[i][j] += amat[i][k]*bmat[k][j];
      }
    }
  }
}

void Cfem::Transition(int row, int col, double** amat, double** bmat)
{
  for(int i=0; i<row; i++)
    for(int j=0; j<col; j++)
      bmat[j][i] = amat[i][j];
}


// make local element matrix (make amat & bmat from dj[j], dPhiS[j])
void Cfem::set_amat_bmat(int element, double dt, double** amat, double** bmat)
{
  // Derivative of shape function in terms of physical coordinate
  int i, j, k, a, gauss;

  set_xyele(element, xyele);  // make xyele[node][pmesh->dimension] (local coordinate) from XY (global coordinate)
  for(gauss=0; gauss<pmesh->node; gauss++){ 
    MultMatMat(pmesh->dimension, pmesh->node, pmesh->dimension, dPhiSt[gauss], xyele, Jacob[gauss]); // make Jacob[node][pmesh->dimension][pmesh->dimension] from xyele and dPhiS
    /*for(i=0; i<pmesh->dimension; i++)
      for(j=0; j<pmesh->dimension; j++)
        PetscPrintf(PETSC_COMM_WORLD, "%f  ", Jacob[i][j]);
    PetscPrintf(PETSC_COMM_WORLD, "\n");*/
    inverseMat(Jacob[gauss], xJ[gauss], &detJ[gauss]); // xJ[node][dim][dim]: inverse Jacobian, detJ[node]: determinants of Jacobian
    Transition(pmesh->dimension, pmesh->dimension, xJ[gauss], xJt[gauss]);
    MultMatMat(pmesh->node, pmesh->dimension, pmesh->dimension, dPhiS[gauss], xJt[gauss], dPhiX[gauss]);
    Transition(pmesh->node, pmesh->dimension, dPhiX[gauss], dPhiXt[gauss]);

    MultMatMat(pmesh->node, pmesh->dimension, pmesh->node, dPhiX[gauss], dPhiXt[gauss], dPhiX2[gauss]);
  }

  for(i=0; i<pmesh->node; i++){
    for(j=0; j<pmesh->node; j++){
      amat[i][j]=0.0; bmat[i][j]=0.0;
    }
  }
  for(gauss=0; gauss<pmesh->node; gauss++){
    for(j=0; j<pmesh->node; j++) {
      for(k=0; k<pmesh->node; k++) {
        amat[j][k] += (Phi[gauss][j]*Phi[gauss][k]*detJ[gauss]/dt + dPhiX2[gauss][j][k]*detJ[gauss]*DE[element]); // conduction term
        bmat[j][k] += (Phi[gauss][j]*Phi[gauss][k]*detJ[gauss]);  // mass term matrix
      }
    }
  }
}

void Cfem::prepare_Amat_Bmat_bb_Voltsc(PetscScalar state0)
{
  int istart, iend, localsize;
  // Create sparse parallel matrix Amat & Bmat.
  MatCreateMPIAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,pmesh->all_nodes,pmesh->all_nodes,60,PETSC_NULL,60,PETSC_NULL,&Amat);
  MatCreateMPIAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,pmesh->all_nodes,pmesh->all_nodes,60,PETSC_NULL,60,PETSC_NULL,&Bmat);
  MatGetOwnershipRange(Amat,&istart,&iend); // range of rows in Amat for the process.
  localsize = iend-istart;
  printf("MatGetOwnershipRange(): istart: %d\t iend: %d\t localsize: %d\n", istart, iend, localsize);

  // Create and Set a size of Vector bb.
  VecCreate(PETSC_COMM_WORLD,&bb); // bb: vector
  VecSetSizes(bb, localsize, pmesh->all_nodes);// localsize: local size, all_node: global size
  VecSetFromOptions(bb); //

  // Duplicate bb to Volt[] & Volt1[] vector.
  for(int i=0; i<vtkout_states_vec.size(); i++)
    VecDuplicate(bb, &Volt[i]);   // create Volt vector as same format as bb
  for(int i=0; i<vtkout_algebraic_vec.size(); i++)
    VecDuplicate(bb, &Volt1[i]);   // create Volt vector as same format as bb

  // initialize vector Volt[0] with p_cell->p[0]
  VecSet(Volt[0], state0);

  // Scatter local vector Volt[] into global vector Voltsc[]
  for(int i=0; i< vtkout_states_vec.size(); i++){
    VecScatterCreateToAll(Volt[i], &ctx[i], &V_SEQ[i]);
    VecScatterBegin(ctx[i], Volt[i], V_SEQ[i], INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(ctx[i], Volt[i], V_SEQ[i], INSERT_VALUES, SCATTER_FORWARD);
    VecGetArray(V_SEQ[i], &Voltsc[i]);
  }
  for(int i=0; i< vtkout_algebraic_vec.size(); i++){
    VecScatterCreateToAll(Volt1[i], &ctx1[i], &V_SEQ1[i]);
    VecScatterBegin(ctx1[i], Volt1[i], V_SEQ1[i], INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(ctx1[i], Volt1[i], V_SEQ1[i], INSERT_VALUES, SCATTER_FORWARD);
    VecGetArray(V_SEQ1[i], &Voltsc1[i]);
  }
  PetscPrintf(PETSC_COMM_WORLD, "\n============vector preparation finished!\n");
//  if(mympi::rank==0) 
//    for(int i=0; i<pmesh->all_nodes; i++) printf("Voltsc[%d]=%lf\n", i, Voltsc[i]);
}

void Cfem::solve_Bmat(double dt)
{
  int its;
  VecSet(bb, 0);
  MatMult(Bmat,Volt[0],bb);  // bb = Bmat*Volt  --> bb

  VecScale(bb,1.0/dt);  // bb = bb*1/dt
  KSPSolve(ksp,bb,Volt[0]);  // ksp*Volt=bb --> Volt
  KSPGetIterationNumber(ksp,&its); // its: number of iteration done.

  VecScatterBegin(ctx[0],Volt[0],V_SEQ[0],INSERT_VALUES,SCATTER_FORWARD);
  VecScatterEnd(ctx[0],Volt[0],V_SEQ[0],INSERT_VALUES,SCATTER_FORWARD);
  VecGetArray(V_SEQ[0],&Voltsc[0]);
}


void Cfem::assemble_Amat_Bmat(double dt)
{
  int i, j, k;
  // derivative of shape function in terms of physical coordinate.
  // local coordinate
  // Set shape function (Phi and dPhiS) based on gauss points
  set_gpt(); PetscPrintf(PETSC_COMM_WORLD, "set_gpt done!\n");
  set_Phi_dPhiS(); PetscPrintf(PETSC_COMM_WORLD, "set_Phi_dPhiS done!\n");
  // ================= Stiffness matrix formulation ============================
  PetscPrintf(PETSC_COMM_WORLD, "amat bmat start!\n");
  for(int element=0; element<pmesh->all_elements; element++){
    if(element_proc[element]==mympi::rank){
      set_amat_bmat(element, dt, amat, bmat);
      // make global matrix
      for(j=0; j<pmesh->node; j++) {
        for(k=0; k<pmesh->node; k++) {
          MatSetValues(Amat, 1, &(pmesh->EC[element][j]), 1, &(pmesh->EC[element][k]), &(amat[j][k]), ADD_VALUES);
          MatSetValues(Bmat, 1, &(pmesh->EC[element][j]), 1, &(pmesh->EC[element][k]), &(bmat[j][k]), ADD_VALUES);
        }
      }
    }
  }
  PetscPrintf(PETSC_COMM_WORLD, "amat bmat done!\n");

  PetscPrintf(PETSC_COMM_WORLD, "Global stiffness matrix Amat, Bmat finished!\n");
  MatAssemblyBegin(Amat,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(Amat,MAT_FINAL_ASSEMBLY);

  PetscPrintf(PETSC_COMM_WORLD, "Amat finished!\n");
  MatAssemblyBegin(Bmat,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(Bmat,MAT_FINAL_ASSEMBLY);
  PetscPrintf(PETSC_COMM_WORLD, "Bmat finished!\n");
  KSPCreate(PETSC_COMM_WORLD,&ksp);
  KSPSetOperators(ksp, Amat, Amat, SAME_NONZERO_PATTERN);
//  KSPSetOperators(ksp, Amat, Amat);
  KSPSetType(ksp, KSPCG);  // linear matrix solving
  PC pc;
  KSPGetPC(ksp,&pc);
  PCSetType(pc,PCJACOBI);
  KSPSetFromOptions(ksp);
  //  KSPSetUp(ksp);
  PetscPrintf(PETSC_COMM_WORLD, "PETSc stiffness matrix Amat & Amat initialization finished\n");
}

void Cfem::assign_nodes_elements_proc()
{
// PETsc: Assign local_nodes[all_nodes/nprocs], node_proc[all_nodes], element_proc[all_elements]
  int proc_num = 0;  // 0~11
  int nodes_process = 0;   // 0~number of node per a process.
  for(int i=0; i<pmesh->all_nodes; i++){
    if(proc_num==mympi::rank)
      local_nodes[i] = nodes_process++;  // local_nodes[]: 0~number of node per a process
    node_proc[i] = proc_num++;   // 0~11
    if(proc_num>=mympi::size) proc_num=0;
  }

  proc_num = 0;
  int element_num = 0;
  for(int i=0; i<pmesh->all_elements; i++){
    if(proc_num==mympi::rank)
      element_num++;
    element_proc[i] = proc_num++;
    if(proc_num>=mympi::size) proc_num=0;
  }
}

void Cfem::assign_conductivity()
{
  int i, j;
 // Assign conductivity to all the nodes.

  for(i=0; i<pmesh->all_nodes; i++){
    isFibroticnode[i] = false;
    if(pmesh->celltype[i]==1) DN[i][0] = DC*pm.sf_diff; // tissue
    else if(pmesh->celltype[i]==4) DN[i][0] = DC*pm.sf_diff*2.0; // pect
    else DN[i][0] = DC*pm.sf_diff*0.5; // pect
  }

  if(strcmp(pm.variant,"Fibrosis")==0 ){   // Fibrosis setting
    std::ifstream ifs_fibrotic_nodes(pm.fibrotic_nodes);
    int node;
    int i = 0;
    while(ifs_fibrotic_nodes >> node){
      isFibroticnode[--node] = true;
      DN[node][0] *= pm.sf_diff_fibro;
      i++;
    }
    PetscPrintf(PETSC_COMM_WORLD, "Number of fibrotic nodes: %d\n", i);
    ifs_fibrotic_nodes.close();
  }


  for(i=0; i<pmesh->all_elements; i++){
    DE[i] = 0.;
    for(j=0; j<pmesh->node; j++)
      DE[i] += DN[pmesh->EC[i][j]][0]/(double)pmesh->node;  // average diffusion coefficient
  }
/*
  if(strcmp(pm.variant,"Fibrosis")==0 ){   // Fibrosis setting
    std::ifstream ifs_fibrotic_elements("fibrotic_nodes.csv");
    int iel;
    while(ifs_fibrotic_elements >> iel){
      isFibroticnode[pmesh->EC[iel][0]] = true;
      DE[iel] *= pm.sf_diff_fibro;
    }
    ifs_fibrotic_elements.close();
  }*/
}

void Cfem::set_gpt()
{
  if(!strcmp(pmesh->type, "tet")) set_gpt_tet();
  else if(!strcmp(pmesh->type, "hex")) set_gpt_hex();
  else if(!strcmp(pmesh->type, "fiber_real")) set_gpt_line();
  else if(!strcmp(pmesh->type, "square")) set_gpt_square();
  else if(!strcmp(pmesh->type, "prism")){
    for(int i=0; i<pmesh->node; i++)
      delete gpt[i];
    delete gpt;
    gpt = new double* [pmesh->node*2];
    for(int i=0; i<pmesh->node*2; i++)
      gpt[i] = new double [pmesh->dimension+1];
    set_gpt_prism();
  }
  printArray(gpt, pmesh->node, pmesh->dimension, "gpt");
}

void Cfem::set_Phi_dPhiS()
{
  if(!strcmp(pmesh->type, "tet")) set_Phi_dPhiS_tet();
  else if(!strcmp(pmesh->type, "hex")) set_Phi_dPhiS_hex();
  else if(!strcmp(pmesh->type, "fiber_real")) set_Phi_dPhiS_line();
  else if(!strcmp(pmesh->type, "square")) set_Phi_dPhiS_square();
  else if(!strcmp(pmesh->type, "prism")) set_Phi_dPhiS_prism();
  else exit(1);
  printArray(Phi, pmesh->node, pmesh->dimension, "Phi");
  printArray(dPhiS, pmesh->node, pmesh->node, pmesh->dimension, "dPhiS");
}



void Cfem::openVTK(FILE* fp)
{
  char tmp[200];
  for(int i=0; i<pmesh->all_nodes+5; i++)
    fgets(tmp, 200, fp);
  for(int i=0; i<pmesh_surf->all_elements*2+2+3; i++)
    fgets(tmp, 200, fp);
  for(int i=0; i<pmesh->all_nodes; i++)
    fscanf(fp, "%lf", &pmesh->celltype[i]);
}

void Cfem::printArray(double*** array, int isize, int jsize, int ksize, char* name)
{
  PetscPrintf(PETSC_COMM_WORLD, "%s[%d][%d][%d]\n", name, isize, jsize, ksize);
  for(int i=0; i<isize; i++){
    for(int j=0; j<jsize; j++){
      for(int k=0; k<ksize; k++)
        PetscPrintf(PETSC_COMM_WORLD, "%lf\t", array[i][j][k]);
      PetscPrintf(PETSC_COMM_WORLD, "\n");
    }
    PetscPrintf(PETSC_COMM_WORLD, "\n");
  }
}

void Cfem::printArray(double** array, int isize, int jsize, char* name)
{
  PetscPrintf(PETSC_COMM_WORLD, "%s[%d][%d]\n", name, isize, jsize);
  for(int i=0; i<isize; i++){
    for(int j=0; j<jsize; j++){
      PetscPrintf(PETSC_COMM_WORLD, "%lf\t", array[i][j]);
    }
    PetscPrintf(PETSC_COMM_WORLD, "\n");
  }
}
