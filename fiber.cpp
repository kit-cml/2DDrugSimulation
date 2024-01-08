#include "fiber.h"
// TN2006 setting
#define DEPOL 30  // mV
#define REPOL -83  // mV

// ORd setting
//#define DEPOL 60  // mV
//#define REPOL -86  // mV

Cfiber::Cfiber(Cmesh* pmeshfiber, param_t parameter): Cfem(pmeshfiber)
{
  pm = parameter;
  setMeshConditions();
}

Cfiber::Cfiber(Cmesh* pmeshfiber, Cbspm* pbspm, param_t parameter) : Cfem(pmeshfiber)
{
  pm = parameter;
  this->pbspm = pbspm;
  setMeshConditions();
}

bool Cfiber::setMeshConditions()
{
  //FILE* fp;
  //char fname[200];
  time=0.;
  ISTIMNODE = new int [pmesh->all_nodes];
  ifERP = new int [pmesh->all_nodes];
  Volt_t = new double [pmesh->all_nodes];
  perwave = new double [pmesh->all_nodes];
  ifastele = new int [pmesh->all_nodes];
  Volt = new double [pmesh->all_nodes];
  for(int i=0; i<pmesh->all_nodes; i++){
    Volt[i]= REPOL;  // 모든노드에 휴지기 전압으로 초기화
    Volt_t[i] = 0.0;  // 모든노드에 타이머 시간을 0으로 초기화
    ifERP[i] = 0;  // 모든노드에 타이머 시간을 0으로 초기화
  }
  // 모든 엘리먼트에 전도유무 0으로 초기화
  for(int i=0; i<pmesh->all_elements; i++) perwave[i]=0.0; 

  PetscPrintf(PETSC_COMM_WORLD, "BEFORE PACINGSITE\n");
  pmesh->AssignPacingSites(pm.fibermesh_dir);
  //pmesh->AssignPacingSites(".");
  PetscPrintf(PETSC_COMM_WORLD, "AFTER PACINGSITE\n");
  PetscPrintf(PETSC_COMM_WORLD, "BEFORE TERMINAL\n");
  pmesh->FindTerminalNodes();
  PetscPrintf(PETSC_COMM_WORLD, "AFTER TERMINAL\n");
  if(pm.is_lbbb || pm.is_rbbb) pmesh->AssignBBBnode(pm.fibermesh_dir, pm.is_rbbb);

  fpvm = fopen("vmcheck_fiber.plt", "wt");
  PetscPrintf(PETSC_COMM_WORLD, "pmesh->all_terminals = %d\n", pmesh->all_terminals);
  return true;
}

void Cfiber::AVpacing()
{
  //=========AV node automaticity==========
 // PetscPrintf(PETSC_COMM_WORLD, "nfoci[0]: %d all_fibers: %d\n", nfoci[0], pmesh->all_fibers);
  for(int i=0; i<pmesh->all_fibers; i++){
    if(int(Volt[pmesh->nfoci[i]])==REPOL){
      Volt[pmesh->nfoci[i]] = DEPOL;
      ifERP[pmesh->nfoci[i]] = 1;
      Volt_t[pmesh->nfoci[i]] = 0.0;
    }
  }
//  for(int i=0; i<all_elements; i++)
//    perwave[i] = 0;
}
               
int Cfiber::processing(int argc, char **args, param_t parameter)
{
  pm = parameter;
  double time = 0.;
  int iter=0;
  int iter_bcl= int(double(pm.bcl_init)/pm.dt);
  int iter_dtw = int(pm.dt_write/pm.dt);


  PetscPrintf(PETSC_COMM_WORLD, "start processing\n");
  while(time< (pm.num_pace1 * pm.bcl_init)/pm.dt ){
    fprintf(fpvm, "%lf %lf %lf %lf\n", time, Volt[0], Volt[1], Volt[2]);
    if((time>=pm.t_write_vtk) && (iter%iter_dtw==0)) writePF(time); 

    if(!(iter%iter_bcl)) AVpacing();
    PropagationPF();

    time += pm.dt;
    iter++;
  }
  PetscPrintf(PETSC_COMM_WORLD, "kmlim finished\n");
  fclose(fpvm);
  PetscPrintf(PETSC_COMM_WORLD, "kmlim finished\n");
  return 1;
}

int Cfiber::PropagationPF()
{
  int iel, i, i1, i2, iw, ierp;
  double dels, percent;
  speeds=pm.sf_diff_fiber;
  double speed000 = speeds;  // 빠른파이버이면 3배
  for(iel=0; iel<pmesh->all_elements; iel++){ // 퍼킨지 모든 엘리먼트에 대하여
    i1 = pmesh->EC[iel][0]; // node 1
    i2 = pmesh->EC[iel][1]; // node 2
    if(pmesh->nbbb==i1 || pmesh->nbbb==i2) continue;  // for Bundle Branch Block
    iw = int(Volt[i1] +Volt[i2]); // 양쪽 노드의 흥분성을 더함
    ierp = ifERP[i1] + ifERP[i2];
    if(iw==REPOL+DEPOL && ierp==1){ // ! only 1 node has been activated : 한쪽노드만 흥분했을때.
      //speed000=speeds/3.0;  // 느린파이버이면 3배느리게
      dels = sqrt(pow((pmesh->XY[i1][x]-pmesh->XY[i2][x]),2)+pow((pmesh->XY[i1][y]-pmesh->XY[i2][y]), 2)+pow((pmesh->XY[i1][z]-pmesh->XY[i2][z]),2)); // 엘리먼트길이
      percent = speed000*pm.dt/dels; // 엘리먼트에서 흥분파의 진행 퍼센트
      if(percent>=0.3) { // 한번의 타임스텝에 엘리먼트의 30%만큼 전도가 일어나면
        PetscPrintf(PETSC_COMM_WORLD, "conduction percent in element is more than 0.3  --> percent=%lf\n", percent);
        PetscPrintf(PETSC_COMM_WORLD, "fiber.cpp : decrease Dscale_fiber or dt\n");
        exit(0);
      }
      perwave[iel] += percent; // 엘리먼트에서 흥분파의 진행 퍼센트
      if(perwave[iel]>=0.99){ // 99퍼센트가 넘으면
        perwave[iel]=0.;
        if(int(Volt[i1])==DEPOL && !ifERP[i2]){ // 왼쪽노드가 말단노드가아니면서 흥분했을 경우
          Volt[i2] = DEPOL;  // 오른쪽노드의 AP 를 30mV 로 세팅함
          ifERP[i2] = 1;  // 오른쪽노드의 AP 를 30mV 로 세팅함
          Volt_t[i2] = 0.0; // 오른쪽느드 타이머 작동  
          //PetscPrintf( PETSC_COMM_WORLD, "ERP Left\n" );
        }
        else if(int(Volt[i2])==DEPOL && !ifERP[i1]){  // 오른쪽 노드가 말단노드가 아니면서 흥분했을 경우
          Volt[i1] = DEPOL; // 왼쪽노드를 흥분시킴
          ifERP[i1] = 1;  // 
          Volt_t[i1]=0.0; // 왼쪽노드 타이머 작동
          //PetscPrintf( PETSC_COMM_WORLD, "ERP Right\n" );
        }
        else{
          printf("Cfiber::PropagationPF(): Error----Next node is still refractory period\n\n");
        }
      }
    }
  }

  for(i=0; i<pmesh->all_nodes; i++){
    if(Volt_t[i]>pm.apd_fiber)  Volt[i] = REPOL; //  50ms APD for purkinje
    if(Volt_t[i]>pm.erp_fiber) ifERP[i] = 0;  // refractory period
    if(ifERP[i]) Volt_t[i] += pm.dt;
    //Volt_t[i] += pm.dt;
  }
  return 1;
}


void Cfiber::writePF(double time)
{
  FILE* fppf;
  char fname[100];
  int i;

  if(!strcmp(pm.mesh_type, "fiber")) 
    sprintf(fname,"V_%d.vtk",(int)floor(time+0.5));
  else 
    sprintf(fname,"PF_%d.vtk",(int)floor(time+0.5));

//  sprintf(fname,"PF_%f.vtk",time);
  fppf=fopen(fname,"wt");

  fprintf(fppf, "# vtk DataFile Version 3.0\n");
  fprintf(fppf, "vtk output\n");
  fprintf(fppf, "ASCII\n");
  fprintf(fppf, "DATASET UNSTRUCTURED_GRID\n");
  fprintf(fppf, "POINTS %d float\n", pmesh->all_nodes);
  for(i=0; i<pmesh->all_nodes; i++){
    fprintf(fppf, "%lf %lf %lf\n", pmesh->XY[i][x], pmesh->XY[i][y], pmesh->XY[i][z]); // node2
  }

  fprintf(fppf, "CELLS %d %d\n", pmesh->all_elements, pmesh->all_elements*3);
  for(i=0; i<pmesh->all_elements; i++){
    fprintf(fppf, "2 %d %d\n", pmesh->EC[i][0], pmesh->EC[i][1]);
  }

  fprintf(fppf, "CELL_TYPES %d\n", pmesh->all_elements);
  for(i=0; i<pmesh->all_elements; i++)
    fprintf(fppf, "3\n");  // 3: fiber element
    
  fprintf(fppf, "POINT_DATA %d\n", pmesh->all_nodes);
  fprintf(fppf, "SCALARS Vm float\n");
  fprintf(fppf, "LOOKUP_TABLE default\n");
  for(i=0; i<pmesh->all_nodes; i++)
    fprintf(fppf, "%lf\n", Volt[i]);  // Vm
  
  fprintf(fppf, "SCALARS idoub float\n");
  fprintf(fppf, "LOOKUP_TABLE default\n");
  for(i=0; i<pmesh->all_nodes; i++){
    fprintf(fppf, "%d\n", pmesh->idoub[i]);  // Vm
  }

  fclose(fppf);
}

void Cfiber::writeVTK()
{}

void Cfiber::writeoutput()
{}


void Cfiber::openVTK(FILE* fp)
{
   char tmp[200];
  for(int i=0; i<pmesh->all_nodes+5; i++)
    fgets(tmp, 200, fp);

  for(int i=0; i<pmesh->all_elements*2+2+3; i++)
    fgets(tmp, 200, fp);

  for(int i=0; i<pmesh->all_nodes; i++){
    fscanf(fp, "%d", &pmesh->celltype[i]);
  } 
}

void Cfiber::computeBSPM()
{
  for(int i=0; i<pbspm->pmesh->all_nodes; i++) // torso nodes
    pbspm->pmesh->celltype[i] = Radiate(pmesh->XY[i]);
}

void Cfiber::computeECG()
{
  for(int i=0; i<9; i++)
    pbspm->ecgl.vnode[i] = Radiate(pbspm->pmesh->XY[pm.ecgnode[i]-1]);
}

double Cfiber::Radiate(double* elecrode)
{
  //int i, j;
  return 1.0;
}

