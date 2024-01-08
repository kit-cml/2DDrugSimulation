#include "mesh.h"
#include "stdio.h"


Cmesh::Cmesh(char* filename)
{
  loadInp(filename);
}

void Cmesh::AssignBBBnode(char* dir, int ifrbbb)
{
  FILE* fp;
  int itmp;
  char fname[200];
  sprintf(fname, "%s/bbb.dat", dir);
  if(!(fp = fopen(fname, "rt"))) {PetscPrintf(PETSC_COMM_WORLD, "Cfiber::setMeshConditions(): failed to open %s\n", fname); exit(-1);}
  fscanf(fp, "%d\n", &itmp);
  fscanf(fp, "%d\n", &nbbb);
  if(ifrbbb) {
    fscanf(fp, "%d\n", &nbbb);
    PetscPrintf(PETSC_COMM_WORLD, "rbbb: %d\n", nbbb);
  }
  else PetscPrintf(PETSC_COMM_WORLD, "lbbb: %d\n", nbbb);
  nbbb--;
  fclose(fp);
}

void Cmesh::AssignPacingSites(char* dir)
{
  FILE* fp;
  char fname[200];
  sprintf(fname, "%s/nfoci.dat", dir);
  if(!(fp = fopen(fname, "rt"))) {PetscPrintf(PETSC_COMM_WORLD, "Cfiber::setMeshConditions(): failed to open %s\n", fname); exit(-1);}
  fscanf(fp, "%d\n", &all_fibers);
  PetscPrintf(PETSC_COMM_WORLD, "ALL FIBER: %d\n", all_fibers);
  nfoci = new int [all_fibers];
  for(int i=0; i<all_fibers; i++){
    fscanf(fp, "%d\n", &nfoci[i]);
    nfoci[i]--;
  }
  PetscPrintf(PETSC_COMM_WORLD, "nfoci[0]: %d\n", nfoci[0]);
  fclose(fp);
}

void Cmesh::FindTerminalNodes()
{
  iPFnode = new int [all_nodes];
  idoub = new int [all_nodes];
  for(int i=0; i<all_nodes; i++) idoub[i]=0;  
  // 모든요소에 기여하는 모든노드의 연결성확인, 1이면 말단노드 2이면 말단노>드아님.
  for(int i=0; i<all_elements; i++){
    idoub[EC[i][0]]++; 
    idoub[EC[i][1]]++;
  }
  all_terminals=0; // number of peripheral nodes: 말단노드의 개수
  for(int i=0; i<all_nodes; i++)  // 모든노드에 대해서 말단노드 확인
    if(idoub[i]==1) iPFnode[all_terminals++]=i;     // iPFnode: peripheral node number
}
void Cmesh::loadInp(char* filename)
{
  FILE* fp;
  int i;
  char fname[1000], tmp[1000];
  char* token;
  nbbb = -1;

  //============Open heart mesh========================
  sprintf(fname, "%s", filename);
  if(!(fp = fopen(fname, "rt"))) {PetscPrintf(PETSC_COMM_WORLD, "Cmesh::Cmesh(): failed to open %s\n", fname); exit(-1);}
  fgets(tmp, 1000, fp);
  while(strncmp(tmp, "CML", 3)) fgets(tmp, 1000, fp);  // read until 'CML'
  fscanf(fp, "%d %d %d %d %d\n", &all_nodes, &all_elements, &target[0], &target[1], &target[2]);
  PetscPrintf(PETSC_COMM_WORLD, "%s--> all_nodes: %d, all_elements: %d\n", filename, all_nodes, all_elements);
  target[0]--;   target[1]--;   target[2]--;
  while(strncmp(tmp, "*NODE\n", 5)) fgets(tmp, 1000, fp);// read until '*NODE'


  XY = new double* [all_nodes];
  for(i=0; i<all_nodes; i++)
    XY[i] = new double[3]; // to read x,y,z coordinates.
  //============Read Data (coordiate, element connectivity)========================
  for(i=0; i<all_nodes; i++) fgets(tmp, 1000, fp); // skip all_nodes
  fgets(tmp, 1000, fp);  // read '**HWCOLOR..'
  fgets(tmp, 1000, fp);  // read '*ELEMENT..'i

  token = strtok(tmp, ",");
  token = strtok(NULL, ",");
  if(!strncmp(token, "TYPE=S2", 7)){
    strcpy(type, "fiber");
    node = 2;
    dimension = 1;
  }
  else if(!strncmp(token, "TYPE=S4", 7)){
    strcpy(type, "square");
    node = 4;  
    dimension = 2;
  }
  else if(!strncmp(token, "TYPE=S3", 7)){
    strcpy(type, "tri");
    node = 3;  
    dimension = 3;
  }
  else if(!strncmp(token, "TYPE=C3D4", 9)){
    strcpy(type, "tet");
    node = 4; 
    dimension = 3;
  }
  else if(!strncmp(token, "TYPE=C3D8I", 10)){
    strcpy(type, "hex");
    node = 8;
    dimension = 3;
  }
  else if(!strncmp(token, "TYPE=T3D2", 9)){
    strcpy(type, "fiber");
    node = 2;
    dimension = 3;
  }
  PetscPrintf(PETSC_COMM_WORLD, "Cmesh::Cmesh(): mesh.type: %s\n", type);
  fseek(fp,0 ,SEEK_SET); // move file pointer to beginning
  PetscPrintf(PETSC_COMM_WORLD, "BEFORE NODE %d\n", all_nodes);
  while(strncmp(tmp, "*NODE\n", 5)) {
    fgets(tmp, 1000, fp); // read until '*NODE'
  }
  for(i=0; i<all_nodes; i++){
    fgets(tmp, 1000, fp);
    token = strtok(tmp, ","); // read node number
    for(int j=0; j<3; j++){  // read x,y,z coordinates
      token = strtok(NULL, ",");
      XY[i][j] = atof(token);
    }
  }

  PetscPrintf(PETSC_COMM_WORLD, "FINISHED NODE\n");
  PetscPrintf(PETSC_COMM_WORLD, "BEFORE EC\n");

  EC = new int* [all_elements];
  for(i=0; i<all_elements; i++)
    EC[i] = new int [node];

  PetscPrintf(PETSC_COMM_WORLD, "FINISHED EC\n");

  fgets(tmp, 1000, fp); // read **HWCOLOR COMP ..
  fgets(tmp, 1000, fp); // read *ELEMENT ..
  for(i=0; i<all_elements; i++){
    fgets(tmp, 1000, fp);
    token = strtok(tmp, ","); // read element number
    for(int j=0; j<node; j++){
      token = strtok(NULL, ",");  
      EC[i][j] = atoi(token)-1;
    }

  }
  PetscPrintf(PETSC_COMM_WORLD, "FINISHED HWCOLOR\n");

  fgets(tmp, 1000, fp); // read 'CELLTYPE'
  if(strncmp(tmp, "CELLTYPE", 8) && strncmp(tmp, "eof", 3)) {
    PetscPrintf(PETSC_COMM_WORLD, "Error reading file: %s: no CELLTYPE or no eof!\n",  filename);
    exit(-1);
  }
  PetscPrintf(PETSC_COMM_WORLD, "FINISHED CELLTYPE\n");

  if(!strncmp(tmp, "CELLTYPE", 8)){
    celltype = new int [all_nodes];
    for(i=0; i<all_nodes; i++){
      fgets(tmp, 1000, fp);
      celltype[i] = atoi(tmp);
    }
    fgets(tmp, 1000, fp); // read eof
    if(!strncmp(tmp, "eof", 3)) 
      PetscPrintf(PETSC_COMM_WORLD, "CELLTYPE included: %s -->\n",  filename);
    else {
      PetscPrintf(PETSC_COMM_WORLD, "Error reading file:%s --> no eof!\n", filename);
      exit(-1);
    }
  }

  else if(!strncmp(tmp, "eof", 3))
    PetscPrintf(PETSC_COMM_WORLD, "No CELLTYPE included: %s -->\n", filename);

  else {
    PetscPrintf(PETSC_COMM_WORLD, "Error reading file: %s --> no CELLTYPE!\n", filename);
    exit(-1);
  }
  fclose(fp);
  PetscPrintf(PETSC_COMM_WORLD, "FINISHED MESHING\n");
}

