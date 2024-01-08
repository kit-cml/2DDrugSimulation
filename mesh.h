#ifndef MESH_H
#define MESH_H

#include <petscksp.h>
#include <petscsys.h>
#include <petscpc.h>
#include <petscmat.h>
#include <petsctime.h>
#include "ccalcexception.h"
#include <stdio.h>
#include <fstream>


const int x = 0;
const int y = 1;
const int z = 2;
const int ss = 0;
const int tt = 1;
const int zz = 2;

typedef struct ECGleads_tag
{
  int inode[20];
  double vnode[20];
  double ref;
  double xL[3];
  double av[3];
  double precordial[6];
} ECGLEADS;

class Cmesh
{
public:
  int dimension;
  char type[20];   // mesh type
  int all_nodes;   // number of all nodes
  int all_elements;  // number all elements
  int all_fibers;
  int node;  // number of nodes in one element;
  int target[5];
  double** XY;  // [all_nodes] coordinates of all nodes
  int** EC;  // [all_elements][node] element node connectivity
  int* celltype;  // [all_nodes] -> endo, M, epi
  int* idoub;  // nodes which belong to two elements
  int all_terminals; // terminal node of purkinje
  int nbbb;  // node number of bundle branch block in purkinje
  int* nfoci;  // number of foci
  int* iPFnode; //[NUM_NODE];
  int* inearsurf;  // [all_nodes]
  public:
  Cmesh(char* filename);
  void loadInp(char* filename);
  void FindTerminalNodes();
  void AssignBBBnode(char* dir, int ifrbbb);
  void AssignPacingSites(char* dir);
};

#endif
