#ifndef PARAM_H
#define PARAM_H

typedef struct parameters
{
  double t_max;
  double t_start;
  double t_write_vtk;   // write VTK files  (twrite)
  bool is_using_output;              // (iopt)
  int is_s1s2;              // (iopt)
  bool is_apply_mutation;
  bool is_using_edison;
  double sf_diff;       // (Dscale)
  double sf_diff_fibro;    // (Dscalefib);
  double sf_diff_fiber;   // (Dscale_fiber);
  char femesh_dir[200];
  char bemesh_dir[200];
  char fibrotic_nodes[100]; // (fibrosis)
  char heart_mesh[100]; // (heart)
  char surface_mesh[100];
  char torso_mesh[100]; // (torso)
  char pace1_nodes[100];   // (pace)
  char pace2_nodes[100];  // (pace2)
  char ionmodel[255];
  char variant[20];
//  char ic[100];
  double celltype;  // cell types
  double dt;
  double dt_write;
  double bcl_init;    // (ibclb)
  double bcl_decrement;   // ibcld
  double bcl_end;   // ibcle
  double bcl_pace2;   // s2bcl
  int num_pace1;    // irepeat
  char stim_type[20]; // stimtype
  double stim_dur;    // stimdur
  double stim_amp;  // stim_amp
//  double stim_v;
//  double stim_i;

//  double vtarget;
  double v_apd90;  // membrane potential value at 90%
  double v_apd50;  // membrane potential value at 50%
//  double apdslop;
  char mesh_type[20];
  char output_mesh_type[20]; //outputmesh
  float apd_fiber;
  float erp_fiber;
  bool is_lbbb;   // lbbb
  bool is_rbbb; // rbbb;
  bool is_crt;   // crt
  int t_crt; // crttime;
  bool is_ecg; // ecg
  bool is_bspm;
  int ecgnode[10];
  char openVTK[200];
  char fibermesh_dir[300];
  char fiber[300];
  bool is_bidomain;  // ifbidomain

  char vtkout_states[255];
  char vtkout_algebraic[255];
  int simulation_mode; 
  int scaled_gate;
  double mutation_scaling;

  double gcal_scale;
  double gk1_scale;
  double gkr_scale;
  double gks_scale;
  double gna_scale;
  double gnal_scale;
  /*CiPA params*/
  char hill_file[1024];
  char herg_file[1024];
  char drug_name[100];
  char concs[100];
  double cmax;
  int is_print_vm;
  int dimension;
} param_t;


#endif
