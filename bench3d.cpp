#include <stdio.h>
#include <algorithm> // for std::find
#include <iterator>  // for std::begin, std::end
#include <string>
#include <ctime>
#include <cstring>
// POSIX library
#include <libgen.h>

#include "heart.h"
#include "fiber.h"
#include "mesh.h"
#include "bspm.h"
#include "DrugSimulation/modules/glob_type.hpp"

//#include "paramset.cpp"
#include <omp.h>

#define IS_EDISON

bool use_edison_format = false;

void default_params(param_t *pm);
void assign_params(int argc, char *args[], param_t *pm);
void edison_assign_params(char *fn, param_t *pm);
void edison_read_arguments(int argc, char *args[], param_t *pm);
int exists(const char *fname);

// int main(int argc, char **args)
int main(int argc, char *args[])
{

  clock_t start;
  double duration;

  // Pestsc initialization
  PetscInitialize(&argc, &args, PETSC_NULL, PETSC_NULL);
  MPI_Comm_rank(PETSC_COMM_WORLD, &mympi::rank);
  MPI_Comm_size(PETSC_COMM_WORLD, &mympi::size);
  printf("bench3d: mympi::rank: %d mympi::size: %d\n", mympi::rank, mympi::size);

  param_t pm;
  default_params(&pm);

  // ============== specific case for Edison ===========================
#ifdef IS_EDISON
    use_edison_format = true;

    PetscPrintf(PETSC_COMM_WORLD, "EDISON ENTER\n"  );
    edison_read_arguments(argc, args, &pm);
//#else
//  assign_params(argc, args, &pm);
#endif
  // ============== end of specific case for Edison ====================


  PetscPrintf(PETSC_COMM_WORLD, "Meshtype: %s\n", pm.mesh_type);
  start = std::clock();
  if(!strcmp(pm.mesh_type, "fiber")){
    Cmesh* pmeshfiber = new Cmesh(pm.fiber);
    Cfiber* pfiber = new Cfiber(pmeshfiber, pm);
    pfiber->processing(argc, args, pm);
  }

  else if(!strcmp(pm.mesh_type, "tet")||!strcmp(pm.mesh_type, "square")||!strcmp(pm.mesh_type, "hex")){
    PetscPrintf(PETSC_COMM_WORLD, "ENTERING TET/SQUARE!!\n");
    Cmesh* pmeshheart = new Cmesh(pm.heart_mesh);
    pmeshheart->dimension = pm.dimension;
    Cfem* pheart = new Cheart(pmeshheart);
    pheart->processing(argc, args, pm);
  }
  else if(!strcmp(pm.mesh_type, "tetfiber")){
    Cmesh* pmeshheart = new Cmesh(pm.heart_mesh);
    Cmesh* pmeshfiber = new Cmesh(pm.fiber);
    Cfiber* pfiber = new Cfiber(pmeshfiber, pm);
    Cfem* pheart = new Cheart(pmeshheart, pfiber);
    pheart->processing(argc, args, pm);
  }
  else if(!strcmp(pm.mesh_type, "tetbspm")){
    Cmesh* pmeshheart = new Cmesh(pm.heart_mesh);
    Cmesh* pmeshtorso = new Cmesh(pm.torso_mesh);
    Cbspm* pbspm = new Cbspm(pmeshtorso);
    Cfem* pheart = new Cheart(pmeshheart, pbspm);
    pheart->processing(argc, args, pm);
  }
  else if(!strcmp(pm.mesh_type, "fiberbspm")){
    Cmesh* pmeshfiber = new Cmesh(pm.fiber);
    Cmesh* pmeshtorso = new Cmesh(pm.torso_mesh);
    Cbspm* pbspm = new Cbspm(pmeshtorso);
    Cfiber* pfiber = new Cfiber(pmeshfiber, pbspm, pm);
    pfiber->processing(argc, args, pm);
  }
  else if(!strcmp(pm.mesh_type, "tetfiberbspm")){
    Cmesh* pmeshheart = new Cmesh(pm.heart_mesh);
    Cmesh* pmeshfiber = new Cmesh(pm.fiber);
    Cmesh* pmeshtorso = new Cmesh(pm.torso_mesh);
    Cbspm* pbspm = new Cbspm(pmeshtorso);
    Cfiber* pfiber = new Cfiber(pmeshfiber, pm);
    Cfem* pheart = new Cheart(pmeshheart, pfiber, pbspm);
    pheart->processing(argc, args, pm);
  }
  else printf("meshtype selection wrong!\n");
#ifdef OMPBSPM
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  // at the end of the simulation, we move all result files to result folder
  // ensure only 1 process do this
#ifdef MESH_TOO_BIG_SO_NO_ZIPPING
  if (mympi::rank == 0) {
    // only do this in Edison
    if (use_edison_format) {
      system("rm -rf result");
      system("mkdir -p result/vtk");

      // add header to vmcheck.plt file
      system("sed -i '1 i\\Time stim_node target[0] target[1] target[2]' vmcheck.plt");
      system("./convert-to-edison-format.sh vmcheck.plt ms > vmcheck_edison.plt");
      //system("zip result_electrical.zip *.vtk *.scalar");
      // move all result files
      //system("mv *.plt result");
      //system("mv *.vtk result/vtk");
      //system("mv *.scalar result/vtk");
      //system("mv *.zip result");
      //system("mv output.* result");
    }
  }
#endif
  duration = ((std::clock() - start) / (double)CLOCKS_PER_SEC) / 60.;
  PetscPrintf(PETSC_COMM_WORLD, "simulation finished perfectly for %lf minutes\n", duration);

  PetscFinalize();

  return EXIT_SUCCESS;
}

void default_params(param_t *pm)
{
  // General options
  pm->t_max = 10000;      // s1-s2 protocol or not
  pm->is_s1s2 = 0;      // s1-s2 protocol or not
  pm->t_write_vtk       = 0;      // time to write vtk files
  pm->dt           = 0.1; // delta time
  pm->dt_write          = 10.0; // delta time to write vtk files.
  pm->stim_dur      = 1.0; // stimulation duration in apex, purkije terminal, or fiber at first node.
  pm->stim_amp      = 5.0; // stimulation amplitude in apex, purkije terminal, or fiber at first node.
  strcpy(pm->mesh_type, "tet"); // stimulation type (it can be voltage or current)
  strcpy(pm->stim_type, "voltage"); // stimulation type (it can be voltage or current)
  pm->is_using_output         = 0; // if you simulate from beginning iopt=0, if you use IC with output.dat iopt=1.

  // Heart single cell simulation
  strcpy(pm->ionmodel, "tn2006epi");  // cell type
  strcpy(pm->variant, "ori"); // mutation type
  pm->num_pace1      = 3; // Number of S1 stimulus
  pm->bcl_init        = 600; // BCL at beginning
  pm->bcl_decrement        = 10;  // BCL decresement (when you make APDR curve)
  pm->bcl_end        = 88000010;  // BCL at the end (when you make APDR curve) Set this higher than ibclb when you don't attend to make APDR curve
  pm->bcl_pace2        = 700; // BCL of S2 stimulation.

  // Heart 2D simulation
  pm->sf_diff       = 1.5; // scale factor of conductivity.
  strcpy(pm->fibrotic_nodes, "fibrosis/5.dat"); // random node numbers of fibrosis
  pm->sf_diff_fibro    = 0.7; // scale factor of concuctivity under fibrosis.

  // Heart 3D simulation
  strcpy(pm->femesh_dir, "/scratch4/maincode/mesh_elect/vent/human_hf_214319"); // mesh directory
  // strcpy(pm->femesh_dir, "./hhf_214319"); // mesh directory
  strcpy(pm->bemesh_dir, "/scratch4/maincode/mesh_elect/torso"); // mesh directory
  strcpy(pm->heart_mesh, "heart_het.inp"); // mesh type of heart
  strcpy(pm->pace1_nodes, "heart_apex.dat"); // node numbers of apex pacing
  strcpy(pm->pace2_nodes, "heart_L.dat"); // node numbers of LV freewall for S2 stimulus.
  strcpy(pm->output_mesh_type, "surface"); // output mesh type (it can be solid or surface)
  pm->is_lbbb         = 0; // LBBB option
  pm->is_rbbb         = 0; // RBBB option
  pm->is_crt          = 0; // CRt opttion
  pm->t_crt      = 2; // CRT timing (ms)

  // Heart 3D with BSPM simulation
  pm->is_ecg         = 0; // ECG generation option during heart simulation
  pm->is_bspm        = 0; // ECG generation option during heart simulation
  strcpy(pm->torso_mesh, "p+torso48_closed.inp"); // mesh for torso
  pm->ecgnode[0]   = 472; // ECG nodes
  pm->ecgnode[1]   = 393;
  pm->ecgnode[2]   = 239;
  pm->ecgnode[3]   = 28;
  pm->ecgnode[4]   = 29;
  pm->ecgnode[5]   = 30;
  pm->ecgnode[6]   = 480;
  pm->ecgnode[7]   = 254;
  pm->ecgnode[8]   = 252;

  // Nerve fiber or Purkinje simulation
  pm->sf_diff_fiber = 0.4; // conduction velocity of neuron
  pm->apd_fiber    = 300; // APD of neuron
  pm->erp_fiber    = 350; // ERP of neuron

  strcpy( pm->vtkout_algebraic, "" );
  strcpy( pm->vtkout_states, "V" );

  // Conductance scaleiplier
  pm->gcal_scale = 1.;
  pm->gk1_scale = 1.;
  pm->gkr_scale = 1.;
  pm->gks_scale = 1.;
  pm->gna_scale = 1.;
  pm->gnal_scale = 1.;

  // BSPM or EEG simulation.
  strcpy(pm->openVTK, "/scratch4/kmlim/ecg4"); // Results of heart or fibers (for BSPM or EEG)
  strcpy(pm->fiber, "line.inp"); 
  //strcpy(pm->fibermesh_dir, "/scratch4/maincode/mesh_elect/vent/human_hf_214319/purkinje"); // Results of heart or fibers (for BSPM or EEG)
  pm->is_bidomain    = 0; // monodomain vs. bidomain

  // AP Slope default
  pm->ap_slope = 1.1;
  pm->celltype = 0.;

}

void assign_params(int argc, char *args[], param_t *pm)
{
  if(argc>1){
    for(int count = 1; count < argc; count += 2){
      if(!strcmp(args[count], "-twrite"))    pm->t_write_vtk       = atof(args[count + 1]);
      //else if(!strcmp(args[count], "-dimension"))  strcpy(pm->dimension, args[count + 1]);
      else if(!strcmp(args[count], "-ionmodel"))  strcpy(pm->ionmodel, args[count + 1]);
      else if(!strcmp(args[count], "-variant"))   strcpy(pm->variant, args[count + 1]);
      else if(!strcmp(args[count], "-tmax"))      pm->t_max           = atof(args[count + 1]);
      else if(!strcmp(args[count], "-dt"))        pm->dt           = atof(args[count + 1]);
      else if(!strcmp(args[count], "-irepeat"))   pm->num_pace1      = atoi(args[count + 1]);
      else if(!strcmp(args[count], "-is_s1s2"))   pm->is_s1s2      = atoi(args[count + 1]);
      else if(!strcmp(args[count], "-ibclb"))     pm->bcl_init        = atof(args[count + 1]);
      else if(!strcmp(args[count], "-ibcld"))     pm->bcl_decrement        = atoi(args[count + 1]);
      else if(!strcmp(args[count], "-ibcle"))     pm->bcl_end        = atoi(args[count + 1]);
      else if(!strcmp(args[count], "-s2bcl"))     pm->bcl_pace2        = atoi(args[count + 1]);
      else if(!strcmp(args[count], "-stimdur"))   pm->stim_dur      = atof(args[count + 1]);
      else if(!strcmp(args[count], "-stimamp"))   pm->stim_amp      = atof(args[count + 1]);
      else if(!strcmp(args[count], "-meshtype"))  strcpy(pm->mesh_type, args[count + 1]);
      else if(!strcmp(args[count], "-stimtype"))  strcpy(pm->stim_type, args[count + 1]);
      else if(!strcmp(args[count], "-iopt"))      pm->is_using_output         = atoi(args[count + 1]);
      else if(!strcmp(args[count], "-Dscale"))    pm->sf_diff       = atof(args[count + 1]);
      else if(!strcmp(args[count], "-Dscalefib")) pm->sf_diff_fibro    = atof(args[count + 1]);
      else if(!strcmp(args[count], "-CV_fiber"))      pm->sf_diff_fiber = atof(args[count + 1]);
      else if(!strcmp(args[count], "-fibrotic_nodes")) strcpy(pm->fibrotic_nodes, args[count + 1]);
      else if(!strcmp(args[count], "-femesh_dir"))    strcpy(pm->femesh_dir, args[count + 1]);
      else if(!strcmp(args[count], "-bemesh_dir"))    strcpy(pm->bemesh_dir, args[count + 1]);
      else if(!strcmp(args[count], "-heart"))         strcpy(pm->heart_mesh, args[count + 1]);
      else if(!strcmp(args[count], "-pace1"))         strcpy(pm->pace1_nodes, args[count + 1]);
      else if(!strcmp(args[count], "-pace2"))         strcpy(pm->pace2_nodes, args[count + 1]);
      else if(!strcmp(args[count], "-dtw"))           pm->dt_write          = atof(args[count + 1]);
      else if(!strcmp(args[count], "-outputmesh"))    strcpy(pm->output_mesh_type, args[count + 1]);
      else if(!strcmp(args[count], "-apd_fiber"))     pm->apd_fiber    = atof(args[count + 1]);
      else if(!strcmp(args[count], "-erp_fiber"))     pm->erp_fiber    = atof(args[count + 1]);
      else if(!strcmp(args[count], "-lbbb"))          pm->is_lbbb         = atoi(args[count + 1]);
      else if(!strcmp(args[count], "-rbbb"))          pm->is_rbbb         = atoi(args[count + 1]);
      else if(!strcmp(args[count], "-crt"))           pm->is_crt          = atoi(args[count + 1]);
      else if(!strcmp(args[count], "-crttime"))       pm->t_crt      = atoi(args[count + 1]);
      else if(!strcmp(args[count], "-torso"))         strcpy(pm->torso_mesh, args[count + 1]);
      else if(!strcmp(args[count], "-iecg"))          pm->is_ecg         = atoi(args[count + 1]);
      else if(!strcmp(args[count], "-ibspm"))         pm->is_bspm         = atoi(args[count + 1]);
      else if(!strcmp(args[count], "-ecgnode[0]"))    pm->ecgnode[0]   = atoi(args[count + 1]);
      else if(!strcmp(args[count], "-ecgnode[1]"))    pm->ecgnode[1]   = atoi(args[count + 1]);
      else if(!strcmp(args[count], "-ecgnode[2]"))    pm->ecgnode[2]   = atoi(args[count + 1]);
      else if(!strcmp(args[count], "-ecgnode[3]"))    pm->ecgnode[3]   = atoi(args[count + 1]);
      else if(!strcmp(args[count], "-ecgnode[4]"))    pm->ecgnode[4]   = atoi(args[count + 1]);
      else if(!strcmp(args[count], "-ecgnode[5]"))    pm->ecgnode[5]   = atoi(args[count + 1]);
      else if(!strcmp(args[count], "-ecgnode[6]"))    pm->ecgnode[6]   = atoi(args[count + 1]);
      else if(!strcmp(args[count], "-ecgnode[7]"))    pm->ecgnode[7]   = atoi(args[count + 1]);
      else if(!strcmp(args[count], "-ecgnode[8]"))    pm->ecgnode[8]   = atoi(args[count + 1]);
      else if(!strcmp(args[count], "-openVTK"))       strcpy(pm->openVTK, args[count + 1]);
      else if(!strcmp(args[count], "-fibermesh_dir")) strcpy(pm->fibermesh_dir, args[count + 1]);
      else if(!strcmp(args[count], "-fiber"))         strcpy(pm->fiber, args[count + 1]);
      else if(!strcmp(args[count], "-ifbidomain"))    pm->is_bidomain   = atoi(args[count + 1]);
      else if(!strcmp(args[count], "-vtkout_algebraic"))  strcpy( pm->vtkout_algebraic ,args[count + 1] );
      else if(!strcmp(args[count], "-vtkout_states"))  strcpy( pm->vtkout_states ,args[count + 1] );


      //else                                            printf("Unknown command option: %s\n", args[count]);
      //PetscPrintf(PETSC_COMM_WORLD, "%s: %s\n", args[count], args[count+1]);
    }
  }
  else{
    printf("Invalid arguments number.\n");
    exit(1);
  }
}

/**
 * @brief Read Edison-related input files and process parameter from input deck.
 * Please note this function will modify several simulation parameters which
 * declared as global variables.
 * 
 * @param_t file path to input deck file
 */
void edison_read_arguments(int argc, char *args[], param_t *pm) {

  char filename[1024] = "";
  char outputzip[1024] = "";
  char heartMesh[1024] = "";
  char pace1Mesh[1024] = "";
  char pace2Mesh[1024] = "";
  char fiberMesh[1024] = "";
  char torsoMesh[1024] = "";
  char ecSurf[1024] = "";
  char fibrotic_nodes[1024] = "";
  char hill_file[1024] = "";
  char herg_file[1024] = "";


  for (int i = 1; i < argc; i += 2) {
    if (!strcmp(args[i], "-input_deck"))
      strcpy(filename, args[i + 1]);
    else if (!strcmp(args[i], "-outputzip"))
      strcpy(outputzip, args[i + 1]);
    else if (!strcmp(args[i], "-heartmesh"))
      strcpy(heartMesh, args[i + 1]);
    else if (!strcmp(args[i], "-pace1mesh"))
      strcpy(pace1Mesh, args[i + 1]);
    else if (!strcmp(args[i], "-pace2mesh"))
      strcpy(pace2Mesh, args[i + 1]);
    else if (!strcmp(args[i], "-fibermesh"))
      strcpy(fiberMesh, args[i + 1]);
    else if(!strcmp(args[i], "-torsomesh"))
      strcpy(torsoMesh, args[i + 1]);
    else if (!strcmp(args[i], "-ecsurf"))
      strcpy(ecSurf, args[i + 1]);
    else if(!strcmp(args[i], "-fibrotic_nodes"))
      strcpy(fibrotic_nodes, args[i + 1]);
    else if(!strcmp(args[i], "-hill_file"))
      strcpy(hill_file, args[i + 1]);
    else if(!strcmp(args[i], "-herg_file"))
      strcpy(herg_file, args[i + 1]);

  }

  // TODO: param_t meshtype tet



  // if user uploaded a custom mesh
  // check wheter he also upload pace1 and pace2 files
  if (!strcmp(pace1Mesh, "") || !strcmp(pace2Mesh, "")) {
    PetscPrintf(PETSC_COMM_WORLD, "If you are using custom (uploaded) mesh, please provide pace1 and pace2 too.\n");
    exit(1);
  }

  // get heart file name
  strcpy(pm->heart_mesh, heartMesh);
  // get heart file name
  strcpy(pm->surface_mesh, ecSurf);
  // get pace1 file name
  strcpy(pm->pace1_nodes, pace1Mesh);
  // get pace2 file name
  strcpy(pm->pace2_nodes, pace2Mesh);

  // get fiber file name
  string temp = fiberMesh;
  char fiberdir[255];
  strcpy(fiberdir, temp.substr(0, strrchr(fiberMesh,'/')-fiberMesh+1).c_str());
  strcpy(pm->fibermesh_dir, fiberdir);
  strcpy(pm->fiber, fiberMesh);

  PetscPrintf(PETSC_COMM_WORLD, "WTF %s\n", pm->fiber);
  strcpy(pm->fibrotic_nodes, fibrotic_nodes);
  strcpy(pm->hill_file, hill_file);
  strcpy(pm->herg_file, herg_file);

  PetscPrintf(PETSC_COMM_WORLD, "Input Deck Filename %s\n", filename);
  PetscPrintf(PETSC_COMM_WORLD, "heart %s\nsurface %s\n pace1 %s\npace2 %s\nfiber %s\n", pm->heart_mesh, pm->surface_mesh, pm->pace1_nodes, pm->pace2_nodes, pm->fiber);
  PetscPrintf(PETSC_COMM_WORLD, "Fibrotic filename: %s\n", pm->fibrotic_nodes );
  PetscPrintf(PETSC_COMM_WORLD, "CiPA Hill sample: %s\n", pm->hill_file );
  PetscPrintf(PETSC_COMM_WORLD, "CiPA hERG sample: %s\n", pm->herg_file );


  edison_assign_params(filename, pm);

  // if user use iopt == 1
  if (pm->is_using_output == 1) {
    ifstream zip(outputzip);
    // if user uploaded the required output.*
    if (zip.good()) {
      char unzipCommand[1024];
      strcpy(unzipCommand, "unzip -o ");
      strcat(unzipCommand, outputzip);
      strcat(unzipCommand, "\n");
      // extract the compressed file in current working directory
      system(unzipCommand);
    } else {
      PetscPrintf(PETSC_COMM_WORLD, "Please upload output.* files as zip file when you set iopt to 1.\n");
      exit(2);
    }
  }
}

void edison_assign_params(char *fn, param_t *pm)
{
  FILE* fp_inputdeck = fopen(fn, "r");
  char buf_char[1024];
  while(1)
  {
    fscanf(fp_inputdeck,"%s", buf_char);
    if(feof(fp_inputdeck)) break;

    if (!strcmp(buf_char, "tmax")) { fscanf(fp_inputdeck, "%*s %lf", &pm->t_max); PetscSynchronizedPrintf(PETSC_COMM_WORLD, "tmax: %lf\n", pm->t_max);}
    else if (!strcmp(buf_char, "tstart")) { fscanf(fp_inputdeck, "%*s %lf", &pm->t_start); PetscSynchronizedPrintf(PETSC_COMM_WORLD, "tstart: %lf\n", pm->t_start);}
    else if (!strcmp(buf_char, "twrite")) { fscanf(fp_inputdeck, "%*s %lf", &pm->t_write_vtk); PetscSynchronizedPrintf(PETSC_COMM_WORLD, "twrite: %lf\n", pm->t_write_vtk);}
    else if (!strcmp(buf_char, "variant")) { fscanf(fp_inputdeck, "%*s %s", &pm->variant); PetscSynchronizedPrintf(PETSC_COMM_WORLD, "variant: %s\n", pm->variant);}
    else if (!strcmp(buf_char, "meshtype")) { fscanf(fp_inputdeck, "%*s %s", &pm->mesh_type); PetscSynchronizedPrintf(PETSC_COMM_WORLD, "meshtype: %s\n", pm->mesh_type);}
    else if (!strcmp(buf_char, "outputmesh")) { fscanf(fp_inputdeck, "%*s %s", &pm->output_mesh_type); PetscSynchronizedPrintf(PETSC_COMM_WORLD, "outputmesh: %s\n", pm->output_mesh_type);}
    else if (!strcmp(buf_char, "dt")) { fscanf(fp_inputdeck, "%*s %lf", &pm->dt); PetscSynchronizedPrintf(PETSC_COMM_WORLD, "dt: %lf\n", pm->dt);}
    else if (!strcmp(buf_char, "irepeat")) { fscanf(fp_inputdeck, "%*s %d", &pm->num_pace1); PetscSynchronizedPrintf(PETSC_COMM_WORLD, "irepeat: %d\n", pm->num_pace1);}
    else if (!strcmp(buf_char, "simulation_mode")) { fscanf(fp_inputdeck, "%*s %d", &pm->simulation_mode); PetscSynchronizedPrintf(PETSC_COMM_WORLD, "simulation_mode: %d\n", pm->simulation_mode);}
    else if (!strcmp(buf_char, "is_s1s2")) { fscanf(fp_inputdeck, "%*s %d", &pm->is_s1s2); PetscSynchronizedPrintf(PETSC_COMM_WORLD, "is_s1s2: %d\n", pm->is_s1s2);}
    else if (!strcmp(buf_char, "ibclb")) { fscanf(fp_inputdeck, "%*s %lf", &pm->bcl_init); PetscSynchronizedPrintf(PETSC_COMM_WORLD, "ibclb: %lf\n", pm->bcl_init);}
    else if (!strcmp(buf_char, "s2bcl")) { fscanf(fp_inputdeck, "%*s %lf", &pm->bcl_pace2); PetscSynchronizedPrintf(PETSC_COMM_WORLD, "s2bcl: %lf\n", pm->bcl_pace2);}
    else if (!strcmp(buf_char, "stimdur")) { fscanf(fp_inputdeck, "%*s %lf", &pm->stim_dur); PetscSynchronizedPrintf(PETSC_COMM_WORLD, "stimdur: %lf\n", pm->stim_dur);}
    else if (!strcmp(buf_char, "stimamp")) { fscanf(fp_inputdeck, "%*s %lf", &pm->stim_amp); PetscSynchronizedPrintf(PETSC_COMM_WORLD, "stimamp: %lf\n", pm->stim_amp);}
    else if (!strcmp(buf_char, "iopt")) { fscanf(fp_inputdeck, "%*s %d", &pm->is_using_output); PetscSynchronizedPrintf(PETSC_COMM_WORLD, "iopt: %d\n", pm->is_using_output);}
    else if (!strcmp(buf_char, "Dscale")) { fscanf(fp_inputdeck, "%*s %lf", &pm->sf_diff); PetscSynchronizedPrintf(PETSC_COMM_WORLD, "Dscale: %lf\n", pm->sf_diff);}
    else if (!strcmp(buf_char, "Dscalefib")) { fscanf(fp_inputdeck, "%*s %lf", &pm->sf_diff_fiber);PetscSynchronizedPrintf(PETSC_COMM_WORLD, "Dscalefib: %lf\n", pm->sf_diff_fiber);}
    else if (!strcmp(buf_char, "dtw")) { fscanf(fp_inputdeck, "%*s %lf", &pm->dt_write);PetscSynchronizedPrintf(PETSC_COMM_WORLD, "dtw: %lf\n", pm->dt_write);}
    else if (!strcmp(buf_char, "vtkout_states")) { fscanf(fp_inputdeck, "%*s %s", &pm->vtkout_states); PetscSynchronizedPrintf(PETSC_COMM_WORLD, "vtkout_states: %s\n", pm->vtkout_states);}
    else if (!strcmp(buf_char, "dimension")) { fscanf(fp_inputdeck, "%*s %d", &pm->dimension); PetscSynchronizedPrintf(PETSC_COMM_WORLD, "dimension: %d\n", pm->dimension);}
    else if (!strcmp(buf_char, "simulation_mode")) { fscanf(fp_inputdeck, "%*s %d", &pm->simulation_mode); PetscSynchronizedPrintf(PETSC_COMM_WORLD, "simulation_mode: %d\n", pm->simulation_mode);}
    else if (!strcmp(buf_char, "drug_name")) { fscanf(fp_inputdeck, "%*s %s", &pm->drug_name); PetscSynchronizedPrintf(PETSC_COMM_WORLD, "drug_name: %s\n", pm->drug_name);}
    else if (!strcmp(buf_char, "concs")) { fscanf(fp_inputdeck, "%*s %s", &pm->concs); PetscSynchronizedPrintf(PETSC_COMM_WORLD, "concs: %s\n", pm->concs);}
    //else if (!strcmp(buf_char, "A1656D_Type")) { fscanf(fp_inputdeck, "%*s %d", &glob_var::A1656D_mode); PetscSynchronizedPrintf(PETSC_COMM_WORLD, "A1656D_mode: %d\n", glob_var::A1656D_mode);}

 
    else if (!strcmp(buf_char, "gks_scale")) { fscanf(fp_inputdeck, "%*s %lf", &pm->gks_scale); PetscSynchronizedPrintf(PETSC_COMM_WORLD, "gks_scale: %lf\n", pm->gks_scale);}
    else if (!strcmp(buf_char, "gkr_scale")) { fscanf(fp_inputdeck, "%*s %lf", &pm->gkr_scale); PetscSynchronizedPrintf(PETSC_COMM_WORLD, "gkr_scale: %lf\n", pm->gkr_scale);}
    else if (!strcmp(buf_char, "gk1_scale")) { fscanf(fp_inputdeck, "%*s %lf", &pm->gk1_scale); PetscSynchronizedPrintf(PETSC_COMM_WORLD, "gk1_scale: %lf\n", pm->gk1_scale);}
    else if (!strcmp(buf_char, "gna_scale")) { fscanf(fp_inputdeck, "%*s %lf", &pm->gna_scale); PetscSynchronizedPrintf(PETSC_COMM_WORLD, "gna_scale: %lf\n", pm->gna_scale);}
    else if (!strcmp(buf_char, "gbna_scale")) { fscanf(fp_inputdeck, "%*s %lf", &pm->gbna_scale); PetscSynchronizedPrintf(PETSC_COMM_WORLD, "gbna_scale: %lf\n", pm->gbna_scale);}
    else if (!strcmp(buf_char, "gcal_scale")) { fscanf(fp_inputdeck, "%*s %lf", &pm->gcal_scale); PetscSynchronizedPrintf(PETSC_COMM_WORLD, "gcal_scale: %lf\n", pm->gcal_scale);}
    else if (!strcmp(buf_char, "gnal_scale")) { fscanf(fp_inputdeck, "%*s %lf", &pm->gnal_scale); PetscSynchronizedPrintf(PETSC_COMM_WORLD, "gnal_scale: %lf\n", pm->gnal_scale);}
    else if (!strcmp(buf_char, "gpca_scale")) { fscanf(fp_inputdeck, "%*s %lf", &pm->gpca_scale); PetscSynchronizedPrintf(PETSC_COMM_WORLD, "gpca_scale: %lf\n", pm->gpca_scale);}
    else if (!strcmp(buf_char, "gpk_scale")) { fscanf(fp_inputdeck, "%*s %lf", &pm->gpk_scale); PetscSynchronizedPrintf(PETSC_COMM_WORLD, "gpk_scale: %lf\n", pm->gpk_scale);}
    else if (!strcmp(buf_char, "gbca_scale")) { fscanf(fp_inputdeck, "%*s %lf", &pm->gbca_scale); PetscSynchronizedPrintf(PETSC_COMM_WORLD, "gbca_scale: %lf\n", pm->gbca_scale);}
    else if (!strcmp(buf_char, "gto_scale")) { fscanf(fp_inputdeck, "%*s %lf", &pm->gto_scale); PetscSynchronizedPrintf(PETSC_COMM_WORLD, "gto_scale: %lf\n", pm->gto_scale);}
    else if (!strcmp(buf_char, "ap_slope")) { fscanf(fp_inputdeck, "%*s %lf", &pm->ap_slope); PetscSynchronizedPrintf(PETSC_COMM_WORLD, "ap_slope: %lf\n", pm->ap_slope);}
    else if (!strcmp(buf_char, "celltype")) { fscanf(fp_inputdeck, "%*s %lf", &pm->celltype); PetscSynchronizedPrintf(PETSC_COMM_WORLD, "celltype: %lf\n", pm->celltype);}


    else if (!strcmp(buf_char, "is_ecg")) { fscanf(fp_inputdeck, "%*s %d", &pm->is_ecg); PetscSynchronizedPrintf(PETSC_COMM_WORLD, "is_ecg: %d\n", pm->is_ecg);}
    else if (!strcmp(buf_char, "is_bspm")) { fscanf(fp_inputdeck, "%*s %d", &pm->is_bspm); PetscSynchronizedPrintf(PETSC_COMM_WORLD, "is_bspm: %d\n", pm->is_bspm);}
    else if (!strcmp(buf_char, "ecgnode[0]")) { fscanf(fp_inputdeck, "%*s %d", &pm->ecgnode[0]); PetscSynchronizedPrintf(PETSC_COMM_WORLD, "ecgnode[0]: %d\n", pm->ecgnode[0]);}
    else if (!strcmp(buf_char, "ecgnode[1]")) { fscanf(fp_inputdeck, "%*s %d", &pm->ecgnode[1]); PetscSynchronizedPrintf(PETSC_COMM_WORLD, "ecgnode[1]: %d\n", pm->ecgnode[1]);}
    else if (!strcmp(buf_char, "ecgnode[2]")) { fscanf(fp_inputdeck, "%*s %d", &pm->ecgnode[2]); PetscSynchronizedPrintf(PETSC_COMM_WORLD, "ecgnode[2]: %d\n", pm->ecgnode[2]);}
    else if (!strcmp(buf_char, "ecgnode[3]")) { fscanf(fp_inputdeck, "%*s %d", &pm->ecgnode[3]); PetscSynchronizedPrintf(PETSC_COMM_WORLD, "ecgnode[3]: %d\n", pm->ecgnode[3]);}
    else if (!strcmp(buf_char, "ecgnode[4]")) { fscanf(fp_inputdeck, "%*s %d", &pm->ecgnode[4]); PetscSynchronizedPrintf(PETSC_COMM_WORLD, "ecgnode[4]: %d\n", pm->ecgnode[4]);}
    else if (!strcmp(buf_char, "ecgnode[5]")) { fscanf(fp_inputdeck, "%*s %d", &pm->ecgnode[5]); PetscSynchronizedPrintf(PETSC_COMM_WORLD, "ecgnode[5]: %d\n", pm->ecgnode[5]);}
    else if (!strcmp(buf_char, "ecgnode[6]")) { fscanf(fp_inputdeck, "%*s %d", &pm->ecgnode[6]); PetscSynchronizedPrintf(PETSC_COMM_WORLD, "ecgnode[6]: %d\n", pm->ecgnode[6]);}
    else if (!strcmp(buf_char, "ecgnode[7]")) { fscanf(fp_inputdeck, "%*s %d", &pm->ecgnode[7]); PetscSynchronizedPrintf(PETSC_COMM_WORLD, "ecgnode[7]: %d\n", pm->ecgnode[7]);}
    else if (!strcmp(buf_char, "ecgnode[8]")) { fscanf(fp_inputdeck, "%*s %d", &pm->ecgnode[8]); PetscSynchronizedPrintf(PETSC_COMM_WORLD, "ecgnode[8]: %d\n", pm->ecgnode[8]);}

    else {
      printf("Error Invalid value name :: %s\n", buf_char);
    }
  }

  fclose(fp_inputdeck);
}





/**
 * @brief Check existance of a file.
 *
 * @param_t fname filename
 * @return int 0 if not exist, 1 otherwise.
 */
int exists(const char *fname)
{
  FILE *file;
  if ((file = fopen(fname, "r"))) {
    fclose(file);
    return 1;
  }
  return 0;
}
