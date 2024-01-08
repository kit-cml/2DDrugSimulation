#ifndef PARAMSET_CPP
#define PARAMSET_CPP

#include "CEP2Cell/cellmodels/cellml.hpp"
#include <cstdio>
#include <string.h>
#include <stdlib.h>

void paramset(char* fn, param_t* pm)
{
  FILE* fp_inputdeck = fopen(fn, "r");
  char buf_char[1024];
  while(1)
  {
    fscanf(fp_inputdeck,"%s", buf_char);
    printf("done reading %s\n", buf_char);
    if(feof(fp_inputdeck)) break;

    //if (!strcmp(buf_char, "tmax")) { fscanf(fp_inputdeck, "%*s %lf", &pm->tmax); }
    if (!strcmp(buf_char, "twrite")) { fscanf(fp_inputdeck, "%*s %lf", &pm->t_write_vtk); }
    //else if (!strcmp(buf_char, "mesh")) { fscanf(fp_inputdeck, "%*s %s", &pm->mesh); }
    else if (!strcmp(buf_char, "ionmodel")) { fscanf(fp_inputdeck, "%*s %s", &pm->ionmodel); }
    else if (!strcmp(buf_char, "variant")) { fscanf(fp_inputdeck, "%*s %s", &pm->variant); }
    else if (!strcmp(buf_char, "dt")) { fscanf(fp_inputdeck, "%*s %lf", &pm->dt); }
    else if (!strcmp(buf_char, "dtw")) { fscanf(fp_inputdeck, "%*s %lf", &pm->dt_write); }
    else if (!strcmp(buf_char, "irepeat")) { fscanf(fp_inputdeck, "%*s %d", &pm->num_pace1); }
    else if (!strcmp(buf_char, "ibclb")) { fscanf(fp_inputdeck, "%*s %lf", &pm->bcl_init); }
    else if (!strcmp(buf_char, "ibcld")) { fscanf(fp_inputdeck, "%*s %d", &pm->bcl_decrement); }
    else if (!strcmp(buf_char, "ibcle")) { fscanf(fp_inputdeck, "%*s %d", &pm->bcl_end); }
    else if (!strcmp(buf_char, "sbcl")) { fscanf(fp_inputdeck, "%*s %d", &pm->bcl_pace2); }
    else if (!strcmp(buf_char, "stimamp")) { fscanf(fp_inputdeck, "%*s %lf", &pm->stim_amp); }
    else if (!strcmp(buf_char, "stimdur")) { fscanf(fp_inputdeck, "%*s %lf", &pm->stim_dur); }
    //else if (!strcmp(buf_char, "vth")) { fscanf(fp_inputdeck, "%*s %lf", &pm->vth); }
    //else if (!strcmp(buf_char, "vtarget")) { fscanf(fp_inputdeck, "%*s %lf", &pm->vtarget); }
    else if (!strcmp(buf_char, "apd90")) { fscanf(fp_inputdeck, "%*s %lf", &pm->v_apd90); }
    //else if (!strcmp(buf_char, "apdslop")) { fscanf(fp_inputdeck, "%*s %lf", &pm->apdslop); }
    else if (!strcmp(buf_char, "stimtype")) { fscanf(fp_inputdeck, "%*s %s", &pm->stim_type); }
    //else if (!strcmp(buf_char, "ifischem")) { fscanf(fp_inputdeck, "%*s %d", &pm->ifischem); }
    //else if (!strcmp(buf_char, "ifbrugada")) { fscanf(fp_inputdeck, "%*s %d", &pm->ifbrugada); }
    else if (!strcmp(buf_char, "iopt")) { fscanf(fp_inputdeck, "%*s %d", &pm->is_using_output); }
    else if (!strcmp(buf_char, "Dscale")) { fscanf(fp_inputdeck, "%*s %lf", &pm->sf_diff); }
    else if (!strcmp(buf_char, "Dscalefib")) { fscanf(fp_inputdeck, "%*s %lf", &pm->sf_diff_fibro); }
    //else if (!strcmp(buf_char, "CV_nerve")) { fscanf(fp_inputdeck, "%*s %lf", &pm->Dscale_nerve); }
    else if (!strcmp(buf_char, "fibrosis")) { fscanf(fp_inputdeck, "%*s %s", &pm->fibrotic_nodes); }
    else if (!strcmp(buf_char, "outputmesh")) { fscanf(fp_inputdeck, "%*s %s", &pm->output_mesh_type); }
    else if (!strcmp(buf_char, "femesh_dir")) { fscanf(fp_inputdeck, "%*s %s", &pm->femesh_dir); }
    else if (!strcmp(buf_char, "heart")) { fscanf(fp_inputdeck, "%*s %s", &pm->heart_mesh); }
    else if (!strcmp(buf_char, "pace1")) { fscanf(fp_inputdeck, "%*s %s", &pm->pace1_nodes); }
    else if (!strcmp(buf_char, "pace2")) { fscanf(fp_inputdeck, "%*s %s", &pm->pace2_nodes); }
    // TODO change this!
    // else if (!strcmp(buf_char, "lbbb")) { fscanf(fp_inputdeck, "%*s %d%*s", &pm->is_lbbb); }
    // else if (!strcmp(buf_char, "rbbb")) { fscanf(fp_inputdeck, "%*s %d%*s", &pm->is_rbbb); }
    else if (!strcmp(buf_char, "lbbb")) { fscanf(fp_inputdeck, "%*s %d", &pm->is_lbbb); }
    else if (!strcmp(buf_char, "rbbb")) { fscanf(fp_inputdeck, "%*s %d", &pm->is_rbbb); }
    else if (!strcmp(buf_char, "crt")) { fscanf(fp_inputdeck, "%*s %d", &pm->is_crt); }
    else if (!strcmp(buf_char, "crttime")) { fscanf(fp_inputdeck, "%*s %d", &pm->t_crt); }
    else if (!strcmp(buf_char, "torso")) { fscanf(fp_inputdeck, "%*s %s", &pm->torso_mesh); }
    else if (!strcmp(buf_char, "iecg")) { fscanf(fp_inputdeck, "%*s %d", &pm->is_ecg); }
    else if (!strcmp(buf_char, "ecgnode0")) { fscanf(fp_inputdeck, "%*s %d", &pm->ecgnode[0]); }
    else if (!strcmp(buf_char, "ecgnode1")) { fscanf(fp_inputdeck, "%*s %d", &pm->ecgnode[1]); }
    else if (!strcmp(buf_char, "ecgnode2")) { fscanf(fp_inputdeck, "%*s %d", &pm->ecgnode[2]); }
    else if (!strcmp(buf_char, "ecgnode3")) { fscanf(fp_inputdeck, "%*s %d", &pm->ecgnode[3]); }
    else if (!strcmp(buf_char, "ecgnode4")) { fscanf(fp_inputdeck, "%*s %d", &pm->ecgnode[4]); }
    else if (!strcmp(buf_char, "ecgnode5")) { fscanf(fp_inputdeck, "%*s %d", &pm->ecgnode[5]); }
    else if (!strcmp(buf_char, "ecgnode6")) { fscanf(fp_inputdeck, "%*s %d", &pm->ecgnode[6]); }
    else if (!strcmp(buf_char, "ecgnode7")) { fscanf(fp_inputdeck, "%*s %d", &pm->ecgnode[7]); }
    else if (!strcmp(buf_char, "ecgnode8")) { fscanf(fp_inputdeck, "%*s %d", &pm->ecgnode[8]); }
    else {
      printf("Error Invalid value name :: %s\n", buf_char);
      exit(1);
    }
  }

  /*
  // empty ic
  strcpy(pm->ic, "");
  strcat(pm->ic, pm->ionmodel);
  strcat(pm->ic, "/init_");
  strcat(pm->ic, pm->variant);
  strcat(pm->ic, ".dat");

  pm->vtarget = -45.0;
  pm->ifischem = 0;
  pm->ifbrugada = -45.0;
  // strcpy(pm->fibrotic_nodes, "fibrosis/5.dat");
  strcpy(pm->femesh_dir, "human_hf_214319");
  // strcpy(pm->femesh_dir, "canine_norm_149344");
  // strcpy(pm->femesh_dir, "canine_hf_241725");
  strcpy(pm->heart_mesh, "heart_het.dat");
  strcpy(pm->pace1_nodes, "heart_apex.dat");
  strcpy(pm->pace2_nodes, "heart_L.dat");
  strcpy(pm->torso_mesh, "torso/torso.dat");
  pm->iecg = 0;
  pm->ecgnode[0] = 1;
  pm->ecgnode[1] = 1;
  pm->ecgnode[2] = 1;
  pm->ecgnode[3] = 1;
  pm->ecgnode[4] = 1;
  pm->ecgnode[5] = 1;
  pm->ecgnode[6] = 1;
  pm->ecgnode[7] = 1;
  pm->ecgnode[8] = 1;
  */

  fclose(fp_inputdeck);
}
#endif
