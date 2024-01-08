#ifndef PARAMSET_CPP
#define PARAMSET_CPP

#include <string.h>
#include <stdlib.h>

void paramset(char* fn, param_t* pm)
{
  FILE* fp = fopen(fn, "rt");
  char str[500];
  char* ptr;
  int i;
  fgets(str, 500, fp); if(str[0]!='*')  {ptr = strtok(str, " "); ptr = strtok(NULL, " "); pm->tmax = atof(ptr);}
  fgets(str, 500, fp); if(str[0]!='*')  {ptr = strtok(str, " "); ptr = strtok(NULL, " "); pm->t_write_vtk = atof(ptr);}
  fgets(str, 500, fp); if(str[0]!='*')  {ptr = strtok(str, " "); ptr = strtok(NULL, " "); strcpy(pm->mesh_type, ptr);}
  fgets(str, 500, fp); if(str[0]!='*')  {ptr = strtok(str, " "); ptr = strtok(NULL, " "); strcpy(pm->ionmodel, ptr);}
  fgets(str, 500, fp); if(str[0]!='*')  {ptr = strtok(str, " "); ptr = strtok(NULL, " "); strcpy(pm->ic, ptr);}
  fgets(str, 500, fp); if(str[0]!='*')  {ptr = strtok(str, " "); ptr = strtok(NULL, " "); strcpy(pm->variant, ptr);}
  fgets(str, 500, fp); if(str[0]!='*')  {ptr = strtok(str, " "); ptr = strtok(NULL, " "); pm->celltype = atoi(ptr);}
  fgets(str, 500, fp); if(str[0]!='*')  {ptr = strtok(str, " "); ptr = strtok(NULL, " "); pm->dt = atof(ptr);}
  fgets(str, 500, fp); if(str[0]!='*')  {ptr = strtok(str, " "); ptr = strtok(NULL, " "); pm->num_pace1 = atoi(ptr);} 
  fgets(str, 500, fp); if(str[0]!='*')  {ptr = strtok(str, " "); ptr = strtok(NULL, " "); pm->ibclb = atoi(ptr);}
  fgets(str, 500, fp); if(str[0]!='*')  {ptr = strtok(str, " "); ptr = strtok(NULL, " "); pm->bcl_decrement = atoi(ptr);}
  fgets(str, 500, fp); if(str[0]!='*')  {ptr = strtok(str, " "); ptr = strtok(NULL, " "); pm->bcl_end = atoi(ptr);} 
  fgets(str, 500, fp); if(str[0]!='*')  {ptr = strtok(str, " "); ptr = strtok(NULL, " "); pm->bcl_pace2 = atoi(ptr);} 
  fgets(str, 500, fp); if(str[0]!='*')  {ptr = strtok(str, " "); ptr = strtok(NULL, " "); pm->stim_dur = atof(ptr);}
  fgets(str, 500, fp); if(str[0]!='*')  {ptr = strtok(str, " "); ptr = strtok(NULL, " "); pm->vth = atof(ptr);}
  fgets(str, 500, fp); if(str[0]!='*')  {ptr = strtok(str, " "); ptr = strtok(NULL, " "); pm->vtarget = atof(ptr);}
  fgets(str, 500, fp); if(str[0]!='*')  {ptr = strtok(str, " "); ptr = strtok(NULL, " "); pm->v_apd90 = atof(ptr);}
  fgets(str, 500, fp); if(str[0]!='*')  {ptr = strtok(str, " "); ptr = strtok(NULL, " "); pm->apdslop = atof(ptr);}
  fgets(str, 500, fp); if(str[0]!='*')  {ptr = strtok(str, " "); ptr = strtok(NULL, " "); strcpy(pm->stim_type, ptr);}
  fgets(str, 500, fp); if(str[0]!='*')  {ptr = strtok(str, " "); ptr = strtok(NULL, " "); pm->ifischem = atoi(ptr);}
  fgets(str, 500, fp); if(str[0]!='*')  {ptr = strtok(str, " "); ptr = strtok(NULL, " "); pm->ifbrugada = atoi(ptr);}
  fgets(str, 500, fp); if(str[0]!='*')  {ptr = strtok(str, " "); ptr = strtok(NULL, " "); pm->is_using_output = atoi(ptr);}
  fgets(str, 500, fp); if(str[0]!='*')  {ptr = strtok(str, " "); ptr = strtok(NULL, " "); pm->sf_diff = atof(ptr);} 
  fgets(str, 500, fp); if(str[0]!='*')  {ptr = strtok(str, " "); ptr = strtok(NULL, " "); pm->sf_diff_fibro = atof(ptr);} 
  fgets(str, 500, fp); if(str[0]!='*')  {ptr = strtok(str, " "); ptr = strtok(NULL, " "); strcpy(pm->fibrotic_nodes, ptr);} 
  fgets(str, 500, fp); if(str[0]!='*')  {ptr = strtok(str, " "); ptr = strtok(NULL, " "); strcpy(pm->femesh_dir, ptr);} 
  fgets(str, 500, fp); if(str[0]!='*')  {ptr = strtok(str, " "); ptr = strtok(NULL, " "); strcpy(pm->heart_mesh, ptr);} 
  fgets(str, 500, fp); if(str[0]!='*')  {ptr = strtok(str, " "); ptr = strtok(NULL, " "); strcpy(pm->pace1_nodes, ptr);} 
  fgets(str, 500, fp); if(str[0]!='*')  {ptr = strtok(str, " "); ptr = strtok(NULL, " "); strcpy(pm->pace2_nodes, ptr);} 
  fgets(str, 500, fp); if(str[0]!='*')  {ptr = strtok(str, " "); ptr = strtok(NULL, " "); pm->dt_write = atof(ptr);}
  fgets(str, 500, fp); if(str[0]!='*')  {ptr = strtok(str, " "); ptr = strtok(NULL, " "); strcpy(pm->output_mesh_type, ptr);} 
  fgets(str, 500, fp); if(str[0]!='*')  {ptr = strtok(str, " "); ptr = strtok(NULL, " "); pm->apd_fiber = atof(ptr);}
  fgets(str, 500, fp); if(str[0]!='*')  {ptr = strtok(str, " "); ptr = strtok(NULL, " "); pm->erp_fiber = atof(ptr);}
  fgets(str, 500, fp); if(str[0]!='*')  {ptr = strtok(str, " "); ptr = strtok(NULL, " "); pm->ibcl_fiber = atoi(ptr);} 
  fgets(str, 500, fp); if(str[0]!='*')  {ptr = strtok(str, " "); ptr = strtok(NULL, " "); pm->is_lbbb = atoi(ptr);}
  fgets(str, 500, fp); if(str[0]!='*')  {ptr = strtok(str, " "); ptr = strtok(NULL, " "); pm->is_rbbb = atoi(ptr);}
  fgets(str, 500, fp); if(str[0]!='*')  {ptr = strtok(str, " "); ptr = strtok(NULL, " "); pm->is_crt = atoi(ptr);}
  fgets(str, 500, fp); if(str[0]!='*')  {ptr = strtok(str, " "); ptr = strtok(NULL, " "); pm->t_crt = atoi(ptr);}
  fgets(str, 500, fp); if(str[0]!='*')  {ptr = strtok(str, " "); ptr = strtok(NULL, " "); strcpy(pm->torso_mesh, ptr);} 
  fgets(str, 500, fp); if(str[0]!='*')  {ptr = strtok(str, " "); ptr = strtok(NULL, " "); pm->iecg = atoi(ptr);}
  for(i=0; i<9; i++)
    fgets(str, 500, fp); if(str[0]!='*')  {ptr = strtok(str, " "); ptr = strtok(NULL, " "); pm->ecgnode[i] = atoi(ptr);}
  fgets(str, 500, fp); if(str[0]!='*')  {ptr = strtok(str, " "); ptr = strtok(NULL, " "); strcpy(pm->openVTK, ptr);} 
  fclose(fp);

/*  PetscPrintf(PETSC_COMM_WORLD, "-------------ParamSet------------\n");
  PetscPrintf(PETSC_COMM_WORLD, "tmax = %lf\n", pm->tmax);
  PetscPrintf(PETSC_COMM_WORLD, "twrite = %lf\n", pm->t_write_vtk);
  PetscPrintf(PETSC_COMM_WORLD, "mesh = %s\n", pm->mesh);
  PetscPrintf(PETSC_COMM_WORLD, "ionmodel = %s\n", pm->ionmodel);
  PetscPrintf(PETSC_COMM_WORLD, "ic = %s\n", pm->ic);
  PetscPrintf(PETSC_COMM_WORLD, "variant = %s\n", pm->variant);
  PetscPrintf(PETSC_COMM_WORLD, "celltype = %d\n", pm->celltype);
  PetscPrintf(PETSC_COMM_WORLD, "dt = %lf\n", pm->dt);
  PetscPrintf(PETSC_COMM_WORLD, "irepeat = %d\n", pm->num_pace1);
  PetscPrintf(PETSC_COMM_WORLD, "ibclb = %d\n", pm->ibclb);
  PetscPrintf(PETSC_COMM_WORLD, "ibcld = %d\n", pm->bcl_decrement);
  PetscPrintf(PETSC_COMM_WORLD, "ibcle = %d\n", pm->bcl_end);
  PetscPrintf(PETSC_COMM_WORLD, "stimdur = %lf\n", pm->stim_dur);
  PetscPrintf(PETSC_COMM_WORLD, "vth = %lf\n", pm->vth);
  PetscPrintf(PETSC_COMM_WORLD, "vtarget = %lf\n", pm->vtarget);
  PetscPrintf(PETSC_COMM_WORLD, "apd90 = %lf\n", pm->v_apd90);
  PetscPrintf(PETSC_COMM_WORLD, "apdslop = %lf\n", pm->apdslop);
  PetscPrintf(PETSC_COMM_WORLD, "stimtype = %s\n", pm->stim_type);
  PetscPrintf(PETSC_COMM_WORLD, "ifischem = %d\n", pm->ifischem);
  PetscPrintf(PETSC_COMM_WORLD, "ifbrugada = %d\n", pm->ifbrugada);
  PetscPrintf(PETSC_COMM_WORLD, "iopt = %d\n", pm->is_using_output);
  PetscPrintf(PETSC_COMM_WORLD, "Dscale = %lf\n", pm->sf_diff);
  PetscPrintf(PETSC_COMM_WORLD, "Dscalefib = %lf\n", pm->sf_diff_fibro);
  PetscPrintf(PETSC_COMM_WORLD, "fibrosisratio = %s\n", pm->fibrotic_nodes);
  PetscPrintf(PETSC_COMM_WORLD, "femesh_dir = %s\n", pm->femesh_dir);
  PetscPrintf(PETSC_COMM_WORLD, "heart = %s\n", pm->heart_mesh);
  PetscPrintf(PETSC_COMM_WORLD, "pace= %s\n", pm->pace1_nodes);
  PetscPrintf(PETSC_COMM_WORLD, "pace2 = %s\n", pm->pace2_nodes);
  PetscPrintf(PETSC_COMM_WORLD, "xlen = %lf\n", pm->xlen);
  PetscPrintf(PETSC_COMM_WORLD, "ylen = %lf\n", pm->ylen);
  PetscPrintf(PETSC_COMM_WORLD, "dx = %lf\n", pm->dx);
  PetscPrintf(PETSC_COMM_WORLD, "dy = %lf\n", pm->dy);
  PetscPrintf(PETSC_COMM_WORLD, "dtw = %lf\n", pm->dt_write);
  PetscPrintf(PETSC_COMM_WORLD, "outputmesh = %s\n\n", pm->output_mesh_type);
  PetscPrintf(PETSC_COMM_WORLD, "apd_nerve = %lf\n\n", pm->apd_nerve);
  PetscPrintf(PETSC_COMM_WORLD, "erp_nerve = %lf\n\n", pm->erp_nerve);
  PetscPrintf(PETSC_COMM_WORLD, "ibcl_nerve = %d\n\n", pm->ibcl_nerve);
  PetscPrintf(PETSC_COMM_WORLD, "lbbb = %d\n\n", pm->is_lbbb);
  PetscPrintf(PETSC_COMM_WORLD, "rbbb = %d\n\n", pm->is_rbbb);
  PetscPrintf(PETSC_COMM_WORLD, "crt = %d\n\n", pm->is_crt);
  PetscPrintf(PETSC_COMM_WORLD, "crttime = %d\n\n", pm->t_crt);
  PetscPrintf(PETSC_COMM_WORLD, "torso = %s\n", pm->torso_mesh);
  PetscPrintf(PETSC_COMM_WORLD, "iecg = %d\n", pm->iecg);
*/
/*
  printf("tmax = %lf\n", pm->tmax);
  printf("twrite = %lf\n", pm->t_write_vtk);
  printf("mesh = %s\n", pm->mesh);
  printf("ionmodel = %s\n", pm->ionmodel);
  printf("ic = %s\n", pm->ic);
  printf("variant = %s\n", pm->variant);
  printf("celltype = %d\n", pm->celltype);
  printf("dt = %lf\n", pm->dt);
  printf("irepeat = %d\n", pm->num_pace1);
  printf("ibclb = %d\n", pm->ibclb);
  printf("ibcld = %d\n", pm->bcl_decrement);
  printf("ibcle = %d\n", pm->bcl_end);
  printf("stimdur = %lf\n", pm->stim_dur);
  printf("vth = %lf\n", pm->vth);
  printf("vtarget = %lf\n", pm->vtarget);
  printf("apd90 = %lf\n", pm->v_apd90);
  printf("apdslop = %lf\n", pm->apdslop);
  printf("stimtype = %s\n", pm->stim_type);
  printf("ifischem = %d\n", pm->ifischem);
  printf("ifbrugada = %d\n", pm->ifbrugada);
  printf("iopt = %d\n", pm->is_using_output);
  printf("Dscale = %lf\n", pm->sf_diff);
  printf("Dscalefib = %lf\n", pm->sf_diff_fibro);
  printf("fibrosisratio = %s\n", pm->fibrotic_nodes);
  printf("femesh_dir = %s\n", pm->femesh_dir);
  printf("heart = %s\n", pm->heart_mesh);
  printf("pace= %s\n", pm->pace1_nodes);
  printf("pace2 = %s\n", pm->pace2_nodes);
  printf("xlen = %lf\n", pm->xlen);
  printf("ylen = %lf\n", pm->ylen);
  printf("dx = %lf\n", pm->dx);
  printf("dy = %lf\n", pm->dy);
  printf("dtw = %lf\n", pm->dt_write);
  printf("outputmesh = %s\n", pm->output_mesh_type);
*/
}
#endif
