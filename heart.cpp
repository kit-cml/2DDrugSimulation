#include "heart.h"

#include <algorithm>
#include <fstream>
#include <iterator>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <string>
#include <queue>
#include <omp.h>

#include "DrugSimulation/modules/commons.hpp"

bool print = false;
double max_time_step = 1.0;
double time_point = 25.0;

// constants to avoid magic values
static const char *RESULT_FOLDER_PATH = "result";


double get_cmax(char *drug_name);


Cellmodel *p_cell;  // cell model object

Cheart::Cheart(Cmesh* pmesh) : Cfem(pmesh)  // single
{
  memAllocationHeart();
}

Cheart::Cheart(Cmesh* pmesh, Cfiber* pfiber) : Cfem(pmesh)
{
  this->pfiber = pfiber;
  memAllocationHeart();
}

Cheart::Cheart(Cmesh* pmesh, Cbspm* pbspm) : Cfem(pmesh)
{
  this->pbspm = pbspm;
  memAllocationHeart();
}

Cheart::Cheart(Cmesh* pmesh, Cfiber* pfiber, Cbspm* pbspm) : Cfem(pmesh)
{
  this->pfiber = (Cfiber*) pfiber;
  this->pbspm = pbspm;
  memAllocationHeart();
}

void Cheart::memAllocationHeart()
{
  vm_glob = new double [pmesh->all_nodes];
  s1nodes = new bool [pmesh->all_nodes];
  s2nodes = new bool [pmesh->all_nodes];
  //inet = new double[pmesh->all_nodes];
  //qnet = new double[pmesh->all_nodes];
  //vm_peak = new double[pmesh->all_nodes];
  ectopic_nodes = new bool [pmesh->all_nodes];
  for(int i=0; i<pmesh->all_nodes; i++){
    s1nodes[i] = false;
    s2nodes[i]  = false;
    ectopic_nodes[i] = false;
    //inet[i] = 0.;
    //qnet[i] = 0.;
  }
  qnet = 0.;
  temp_result.vm_peak = -88.;

  Glob_cell = new double* [NUM_CELLPARAM];
  for(int i=0; i<NUM_CELLPARAM; i++)
    Glob_cell[i] = new double [pmesh->all_nodes];
    //Glob_cell[i] = new double [pmesh->all_nodes/mympi::size +1]; 
  PetscPrintf(PETSC_COMM_WORLD, "memAllocationHeart\n");
}

int Cheart::processing(int argc,char **args, param_t parameter)
{
  PetscPrintf(PETSC_COMM_WORLD, "Processing start\n");
  pm = parameter;
  dt = pm.dt;
  PetscPrintf(PETSC_COMM_WORLD, "dt: %lf\n", dt);
  int iter_dtw_vm = int(1.0/dt);  // every 5 ms
  int iter_dtw_vtk = int(pm.dt_write/dt);  // every dtw ms
  int iter_dtw_output = int(1000.0/dt);  // every 1000 ms
  int apd_count = 0;

  // get concentration from input
  char *token = strtok( pm.concs, "," );
  int idx = 0;
  while( token != NULL )
  { // begin data tokenizing
    concs[idx] = strtod(token, NULL);
    token = strtok(NULL, ",");
    PetscPrintf(PETSC_COMM_WORLD, "Concs[%d]: %lf\n", idx, concs[idx]);
    idx++;
  } // end data tokenizing
  snprintf(drug_name, sizeof(drug_name),"%s",pm.drug_name);
  PetscPrintf(PETSC_COMM_WORLD, "Drug name: %s\n", drug_name);

  // make a directory for each concentration
  // and wait till all directories have been created
  if(mympi::rank == 0){
    make_directory(RESULT_FOLDER_PATH);
    for( idx = 0; idx < 4; idx++ )
    { // begin concentration loop
      snprintf( buffer, sizeof(buffer), "%s/%.2lf", RESULT_FOLDER_PATH, concs[idx] );
      if(is_file_existed(RESULT_FOLDER_PATH) == 0) make_directory(buffer);
    } // end concentration loop
  }


for( idx = 0; idx < 4; idx++ ){
  conc = concs[idx];
  pace_count = 0;
  preprocess(); PetscPrintf(PETSC_COMM_WORLD, "preprocess\n");
  prepare_Amat_Bmat_bb_Voltsc(p_cell->STATES[0]); PetscPrintf(PETSC_COMM_WORLD, "prepare --- %lf\n", p_cell->STATES[0]);
  assemble_Amat_Bmat(dt); PetscPrintf(PETSC_COMM_WORLD, "assemble_Amat_Bmat done!\n");
  initialize_Glob_cell(pm.is_using_output); PetscPrintf(PETSC_COMM_WORLD, "initialize_Glob_cell done!\n");

  // =================================================================
  iter_bcl= int(double(pm.bcl_init)/dt);  // num of iteration for bcl
  time = 0.0;
  iter_count=0;

  while(time< pm.t_max){
#ifndef BSPMONLY
    assemble_Voltsc();
		PetscPrintf(PETSC_COMM_WORLD, "entering cellular_reaction()!\n");
    cellular_reaction();
		PetscPrintf(PETSC_COMM_WORLD, "exiting cellular_reaction()!\n");
    assemble_Voltsc();
    solve_Bmat(dt);
		PetscPrintf(PETSC_COMM_WORLD, "entering writing time-series result!\n");
    if(mympi::rank==0) write_result(iter_dtw_vm, iter_dtw_output, iter_dtw_vtk);
		PetscPrintf(PETSC_COMM_WORLD, "entering writing VTK!\n");
    if(mympi::rank==0 && pace_count >= pm.num_pace1-2 ) write_vtk(iter_dtw_vm, iter_dtw_output, iter_dtw_vtk);
		PetscPrintf(PETSC_COMM_WORLD, "exiting writing VTK!\n");
		PetscPrintf(PETSC_COMM_WORLD, "entering writing last output!\n");
    if(!(iter_count%iter_dtw_output)) writeoutput();
		PetscPrintf(PETSC_COMM_WORLD, "exiting writing last output!\n");
#else
    if(pm.is_bspm||pm.is_ecg) solveBSPM(time);
#endif
    iter_count++;
    // at the end of each BCL
    if(iter_count % (int)(pm.bcl_init/dt) == 0 && iter_count > 0){
      // make VTK for the feature
      writeVTKFeature();
      //fprintf(fpvmpeak, "%d %d %lf\n", mympi::rank+(pace_count*mympi::size),pace_count, temp_result.vm_peak);
      pace_count++;
      PetscPrintf(PETSC_COMM_WORLD, "Last qnet values %lf in node %d reset at %lf msec\n", temp_result.qnet, pmesh->target[1], time);
      //PetscPrintf(PETSC_COMM_WORLD, "Last inet and qnet values %lf %lf in node %d reset at %lf msec\n", inet, qnet, pmesh->target[1], time);
      
      for(int idx=0; idx<pmesh->all_nodes; idx++) {
       // inet[idx] = 0.;
       //qnet[idx] = 0.;
       // vm_peak[idx] = initStates[0];
      };
      qnet = 0.;
      temp_result.init( initStates[0], initStates[9] );
    }
    time += dt;
  } // end of while
  PetscPrintf(PETSC_COMM_WORLD, "FINISHED LOOPING\n");
  if (mympi::rank == 0) {
    fclose(fpvm);
    fclose(fpqnet);
    //fclose(fpcai);
    //fclose(fpvmpeak);
    PetscPrintf(PETSC_COMM_WORLD, "FINISHED 3D SIM!!!\n");
  }
    delete []initStates;
}
  return 0;
}

void Cheart::preprocess()
{
  idid=false, idid2=false, idid3=false;
  PetscPrintf( PETSC_COMM_WORLD, "Start preprocssing\n");
  loadPacingsites(pm.pace1_nodes, pm.pace2_nodes);  // S1S2 protocol
  pmesh->dimension = pm.dimension;
  // Loading surface mesh 
  if(!strcmp(pmesh->type, "tet") || !strcmp(pmesh->type, "hex"))
    if(!strcmp(pm.output_mesh_type, "surface")) 
      pmesh_surf = new Cmesh(pm.surface_mesh);

  // CellML model initialization -> pcell

#ifdef TOMEK_2019
  PetscPrintf(PETSC_COMM_WORLD, "Using Tomek cell model\n");
  p_cell = new Tomek_model();
#else
  PetscPrintf(PETSC_COMM_WORLD, "Using O'Hara Rudy cell model\n");
  p_cell = new Ohara_Rudy_2011();
#endif
  if( p_cell == NULL ){
    PetscPrintf(PETSC_COMM_WORLD, "Problem when initializing cellmodel\n");
  }

  PetscPrintf( PETSC_COMM_WORLD, "Test drug file name: %s\n", pm.hill_file );
  ic50 = get_IC50_data_from_file(pm.hill_file);
  sample_id = 0;
  for( int idx = 0; idx < 14; idx++ ){
    PetscPrintf(PETSC_COMM_WORLD, "%lf|", (ic50[sample_id]).data[idx]);
  }
  PetscPrintf(PETSC_COMM_WORLD,"\n");
#ifdef TOMEK_2019
  p_cell->initConsts( pm.celltype, conc, ic50[sample_id].data);
#else
  p_cell->initConsts( pm.celltype, conc, ic50[sample_id].data, true);
#endif
  PetscPrintf(PETSC_COMM_WORLD, "BEFORE: Stimamp: %lf Stimdur: %lf\n", 
                                p_cell->CONSTANTS[amp], 
                                p_cell->CONSTANTS[duration]);
  PetscPrintf(PETSC_COMM_WORLD, "BEFORE: GNa: %lf\n", p_cell->CONSTANTS[GNa]);
  p_cell->CONSTANTS[duration] = pm.stim_dur;
  p_cell->CONSTANTS[amp] *= pm.stim_amp;
  p_cell->CONSTANTS[BCL] = pm.bcl_init;
  p_cell->CONSTANTS[GNa] *= pm.gna_scale;
  PetscPrintf(PETSC_COMM_WORLD, "AFTER: Stimamp: %lf Stimdur: %lf\n", 
                                p_cell->CONSTANTS[amp], 
                                p_cell->CONSTANTS[duration]);
  PetscPrintf(PETSC_COMM_WORLD, "AFTER: GNa: %lf\n", p_cell->CONSTANTS[GNa]);
  PetscPrintf(PETSC_COMM_WORLD, "BCL: %lf\n", pm.bcl_init);

  // Save intitial States[]
  initStates = new double[p_cell->states_size];
  for(int i=0; i<p_cell->states_size; i++) initStates[i] = p_cell->STATES[i];

  assign_conductivity();  // concuctivity

  // Determine what scalars to vtk output files.
  char *token, *buff;
  token = strtok_r(pm.vtkout_states, "," , &buff);
  while( token != NULL ){
      vtkout_states_vec.push_back( token );
      PetscPrintf( PETSC_COMM_WORLD, "%s\n", token );
      token  = strtok_r(NULL, "," , &buff);
  }
  token = strtok_r(pm.vtkout_algebraic, "," , &buff);
  while( token != NULL ){
      vtkout_algebraic_vec.push_back( token );
      PetscPrintf( PETSC_COMM_WORLD, "%s\n", token );
      token  = strtok_r(NULL, "," , &buff);
  }
  PetscPrintf( PETSC_COMM_WORLD, "STATES VTKOUT: %d\n", vtkout_states_vec.size());
  PetscPrintf( PETSC_COMM_WORLD, "ALGEBRAIC VTKOUT: %d\n", vtkout_algebraic_vec.size());
  // Find purkinje synampses, when sinus rhythm.
  if(!strcmp(pm.mesh_type, "tetfiber")){
    FindSynapses();
    if(pm.is_crt == 0) idid= true; // idid:1 --> No apex stimmulation
  }
  // Open plt files.
  if(mympi::rank==0) {
    snprintf(buffer, sizeof(buffer), "%s/%.2f/vmcheck_%s_%.2lf.plt", RESULT_FOLDER_PATH, conc, pm.drug_name, conc);
    fpvm = fopen(buffer, "wt");
    snprintf(buffer, sizeof(buffer), "%s/%.2f/qnet_%s_%.2lf.plt", RESULT_FOLDER_PATH, conc, pm.drug_name, conc);
    fpqnet = fopen(buffer, "wt");
    snprintf(buffer, sizeof(buffer), "%s/%.2f/cai_%s_%.2lf.plt", RESULT_FOLDER_PATH, conc, pm.drug_name, conc);
    //fpcai = fopen(buffer, "wt");
    snprintf(buffer, sizeof(buffer), "%s/%.2f/vmpeak_%s_%.2lf.plt", RESULT_FOLDER_PATH, conc, pm.drug_name, conc);
    //fpvmpeak = fopen(buffer, "wt");
 }
  //snprintf(buffer, sizeof(buffer), "%s/%.2f/vmpeak_%s_%.2lf_proc%d.plt", RESULT_FOLDER_PATH, conc, pm.drug_name, conc, mympi::rank);
  //fpvmpeak = fopen(buffer, "wt");

  assign_nodes_elements_proc();

}

void Cheart::cellular_reaction()
{
  int i, j, iia;
  

  if(!strcmp(pm.mesh_type, "tetfiber")) { // Sinus rhythm
    if(!(iter_count%iter_bcl)) {
      pfiber->AVpacing();
//      PetscPrintf(PETSC_COMM_WORLD, "AVPacing executed at: %lf\n", time );
    }
    pfiber->PropagationPF();
    StimTissue();
  }
  else if( (!strcmp(pm.mesh_type, "tet") || 
           !strcmp(pm.mesh_type, "fiber_real")) &&
           idid == false ) { // apex pacing for S1-S2 protocol
    if( (time -  floor(time/p_cell->CONSTANTS[BCL]) * p_cell->CONSTANTS[BCL]>=p_cell->CONSTANTS[stim_start]) && (time -  floor(time/p_cell->CONSTANTS[BCL])*p_cell->CONSTANTS[BCL]<=p_cell->CONSTANTS[stim_start]+p_cell->CONSTANTS[duration]) ){
      for(i=0; i<num_s1nodes; i++){
        s1nodes[s1_num[i]] = true;
      }
    }
    else{
      for(i=0; i<num_s1nodes; i++){
        s1nodes[s1_num[i]] = false;
      }
    }
  }
  else if( (!strcmp(pm.mesh_type, "square") || !strcmp(pm.mesh_type, "fiber_real")) && idid == false ) {
    for(i=0; i<num_s1nodes; i++){
      s1nodes[s1_num[i]] = false;
    }
    iia=(int)(time/p_cell->CONSTANTS[BCL]);
    if( time < (double)iia*p_cell->CONSTANTS[BCL] + p_cell->CONSTANTS[duration] + pm.t_crt &&
        time > (double)iia*p_cell->CONSTANTS[BCL] + pm.t_crt ){
      for(i=0; i<num_s1nodes; i++){
        s1nodes[s1_num[i]] = true;
      }
    }
  }


  if(int(time/pm.bcl_init) >= pm.num_pace1){
    for(int i=0; i<num_s1nodes; i++)
      s1nodes[s1_num[i]] = false;
  }


  if(int(time/pm.bcl_init) == pm.num_pace1-1)
    if(pm.is_s1s2==1) apply_s2();
    else if(pm.is_s1s2==2) apply_ectopic();


	PetscPrintf( PETSC_COMM_WORLD, "Preparing cellular node loop\n");
  for(i=0; i<pmesh->all_nodes; i++){
    if(node_proc[i]==mympi::rank){
      ilocal=local_nodes[i];
      for(j=0; j< p_cell->states_size; j++){
        p_cell->STATES[j] =Glob_cell[j][ilocal];
      }
      p_cell->STATES[0]=Voltsc[0][i];

//PetscPrintf( PETSC_COMM_WORLD, "Debug point 0\n"); //PASSED
//exit(0);         

      p_cell->isS1 = s1nodes[i];
      p_cell->isEctopic = ectopic_nodes[i];
// call the ODE function and solve it with analytical method
      p_cell->computeRates(time,
                           p_cell->CONSTANTS,
                           p_cell->RATES,
                           p_cell->STATES,
                           p_cell->ALGEBRAIC
                          );
// compute accepted timestep
#ifdef ON_CONSTRUCTION
      dt_set = Ohara_Rudy_2011::set_time_step(time,
                           time_point,
                           max_time_step,
                           p_cell->CONSTANTS,
                           p_cell->RATES,
                           p_cell->STATES,
                           p_cell->ALGEBRAIC);

//PetscPrintf( PETSC_COMM_WORLD, "Debug point 1\n");
//exit(0);


      if (floor((time + dt_set) / pm.bcl_init) == floor(time / pm.bcl_init)) {
        dt = dt_set;
      }
      else {
        dt = (floor(time / pm.bcl_init) + 1) * pm.bcl_init - time;
      }
#endif
      p_cell->solveAnalytical(dt);
      if(p_cell->STATES[V] > -88.0){
      qnet += (p_cell->ALGEBRAIC[INaL]+p_cell->ALGEBRAIC[ICaL]+p_cell->ALGEBRAIC[Ito]+p_cell->ALGEBRAIC[IKr]+p_cell->ALGEBRAIC[IKs]+p_cell->ALGEBRAIC[IK1])*dt;
      }

//PetscPrintf( PETSC_COMM_WORLD, "Debug point 2\n"); // PASSED
//exit(0);

/*
      // applying other 14 features
      // (move single cell code parts to here) 

      // vmpeak test
      // get the peak Vm 6 secs after depolarization (when Na channel just closed after bursting)
       if( iter_count % (int)(((p_cell->CONSTANTS[BCL]*pace_count)+(p_cell->CONSTANTS[stim_start]+6.)) / dt) == 0 && temp_result.vm_peak < p_cell->STATES[V] ){
        temp_result.vm_peak = p_cell->STATES[V];
      }
      if( temp_result.ca_peak < p_cell->STATES[cai] ){
          temp_result.ca_peak = p_cell->STATES[cai];
          ca_amp50 = temp_result.ca_peak - (0.5 * (temp_result.ca_peak - temp_result.ca_valley));
          ca_amp90 = temp_result.ca_peak - (0.9 * (temp_result.ca_peak - temp_result.ca_valley));
          t_ca_peak = time;

      }
*/     
      
      for(j=0; j< p_cell->states_size ; j++){
        Glob_cell[j][ilocal]= p_cell->STATES[j];
      }

//PetscPrintf( PETSC_COMM_WORLD, "Debug point 3\n"); //PASSED
//exit(0);

// SUSPECTED ERRORS AROUND HERE BEGIN
      VecSetValues(Volt[0],1,&i,&p_cell->STATES[V],INSERT_VALUES);
      temp_result.qnet = qnet/1000.;
      VecSetValues(Volt[1],1,&i,&temp_result.qnet,INSERT_VALUES);
#ifdef TEMP
      VecSetValues(Volt[2],1,&i,&temp_result.vm_peak,INSERT_VALUES);
      VecSetValues(Volt[3],1,&i,&temp_result.ca_peak,INSERT_VALUES);
      if(iter_count % (int)(pm.bcl_init/dt) == 0 && iter_count > 0){
        temp_result.vm_dia = p_cell->STATES[V];
        temp_result.ca_dia = p_cell->STATES[cai];
        VecSetValues(Volt[4],1,&i,&temp_result.apd90,INSERT_VALUES);
        VecSetValues(Volt[5],1,&i,&temp_result.apd50,INSERT_VALUES);
      }
#endif
// SUSPECTED ERRORS AROUND HERE END
//PetscPrintf( PETSC_COMM_WORLD, "Debug point 4\n");
//exit(0);


      //VecSetValues(Volt[1],1,&i,&p_cell->STATES[cai],INSERT_VALUES);
#ifdef TEMP
      for(j=0; j<vtkout_states_vec.size(); j++) 
        VecSetValues(Volt[j],1,&i,&p_cell->STATES[p_cell->statesMap[vtkout_states_vec[j]]],INSERT_VALUES);
      for(j=0; j<vtkout_algebraic_vec.size(); j++) 
        VecSetValues(Volt1[j],1,&i,&p_cell->ALGEBRAIC[p_cell->algebraicMap[vtkout_algebraic_vec[j]]],INSERT_VALUES);
#endif
    }
  }
	PetscPrintf( PETSC_COMM_WORLD, "Finished cellular loop\n");
}


void Cheart::StimTissue()
{
  int i, ipp;
  for(i=0; i<pmesh->all_nodes; i++) s1nodes[i]=false;
  for(i=0; i<pfiber->pmesh->all_terminals; i++){  // 모든 말단노드에 대하여
    ipp=pfiber->pmesh->iPFnode[i];  // i번째 말단노드의 퍼킨지메쉬 전역노드번호
    if(pfiber->Volt[ipp]==Cfiber::DEPOL){ // i번째 말단노드가 흥분하였다면.
      if(pfiber->Volt_t[ipp] < p_cell->CONSTANTS[duration]){ // 말단노드가 흥분한지 stimdur 이내라면,
//          PetscPrintf(PETSC_COMM_WORLD, "TEST STIMTISSUE %lf %lf\n", pfiber->Volt_t[ipp], time);
        s1nodes[pfiber->pmesh->inearsurf[i]] = true;  // 심실표면에서 가장가까운 노드에 자>극을 줌 stimulus in mV
        //PetscPrintf(PETSC_COMM_WORLD, "StimTissue executed at %d\n", pfiber->pmesh->inearsurf[i]);
      }
    }
  }
}




void Cheart::apply_ectopic()
{
  int i, j;
  if(Voltsc[0][s2_num[0]] > -75) {
    idid2 = true; // to depol
    PetscPrintf(PETSC_COMM_WORLD, "time: %lf idid2 true!\n", time);
  }
  if(idid2==true && Voltsc[0][s2_num[0]]< -75){
    idid3 = true; // to repol
    PetscPrintf(PETSC_COMM_WORLD, "time: %lf idid3 true!\n", time);
  }
  if(idid3==true && idid==false) { 
    ectopic_start_time = time;
    for(int i=0; i<num_s2nodes; i++) ectopic_nodes[s2_num[i]] = true; 
    idid = true;
    PetscPrintf(PETSC_COMM_WORLD, "time: %lf idid true  ectopic true!\n", time);
  }
  if(idid && time - ectopic_start_time > p_cell->CONSTANTS[duration]){
    for(int i=0; i<num_s2nodes; i++) ectopic_nodes[s2_num[i]] = false; 
    PetscPrintf(PETSC_COMM_WORLD, "time: %lf ectopic false again!\n", time);
  }
}

void Cheart::apply_s2()
{
  // stim start, period, duration, amplitude - grandi ?
  int i, j;
  if(Voltsc[0][pmesh->target[1]] > -45) idid2 = true;
  if(idid2==true && Voltsc[0][pmesh->target[1]]< -45) idid3 = true;

  if(idid3==true && idid==false) { // S2 protocol
    PetscPrintf(PETSC_COMM_WORLD, "s2 start~~!\t\ttime: %lf   idid==false\n", time);
    for (i=0; i<pmesh->all_nodes; i++) {
      if (s2nodes[i]){
        if (node_proc[i] == mympi::rank) {
          ilocal = local_nodes[i];
          for (j = 0; j < p_cell->states_size; j++) {
            Glob_cell[j][ilocal] = initStates[j];
          }
        }
        VecSetValues(Volt[0], 1, &i, &(initStates[0]), INSERT_VALUES);
      }
    }
    idid = true;
    // This one is MANDATORY
    // Otherwise, the S2 initial rotation will have strange result
    assemble_Voltsc();
    PetscPrintf(PETSC_COMM_WORLD, "s2 over~~!\t\ttime: %lf   idid==true\n", time);
  }
}


void Cheart::write_result(int iter_dtw_vm, int iter_dtw_output, int iter_dtw_vtk)
{
  if(!(iter_count%iter_dtw_vm)){
    //fprintf(fpvm, "%lf %lf %lf %lf %lf %lf\n", time, Voltsc[0][s1_num[0]], Voltsc[0][pmesh->target[0]], Voltsc[0][pmesh->target[1]], Voltsc[0][pmesh->target[2]], Voltsc[0][s2_num[0]]);
    fprintf(fpvm, "%lf %lf %lf %lf %lf\n", time, Voltsc[0][s1_num[0]], Voltsc[0][pmesh->target[0]], Voltsc[0][pmesh->target[1]], Voltsc[0][pmesh->target[2]]);
    fprintf(fpqnet, "%lf %lf %lf %lf %lf\n", time, Voltsc[1][s1_num[0]], Voltsc[1][pmesh->target[0]], Voltsc[1][pmesh->target[1]], Voltsc[1][pmesh->target[2]]);
    // fprintf(fpvmpeak, "%lf %lf %lf %lf %lf\n", time, Voltsc[2][s1_num[0]], Voltsc[2][pmesh->target[0]], Voltsc[2][pmesh->target[1]], Voltsc[2][pmesh->target[2]]);
   fflush(fpvm);
    fflush(fpqnet);
    //fflush(fpvmpeak);
    //fflush(fpcai);
    if(!strcmp(pm.mesh_type, "tetfiber")) fprintf(pfiber->fpvm, "%lf %lf\n", time, pfiber->Volt[10]);
  }
}



void Cheart::write_vtk(int iter_dtw_vm, int iter_dtw_output, int iter_dtw_vtk)
{

#if defined VTK_COMBINED
  if(time >= pm.t_write_vtk && !(iter_count%iter_dtw_vtk)) {
    writeVTK();
  }
#else
  if( time <= 0.0 ){
    writeVTK();
  }
  if(time >= pm.t_write_vtk && !(iter_count%iter_dtw_vtk)) {
    PetscPrintf(PETSC_COMM_WORLD, "Writing at time %lf and iter_count%%iter_dtw_vtk= %d %% %d = %d!!!\n", 
          time, iter_count, iter_dtw_vtk, iter_count%iter_dtw_vtk);
    writeScalar();
  }
#endif
}


void Cheart::solveBSPM(double time)
{
   // read V_001.vtk  based on time
   // save scalar of V_001.vtk into vm_glob[]
  char vtk_name[255];
  char buff[255];
  bool is_voltage_found = false;
  double buff_val;
  sprintf(vtk_name, "./heart_results_sinus/V_%.0lf.vtk", time);
  PetscPrintf(PETSC_COMM_WORLD, "read %s\n", vtk_name);
  FILE *fp_vtk = fopen( vtk_name, "r" );
  while( fgets( buff, 255,fp_vtk) != NULL )
    if( strstr( buff, "LOOKUP" ) != NULL) break;
  for(int i=0; i<pmesh->all_nodes; i++){
    fscanf(fp_vtk, "%lf\n", &vm_glob[i]);
  }
  fclose(fp_vtk);
  PetscPrintf(PETSC_COMM_WORLD, "finished reading %s\n", vtk_name);
  // ========== compute ECG at the time==============
  PetscPrintf( PETSC_COMM_WORLD, "time=%lf\n", time);

  if(pm.is_ecg){
    for(int i=0; i<9; i++){
      pbspm->ecgl.vnode[i] = Radiate(pbspm->pmesh->XY[pm.ecgnode[i]-1]);
    }
    //if(mympi::rank==0) 
    pbspm->writeECG(time);
  }
  // =========== compute BSPM at the time ============
  if(pm.is_bspm){
    PetscPrintf(PETSC_COMM_WORLD, "BSPM being solved (time: %lf)\n", time);
    PetscPrintf(PETSC_COMM_WORLD, "BSPM all_nodes: %d\n", pbspm->pmesh->all_nodes);
    for(int i=0; i<pbspm->pmesh->all_nodes; i++){// torso nodes
      pbspm->volt[i] = Radiate(pbspm->pmesh->XY[i]);
    }
    pbspm->writeBSPM(time);
  }
}

void Cheart::initialize_Glob_cell(int iopt)
{
  int i, j;
  // ============================= Newly Added Part ===================================
  if(iopt==0){
    for(i=0; i<pmesh->all_nodes; i++){
      if(node_proc[i]==mympi::rank){
        ilocal=local_nodes[i];
	for(j=0; j< p_cell->states_size; j++){
          Glob_cell[j][ilocal]= p_cell->STATES[j];
        }
      }
    }
    PetscPrintf(PETSC_COMM_WORLD, "iopt:0 applied\n");
  }
  else{
    FILE* fpoutput;
    char fname[200];
    idid=true; // to avoid additional S2 stimulation (if no sinus rhythm)
    sprintf(fname,"output.%d",mympi::rank);
    fpoutput = fopen(fname,"r");
    fscanf(fpoutput,"%lf ",&time);
    for(i=0; i<pmesh->all_nodes; i++){
      if(node_proc[i]==mympi::rank){
        ilocal=local_nodes[i];
        for(j=0; j< p_cell->states_size; j++){
          fscanf(fpoutput,"%lf ", &p_cell->STATES[j] );
          Glob_cell[j][ilocal] = p_cell->STATES[j];
        }
        VecSetValues(Volt[0],1,&i,&(p_cell->STATES[0]),INSERT_VALUES);
        // VecSetValues(Cati,1,&i,&( NV_Ith_S(cellModel->getStates(), get_index_from_map( cellModel->getStatesMap(),cellModel->CAI_KEY))  ),INSERT_VALUES);
      }
    }
    fclose(fpoutput);
    assemble_Voltsc();
    PetscPrintf(PETSC_COMM_WORLD, "iopt:1 applied\n");
  }
}

void Cheart::assemble_Voltsc()
{
  for(int i=0; i<vtkout_states_vec.size(); i++){
//  for(int i=0; i < 2; i++){
    VecAssemblyBegin(Volt[i]);
    VecAssemblyEnd(Volt[i]);
    VecScatterBegin(ctx[i],Volt[i],V_SEQ[i],INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(ctx[i],Volt[i],V_SEQ[i],INSERT_VALUES,SCATTER_FORWARD);
    VecGetArray(V_SEQ[i],&Voltsc[i]);
  }
  for(int i=0; i<vtkout_algebraic_vec.size(); i++){
//  for(int i=0; i<0; i++){
    VecAssemblyBegin(Volt1[i]);
    VecAssemblyEnd(Volt1[i]);
    VecScatterBegin(ctx1[i],Volt1[i],V_SEQ1[i],INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(ctx1[i],Volt1[i],V_SEQ1[i],INSERT_VALUES,SCATTER_FORWARD);
    VecGetArray(V_SEQ1[i],&Voltsc1[i]);
  }

}

double Cheart::Radiate(double* electrode)
{
  int i, j, iel, k, ia, gauss;
  double bsp_one, bsp_all=0.0;
  double r[3], centroid[3], fibori[3], RR;
  double Dpar = DC*pm.sf_diff*0.8;
  double Dnor = Dpar; //Dpar/7.0;
  double delD = Dnor-Dpar;
  double w = 1.0/6.0;

  for(i=0; i<3; i++) fibori[i]=0.0;
  for(i=0; i<pmesh->dimension; i++) fibori[i] += 1./float(pmesh->node);

  // For all elements of heart
  for(iel=0; iel<pmesh->all_elements; iel++){ 
    if(element_proc[iel]==mympi::rank){
      bsp_one=0.;
      //==== Get RR (displacement btw heart node and torso node) ==========
      set_xyele(iel, xyele); // set local coordinates
      for(i=0;i<pmesh->dimension;i++) centroid[i]=0.0;
      for(i=0;i<pmesh->node;i++){
        j = pmesh->EC[iel][i];
        vm_ele[i][0] = vm_glob[j];
        //if(iel==1) printf("vm_glob[%d]: %lf\n", j, vm_glob[j]);
        for(k=0; k<pmesh->dimension; k++){
          centroid[k] += xyele[i][k]/float(pmesh->node);
        }
      }
      for(i=0; i<3; i++) r[i] = centroid[i] - electrode[i]; // displacement vector
      RR = sqrt(r[x]*r[x] + r[y]*r[y] + r[z]*r[z]); // distance between point of ventricle and point of torso
      //==== Get RR (displacement btw heart node and torso node) ==========
      //set_Jacob(xyele, dPhiS, Jacob); //make_dj(sourcemesh->node);  // make dj <--- xyele, dn
      for(gauss; gauss<pmesh->node; gauss++){ 
        MultMatMat(pmesh->dimension, pmesh->node, pmesh->dimension, dPhiSt[gauss], xyele, Jacob[gauss]);
        inverseMat(Jacob[gauss], xJ[gauss], &detJ[gauss]); //make_xdj(); // xdj <-- dj
        Transition(pmesh->dimension, pmesh->dimension, xJ[gauss], xJt[gauss]);
        MultMatMat(pmesh->node, pmesh->dimension, pmesh->dimension, dPhiS[gauss], xJt[gauss], dPhiX[gauss]);

        //MultMatVec(dPhiXt[gauss], vm_ele, vd);
        MultMatMat(pmesh->node, pmesh->dimension, 1, dPhiXt[gauss], vm_ele, vd);
        //MultMatVec(Dtensor[iel], vd, vel);  // vel = D*vd
        MultMatMat(pmesh->dimension, pmesh->dimension, 1, Dtensor[iel], vd, vel);

        bsp_one = bsp_one - (vel[x][0]*r[x] + vel[y][0]*r[y] + vel[z][0]*r[z])*detJ[gauss]*w/pow(RR,3);
      }
      VolumeConductor -= bsp_one;
      //printf("bsp_one = %lf\n", bsp_one);
    }
  }
  return VolumeConductor;
}

/*
double Cheart::Radiate(double* electrode)
{
  int i, j, iel, k, ia;
  double zout, zzz=0.0, r[3], vel[3], vd[3], centroid[3], fibori[3], displacement;
  double Dpar = DC*pm->sf_diff*0.8;
  double Dnor = Dpar; //Dpar/7.0;
  double delD = Dnor-Dpar;
  double w = 1.0/6.0;
  for(iel=0; iel<pmesh->all_nodes; iel++){  // heart elements
    for(i=0; i<3; i++) r[i] = pmesh->XY[iel][i] - electrode[i];
    displacement = sqrt(r[x]*r[x] + r[y]*r[y] + r[z]*r[z]); // distance between point of ventricle and point of tors    o
    zzz = zzz - (vm_glob[iel]*r[x] + vm_glob[iel]*r[y] + vm_glob[iel]*r[z])*w/(displacement*displacement*displacement);
  }
  return zzz;
}
*/




/*
bool Cheart::setfibrosis(param_t pm)
{
  int num, n1, i, j;
  char fname[200];
  FILE *fpfib;
  // sprintf(fname, "%s/%s", pm.femesh_dir, pm.fibrotic_nodes);
  sprintf(fname, "%s", pm.fibrotic_nodes);
  if(fpfib = fopen(fname, "rt")) PetscPrintf(PETSC_COMM_WORLD, "fibrotic node file has been opened successfully !%s\n", fname);
  else{
    PetscPrintf(PETSC_COMM_WORLD, "setfibrosis(): failed to open %s\n", fname);
    exit(-1);
  }
  fscanf(fpfib, "%d\n", &num);
  PetscPrintf(PETSC_COMM_WORLD, "num = %d\n", num);
  PetscPrintf(PETSC_COMM_WORLD, "dscalefib = %lf\n", pm.sf_diff_fibro);
  for(i=0; i<num; i++){
    fscanf(fpfib, "%d\n", &n1);
    //if(mympi::rank==0) printf("Infarct node: %d\n", n1);
    // printf("Infarct node: %d\n", n1);
    //PetscPrintf(PETSC_COMM_WORLD, "Infarct node: %d\n", n1);
    n1--;
    for(j=0; j<pmesh->node; j++) sigma_in[pmesh->EC[n1][j]][0] *= pm.sf_diff_fibro;
    //sigma_in[n1][0] *= pm.sf_diff_fibro;
    //   p_cell->nablocker=0.01;
  }
  fclose(fpfib);
  return true;
}
*/
void Cheart::writeVTK()
{
  char fname[200];
  double stepper;
  double mul;

  PetscPrintf(PETSC_COMM_WORLD, "ENTER VTK CREATE!!!\n");  
  stepper = dt >= 0.5 ? 0. : 0.5;
  sprintf(fname, "%s/%.2f/V_%s_%.2lf_%d.vtk", RESULT_FOLDER_PATH, conc, drug_name, conc, (int)floor(time+stepper));
    FILE* fpvtk=fopen(fname,"wt");
    int i, j, k;
    fprintf(fpvtk, "# vtk DataFile Version 3.0\n");
    fprintf(fpvtk, "vtk output\n");
    fprintf(fpvtk, "ASCII\n");
    fprintf(fpvtk, "DATASET UNSTRUCTURED_GRID\n");
    fprintf(fpvtk, "POINTS %d float\n", pmesh->all_nodes);
    for(i=0; i<pmesh->all_nodes; i++)
      fprintf(fpvtk, "%lf %lf %lf\n", pmesh->XY[i][x], pmesh->XY[i][y], pmesh->XY[i][z]);
    if(!strcmp(pmesh->type, "square")){
      fprintf(fpvtk, "CELLS %d %d\n", pmesh->all_elements, pmesh->all_elements*5);
      for(i=0; i<pmesh->all_elements; i++)
        fprintf(fpvtk, "4 %d %d %d %d\n", pmesh->EC[i][0], pmesh->EC[i][1], pmesh->EC[i][2], pmesh->EC[i][3]);
      fprintf(fpvtk, "CELL_TYPES %d\n", pmesh->all_elements);
      for(i=0; i<pmesh->all_elements; i++)
        fprintf(fpvtk, "9\n");  // 10 : tetrahedron
    }
    else if(!strcmp(pmesh->type, "tet")){
      if(!strcmp(pm.output_mesh_type, "surface")){
        fprintf(fpvtk, "CELLS %d %d\n", pmesh_surf->all_elements, pmesh_surf->all_elements*4);
        for(i=0; i<pmesh_surf->all_elements; i++)
          fprintf(fpvtk, "%d %d %d %d\n", pmesh_surf->node, pmesh_surf->EC[i][0], pmesh_surf->EC[i][1], pmesh_surf->EC[i][2]);
        fprintf(fpvtk, "CELL_TYPES %d\n", pmesh_surf->all_elements);
        for(i=0; i<pmesh_surf->all_elements; i++)
          fprintf(fpvtk, "5\n");  // 5 : triangle
      }
      else {
        fprintf(fpvtk, "CELLS %d %d\n", pmesh->all_elements, pmesh->all_elements*5);
        for(i=0; i<pmesh->all_elements; i++)
          fprintf(fpvtk, "4 %d %d %d %d\n", pmesh->EC[i][0], pmesh->EC[i][1], pmesh->EC[i][2], pmesh->EC[i][3]);
        fprintf(fpvtk, "CELL_TYPES %d\n", pmesh->all_elements);
        for(i=0; i<pmesh->all_elements; i++)
          fprintf(fpvtk, "10\n");  // 10 : tetrahedron
      }
    }
    else if(!strcmp(pmesh->type, "fiber_real")){
      fprintf(fpvtk, "CELLS %d %d\n", pmesh->all_elements, pmesh->all_elements*3);
      for(i=0; i<pmesh->all_elements; i++){
        fprintf(fpvtk, "2 %d %d\n", pmesh->EC[i][0], pmesh->EC[i][1]);
      }
      fprintf(fpvtk, "CELL_TYPES %d\n", pmesh->all_elements);
      for(i=0; i<pmesh->all_elements; i++)
        fprintf(fpvtk, "3\n");  // 3: fiber element
    }
    fprintf(fpvtk, "POINT_DATA %d\n", pmesh->all_nodes);
#if defined VTK_COMBINED
   for(int i =0; i<1; i++){
     fprintf(fpvtk, "SCALARS %s float\n", vtkout_states_vec[i].c_str());
     fprintf(fpvtk, "LOOKUP_TABLE default\n");
     if(strcmp(vtkout_states_vec[i].c_str(), "V") == 0){
       mul = 10.0;
     }
     else{
       mul = 10000000.0;
     }
     for(j=0; j<pmesh->all_nodes; j++)
      //fprintf(fpvtk, "%d\n", (int)(Voltsc[i][j]*mul));  // Vm
      fprintf(fpvtk, "%lf\n", (Voltsc[i][j]));  // Vm
   }
#ifdef TEMPTEMNP
   for(int i =0; i<vtkout_algebraic_vec.size(); i++){
     fprintf(fpvtk, "SCALARS %s float\n", vtkout_algebraic_vec[i].c_str());
     fprintf(fpvtk, "LOOKUP_TABLE default\n");
     for(j=0; j<pmesh->all_nodes; j++)
       //fprintf(fpvtk, "%d\n", (int)(Voltsc1[i][j]*1000.0));  // Vm
       fprintf(fpvtk, "%lf\n", (Voltsc1[i][j]));  // Vm
   }
#endif
   fclose(fpvtk);
   PetscPrintf(PETSC_COMM_WORLD, "%s has been successfully generated!!!\n", fname);
#endif

  if(!strcmp(pm.mesh_type, "tetfiber")) pfiber->writePF(time);
}

void Cheart::writeVTKFeature()
{
  char fname[200];
  double stepper;
  double mul;

  PetscPrintf(PETSC_COMM_WORLD, "ENTER VTK CREATE!!!\n");  
  stepper = dt >= 0.5 ? 0. : 0.5;
  sprintf(fname, "%s/%.2f/FEATURE_%s_%.2lf_%d.vtk", RESULT_FOLDER_PATH, conc, drug_name, conc, (int)floor(time+stepper));
    FILE* fpvtk=fopen(fname,"wt");
    int i, j, k;
    fprintf(fpvtk, "# vtk DataFile Version 3.0\n");
    fprintf(fpvtk, "vtk output\n");
    fprintf(fpvtk, "ASCII\n");
    fprintf(fpvtk, "DATASET UNSTRUCTURED_GRID\n");
    fprintf(fpvtk, "POINTS %d float\n", pmesh->all_nodes);
    for(i=0; i<pmesh->all_nodes; i++)
      fprintf(fpvtk, "%lf %lf %lf\n", pmesh->XY[i][x], pmesh->XY[i][y], pmesh->XY[i][z]);
    if(!strcmp(pmesh->type, "square")){
      fprintf(fpvtk, "CELLS %d %d\n", pmesh->all_elements, pmesh->all_elements*5);
      for(i=0; i<pmesh->all_elements; i++)
        fprintf(fpvtk, "4 %d %d %d %d\n", pmesh->EC[i][0], pmesh->EC[i][1], pmesh->EC[i][2], pmesh->EC[i][3]);
      fprintf(fpvtk, "CELL_TYPES %d\n", pmesh->all_elements);
      for(i=0; i<pmesh->all_elements; i++)
        fprintf(fpvtk, "9\n");  // 10 : tetrahedron
    }
    else if(!strcmp(pmesh->type, "tet")){
      if(!strcmp(pm.output_mesh_type, "surface")){
        fprintf(fpvtk, "CELLS %d %d\n", pmesh_surf->all_elements, pmesh_surf->all_elements*4);
        for(i=0; i<pmesh_surf->all_elements; i++)
          fprintf(fpvtk, "%d %d %d %d\n", pmesh_surf->node, pmesh_surf->EC[i][0], pmesh_surf->EC[i][1], pmesh_surf->EC[i][2]);
        fprintf(fpvtk, "CELL_TYPES %d\n", pmesh_surf->all_elements);
        for(i=0; i<pmesh_surf->all_elements; i++)
          fprintf(fpvtk, "5\n");  // 5 : triangle
      }
      else {
        fprintf(fpvtk, "CELLS %d %d\n", pmesh->all_elements, pmesh->all_elements*5);
        for(i=0; i<pmesh->all_elements; i++)
          fprintf(fpvtk, "4 %d %d %d %d\n", pmesh->EC[i][0], pmesh->EC[i][1], pmesh->EC[i][2], pmesh->EC[i][3]);
        fprintf(fpvtk, "CELL_TYPES %d\n", pmesh->all_elements);
        for(i=0; i<pmesh->all_elements; i++)
          fprintf(fpvtk, "10\n");  // 10 : tetrahedron
      }
    }
    else if(!strcmp(pmesh->type, "fiber_real")){
      fprintf(fpvtk, "CELLS %d %d\n", pmesh->all_elements, pmesh->all_elements*3);
      for(i=0; i<pmesh->all_elements; i++){
        fprintf(fpvtk, "2 %d %d\n", pmesh->EC[i][0], pmesh->EC[i][1]);
      }
      fprintf(fpvtk, "CELL_TYPES %d\n", pmesh->all_elements);
      for(i=0; i<pmesh->all_elements; i++)
        fprintf(fpvtk, "3\n");  // 3: fiber element
    }
    fprintf(fpvtk, "POINT_DATA %d\n", pmesh->all_nodes);
#if defined VTK_COMBINED
   for(int i =1; i<vtkout_states_vec.size(); i++){
     fprintf(fpvtk, "SCALARS %s float\n", vtkout_states_vec[i].c_str());
     fprintf(fpvtk, "LOOKUP_TABLE default\n");
     if(strcmp(vtkout_states_vec[i].c_str(), "V") == 0){
       mul = 10.0;
     }
     else{
       mul = 10000000.0;
     }
     for(j=0; j<pmesh->all_nodes; j++)
      //fprintf(fpvtk, "%d\n", (int)(Voltsc[i][j]*mul));  // Vm
      fprintf(fpvtk, "%lf\n", (Voltsc[i][j]));  // Vm
   }
   for(int i =0; i<vtkout_algebraic_vec.size(); i++){
     fprintf(fpvtk, "SCALARS %s float\n", vtkout_algebraic_vec[i].c_str());
     fprintf(fpvtk, "LOOKUP_TABLE default\n");
     for(j=0; j<pmesh->all_nodes; j++)
       //fprintf(fpvtk, "%d\n", (int)(Voltsc1[i][j]*1000.0));  // Vm
       fprintf(fpvtk, "%lf\n", (Voltsc1[i][j]));  // Vm
   }
   fclose(fpvtk);
   PetscPrintf(PETSC_COMM_WORLD, "%s has been successfully generated!!!\n", fname);
#endif

  if(!strcmp(pm.mesh_type, "tetfiber")) pfiber->writePF(time);
}


void Cheart::writeScalar()
{
  char fname[200];
  int i,j;
  double mul;
  double stepper;

  stepper = dt >= 0.5 ? 0. : 0.5;
  sprintf(fname, "V_%d.scalar", (int)floor(time+stepper));
  FILE* fp_scalar=fopen(fname,"wt");
  for(int i =0; i<vtkout_states_vec.size(); i++){
    fprintf(fp_scalar, "SCALARS %s float\n", vtkout_states_vec[i].c_str());
    fprintf(fp_scalar, "LOOKUP_TABLE default\n");
    if(strcmp(vtkout_states_vec[i].c_str(), "V") == 0){
      mul = 10.0;
    }
    else{
      mul = 10000000.0;
    }
    for(j=0; j<pmesh->all_nodes; j++)
      fprintf(fp_scalar, "%d\n", (int)(Voltsc[i][j]*mul));  // Vm
  }
  for(int i =0; i<vtkout_algebraic_vec.size(); i++){
    fprintf(fp_scalar, "SCALARS %s float\n", vtkout_algebraic_vec[i].c_str());
    fprintf(fp_scalar, "LOOKUP_TABLE default\n");
    for(j=0; j<pmesh->all_nodes; j++)
      fprintf(fp_scalar, "%d\n", (int)(Voltsc1[i][j]*1000.0));  // Vm
  }
  fclose(fp_scalar);
  PetscPrintf(PETSC_COMM_WORLD, "%s has been successfully generated!!!\n", fname);
  
}

void Cheart::writeoutput()
{
  char fname[200];
  int i, j;
  FILE* fpoutput;
  if((iter_count != 0)) {
    sprintf(fname,"output.%d",mympi::rank);
    fpoutput = fopen(fname,"w");
    PetscFPrintf(PETSC_COMM_SELF,fpoutput,"%f \n",time);
    for(i=0; i<pmesh->all_nodes; i++){
      if(node_proc[i]==mympi::rank){
        ilocal=local_nodes[i];
	for(j=0; j < p_cell->states_size ; j++){
          PetscFPrintf(PETSC_COMM_SELF,fpoutput,"%f ",Glob_cell[j][ilocal]);
        }
        PetscFPrintf(PETSC_COMM_SELF,fpoutput,"\n");
      }
    }
    fclose(fpoutput);
  }
}

// assign S1 and S2
bool Cheart::loadPacingsites(char* pace1, char* pace2)
{
  FILE* fp;
  char fname[1000];
  int i, s2nd; 
  // load S1 
  sprintf(fname, "%s", pace1);
  if(!(fp = fopen(fname, "rt"))) {PetscPrintf(PETSC_COMM_WORLD, " Error to open %s\n", fname); exit(-1);}
  fscanf(fp, "%d", &num_s1nodes);
  PetscPrintf(PETSC_COMM_WORLD, "num_s1nodes: %d\n", num_s1nodes);
  s1_num = new int [num_s1nodes];
  for(i=0; i<num_s1nodes; i++){
    fscanf(fp, "%d", &s1_num[i]);
    s1_num[i]--; // node number 1 in inp is node number 0 in code !
    s1nodes[s1_num[i]] = true;
  }
  fclose(fp);

  if(pm.is_s1s2){ // load S2
    sprintf(fname, "%s", pace2);
    if(!(fp = fopen(fname, "rt"))) {PetscPrintf(PETSC_COMM_WORLD, "Error to open %s\n", fname); exit(-1);}
    fscanf(fp, "%d", &num_s2nodes);
    PetscPrintf(PETSC_COMM_WORLD, "num_s2nodes: %d\n", num_s2nodes);
    s2_num = new int [num_s2nodes];
    for(i=0; i<num_s2nodes; i++){
      fscanf(fp, "%d", &s2_num[i]);
      s2_num[i]--;
      if(pm.is_s1s2==1) s2nodes[s2_num[i]] = true;  // for reset S2
    }
    fclose(fp);
  }
  return true;
}


bool Cheart::FindSynapses()
{
  int i, j, ii, iel, idistmin;
  double x1, y1, z1, xx, yy, zz, distmin, dist;
  pmesh_endosurf = new Cmesh(pm.surface_mesh);
  pfiber->pmesh->inearsurf = new int [pfiber->pmesh->all_terminals];
  //===Find Tissue node nearest to Nerve perpheral node->inearsurf(itot)===
  for(i=0; i<pfiber->pmesh->all_terminals; i++){  // 모든 말단노드에 대하여
    ii = pfiber->pmesh->iPFnode[i]; // i번째 말단노드의 전역노드번호
    x1 = pfiber->pmesh->XY[ii][x]; // x1, y1, z1 : 말단노드의 좌표값
    y1 = pfiber->pmesh->XY[ii][y];
    z1 = pfiber->pmesh->XY[ii][z];
    //PetscPrintf(PETSC_COMM_WORLD, "x1: %lf   y1: %lf   z1: %lf\n", x1, y1, z1);
    distmin = 10000.0;
    for(iel=0; iel<pmesh_endosurf->all_elements; iel++){ // 모든 표면 엘리먼트에 대하여
      for(j=0; j<pmesh_endosurf->node; j++){   // 각 엘리먼트의 모든 노드에 대하여
        ii = pmesh_endosurf->EC[iel][j]; // iel요소의 i번째 노드의 전역노드번호
        xx = pmesh->XY[ii][x];
        yy = pmesh->XY[ii][y];
        zz = pmesh->XY[ii][z];
        dist = sqrt((xx-x1)*(xx-x1) + (yy-y1)*(yy-y1) + (zz-z1)*(zz-z1));
        //printf("xx: %lf   yy: %lf   zz: %lf\n", xx, yy, zz);

        if(dist<distmin){
          distmin=dist;  // 가장가까운 거리
          idistmin=ii;   // 가장가까운 거리에 있는 노드번호
        }
      }
    }
    pfiber->pmesh->inearsurf[i]=idistmin;  // 가장가까운 거리에 있는 노드번호저장
    PetscPrintf(PETSC_COMM_WORLD, "iPFnode[%d] = %d\n", i, pfiber->pmesh->iPFnode[i]);
    PetscPrintf(PETSC_COMM_WORLD, "inearsurf[%d] = %d\n", i, pfiber->pmesh->inearsurf[i]);
  }

  return true;
}

void Cheart::GetTargetnode()
{
  pmesh->target[0] = pfiber->pmesh->inearsurf[0];
  pmesh->target[1] = pfiber->pmesh->inearsurf[1];
  pmesh->target[2] = pfiber->pmesh->inearsurf[2];
}

int Cheart::ifInList(int array[], int size, int n) {
    int j;
    for(j=0; j<size; j++)
        if(array[j]==n)
            return 1;
    return 0;
}

drug_t Cheart::get_IC50_data_from_file(const char* file_name)
{
  FILE *fp_drugs;
  drug_t ic50;
  char *token, buffer[255];
  row_data temp_array;
  unsigned short idx;

  if( (fp_drugs = fopen(file_name, "r")) == NULL){
    printf("Cannot open file %s in %s at rank %d\n",
      file_name, mympi::host_name, mympi::rank);
    return ic50;
  }

  fgets(buffer, sizeof(buffer), fp_drugs); // skip header
  while( fgets(buffer, sizeof(buffer), fp_drugs) != NULL )
  { // begin line reading
    token = strtok( buffer, "," );
    idx = 0;
    while( token != NULL )
    { // begin data tokenizing
      temp_array.data[idx++] = strtod(token, NULL);
      token = strtok(NULL, ",");
    } // end data tokenizing
    ic50.push_back(temp_array);
  } // end line reading

  fclose(fp_drugs);
  return ic50;
}

