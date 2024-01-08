#include "cipa_drug_herg_bench.hpp"
#include "commons.hpp"
#include "../cellmodels/ohara_rudy_cipa_v1_2017.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <cvode/cvode.h>
#include <cvode/cvode_dense.h>
#include <mpi.h>
#include <nvector/nvector_serial.h>
#include <string>

#if defined( __linux__ )
  #include <sys/types.h>
  #include <sys/stat.h>
  #include <unistd.h>
#elif defined( _WIN32 )
  #include <direct.h>
#endif

void cipa_drug_herg_bench( int argc, char **argv, param_t *p_param )
{
  // CellML object pointer
  patch_clamp *p_cell;

  // Supporting variables
  char buffer[150] = {'\0'};
  char *pch;
  double **hills;
  double *concs;
  FILE *fp_hill;
  FILE *fp_herg;
  int idx, group_id, sim_counter, iout, imax;
  int hill_size, herg_size, sample_size, concs_size;
  int print_freq;
  struct stat st = {0};

  // CiPA related variables
  cipa_t *cipa_result, *temp_result;
  double INaL_auc_control, ICaL_auc_control, INaL_auc_drug, ICaL_auc_drug, steepest_time;
  FILE *fp_vm = NULL, *fp_ca = NULL, *fp_ires = NULL;
  FILE *fp_qnet, *fp_qinward, *fp_ap_profile, *fp_ca_profile, *fp_dvmdt_250;
  int pace, steepest_pace;

  // Calcium related variable
  bool is_ca_amp90_inc, is_ca_amp50_inc, is_cad90_found, is_cad50_found, is_ca_peak_found;
  double t1_cad90, t1_cad50, ca_amp50, ca_amp90, ca_prev;

  // AP related variables
  bool is_eligible_AP;
  double t_depol, vm_repol30, vm_repol90, vm_repol50;

  // SUNDIALs variables
  int cvode_retval;
  double tnext, tcurr = 0.0;
  void *cvode_mem;
  N_Vector states_vec;

  // give limitations to the number of pacing
  if( p_param->num_pace1 > 1000 ){
    printf("Number of pacing cannot exceed 1000! Exiting...\n");
    return;
  }
  else if( p_param->num_pace1 < 250 ){
    printf("Number of pacing should be more than 250! Exiting...\n");
    return;
  }

  // get total of concentrations and 
  // allocate array with that value while giving some limitations
  std::string concs_str(p_param->concs);
  concs_size = std::count( concs_str.begin(), concs_str.end(), ',' ) + 1;
  if( concs_size > 4+1 ){  // +1 because we include concentration 0
    fprintf(stderr, "Only 4 concentrations is allowed in this program! Exiting...\n");
    return;
  }
  concs = new double[concs_size];

  // Split concentration strings into array
  for( idx = 0, pch = strtok( p_param->concs, ", ");  idx < concs_size; idx++, pch = strtok(NULL, ",") ) concs[idx] = strtod( pch, NULL );

  // open the hill sample file
  fp_hill = fopen(p_param->hill_file, "r");
  if(fp_hill == NULL) {
    fprintf(stderr, "Cannot open Hill file %s\n", p_param->hill_file);
    return;
  }
  // get the number of hill sample and give some limitations
  fgets(buffer, sizeof(buffer), fp_hill);
  hill_size = 0;
  while (fgets(buffer, sizeof(buffer), fp_hill) != NULL) hill_size++;

  // open the herg sample file
  fp_herg = fopen(p_param->herg_file, "r");
  if(fp_herg == NULL) {
    fprintf(stderr, "Cannot open hERG file %s\n", p_param->herg_file);
    return;
  }
  // get the number of herg sample and check if hill sample number is equal to herg sample number
  fgets(buffer, sizeof(buffer), fp_herg);
  herg_size = 0;
  while (fgets(buffer, sizeof(buffer), fp_herg) != NULL) herg_size++;

  sample_size = 0;
  if( hill_size == herg_size ){
    sample_size = hill_size;
  }
  else if( hill_size != herg_size){
    fprintf(stderr, "The number of hill sample and hERG sample are not same! Hill: %d, hERG: %d\n", hill_size, herg_size);
    return;
  }

  // allocate hills array
  hills = new double*[sample_size];
  for( idx = 0; idx < sample_size; idx++ ) hills[idx] = new double[20];

  // Reset the FILE pointer to the beginning
  if( fseek(fp_hill, 0L, SEEK_SET) != 0 ){
    fprintf(stderr, "Reset file pointer of Hill file got error\n");
    return;
  }
  if( fseek(fp_herg, 0L, SEEK_SET) != 0 ){
    fprintf(stderr, "Reset file pointer of hERG file got error\n");
    return;
  }

  // copy sample to array
  fgets(buffer, sizeof(buffer), fp_hill);
  group_id = 0;
  while (fgets(buffer, sizeof(buffer), fp_hill) != NULL) {
      for (idx = 0, pch = strtok(buffer, ","); pch != NULL; idx++, pch = strtok(NULL, ",")) {
          hills[group_id][idx] = strtod(pch, NULL);
      }
      group_id++;
  }
  fclose(fp_hill);
  fgets(buffer, sizeof(buffer), fp_herg);
  group_id = 0;
  while (fgets(buffer, sizeof(buffer), fp_herg) != NULL) {
      for (idx = 14, pch = strtok(buffer, ","); pch != NULL; idx++, pch = strtok(NULL, ",")) {
          hills[group_id][idx] = strtod(pch, NULL);
      }
      group_id++;
  }
  fclose(fp_herg);

  p_cell = new ohara_rudy_cipa_v1_2017();

  // begin sample loop
  double begin = MPI_Wtime();
  int sample_group = ( sample_size%mympi::size != 0 ) ? (sample_size/mympi::size)+1 : (sample_size/mympi::size);
  int sample_id;
  for( group_id=0, sim_counter=mympi::rank; 
       group_id < sample_group && sim_counter < sample_size; 
       group_id++, sim_counter += mympi::size)
  { 
    sample_id = (group_id*mympi::size)+mympi::rank;
    printf("Sample_ID:%d  Count:%d Sample Group:%d Rank:%d Group_ID:%d\nHill Data: ", sample_id, sample_size ,sample_group, mympi::rank, group_id );
    for( idx = 0; idx < 14; idx++ ){
      printf("%lf ", hills[sample_id][idx]);
    }
    printf("\nhErg Data: ");
    for( idx = 14; idx < 20; idx++ ){
      printf("%lf ", hills[sample_id][idx]);
    }
    printf("\n");

    for( idx = 0; idx < concs_size; idx++ )
    { // begin concs loop

      // apply some cell initialization
      p_cell->initConsts(p_param->celltype, concs[idx], hills[sample_id]);
      p_cell->CONSTANTS[i_Stim_Period] = p_param->bcl_init;
      p_cell->CONSTANTS[i_Stim_Amplitude] *= p_param->stim_amp;

      // Time and CVode setup
      tnext = p_param->dt;
      tcurr = 0.;
      cvode_mem = CVodeCreate( CV_BDF,CV_NEWTON );
      CVodeSetUserData( cvode_mem, p_cell );
      states_vec = N_VMake_Serial( p_cell->states_size, p_cell->STATES );
      CVodeInit( cvode_mem, rhs_fn, 0.0, states_vec );
      CVodeSetMaxStep( cvode_mem, p_cell->CONSTANTS[i_Stim_PulseDuration] );
      CVDense( cvode_mem, p_cell->states_size );
      CVodeSStolerances( cvode_mem, 1.0e-7, 1.0e-7 );
      print_freq = (1. / p_param->dt) * p_param->dt_write;

      snprintf( buffer, sizeof(buffer), "result/%.2lf", concs[idx] );
#if defined( __linux__ )
      if (stat("result/", &st) == 0) mkdir(buffer, 0775);
#elif defined( _WIN32 )
      _mkdir(buffer);
#endif
      // FILE pointer naming and initialization
      if( p_param->is_print_graph == 1 ){
        snprintf(buffer, sizeof(buffer), "./result/%.2lf/%s_%.2lf_vmcheck_%d.plt", concs[idx], p_param->drug_name, concs[idx], mympi::rank );
        fp_vm = fopen( buffer, "a" );
        //if (fp_vm == NULL) printf("vmcheck at sample %d cannot be created!!\n", sample_id); return;
        snprintf(buffer, sizeof(buffer), "./result/%.2lf/%s_%.2lf_ca_i_%d.plt", concs[idx], p_param->drug_name, concs[idx], mympi::rank );
        fp_ca = fopen( buffer, "a" );
       // if (fp_ca == NULL) printf("ca_i at sample %d cannot be created!!\n", sample_id); return;
        snprintf(buffer, sizeof(buffer), "./result/%.2lf/%s_%.2lf_ires_%d.plt", concs[idx], p_param->drug_name, concs[idx], mympi::rank );
        fp_ires = fopen( buffer, "a" );
        //if (fp_ires == NULL) printf("ires at sample %d cannot be created!!\n", sample_id); return;
      }
      snprintf(buffer, sizeof(buffer), "result/%.2lf/%s_%.2lf_ap_profile_%d.plt", concs[idx], p_param->drug_name, concs[idx], mympi::rank );
      fp_ap_profile = fopen( buffer, "a" );
      snprintf(buffer, sizeof(buffer), "result/%.2lf/%s_%.2lf_ca_profile_%d.plt", concs[idx], p_param->drug_name, concs[idx], mympi::rank );
      fp_ca_profile = fopen( buffer, "a" );
      snprintf(buffer, sizeof(buffer), "result/%.2lf/%s_%.2lf_qnet_%d.plt", concs[idx], p_param->drug_name, concs[idx], mympi::rank );
      fp_qnet = fopen( buffer, "a" );
      snprintf(buffer, sizeof(buffer), "result/%.2lf/%s_%.2lf_qinward_%d.plt", concs[idx], p_param->drug_name, concs[idx], mympi::rank );
      fp_qinward = fopen( buffer, "a" );
      snprintf(buffer, sizeof(buffer), "result/%.2lf/%s_%.2lf_dvmdt_250_%d.plt", concs[idx], p_param->drug_name, concs[idx], mympi::rank );
      fp_dvmdt_250 = fopen( buffer, "a" );

      // Write header for all of the result file
      if( group_id == 0 )
      {
        if( p_param->is_print_graph == 1 ){
          fprintf( fp_vm, "%s %s\n", "TIME", "Vm" );
          fprintf( fp_ca, "%s %s\n", "TIME", "cai" );
          fprintf( fp_ires, "%s %s %s %s %s %s %s %s %s %s\n",
                   "TIME","INa","INaL","ICaL","Ito",
                   "IKr","IKs","IK1","INaCa","INaK" );
        }
        fprintf( fp_ap_profile, "%s %s %s %s %s %s %s %s\n",
                 "Sample_ID", "Max_Dvm/Dt", "Vm_Peak", "Vm_Resting","APD90", "APD50", "APDTri", "Steepest_Pace");
        fprintf( fp_ca_profile, "%s %s %s %s %s %s\n",
                 "Sample_ID", "Ca_Peak", "Ca_Diastole", "CaD90", "CaD50","Catri");
        fprintf( fp_qnet, "%s %s\n","Sample_ID", "Qnet");
        fprintf( fp_qinward, "%s %s\n", "Sample_ID", "QInward");
      }

      // initialize some variables
      temp_result = new cipa_t;
      cipa_result = new cipa_t;
      cipa_result->dvmdt_repol = -1.;
      ca_prev = 0.;
      is_ca_peak_found = false;
      iout = 0;
      imax = ( p_param->num_pace1 *  p_cell->CONSTANTS[i_Stim_Period] ) / p_param->dt;
      pace = steepest_pace = 0;

      // begin simulation loop
      while(iout < imax)
      {
        // solve ODE
        cvode_retval = CVode( cvode_mem, tnext, states_vec, &tcurr, CV_NORMAL  );
        // get the RATES array
        p_cell->computeRates(tnext, p_cell->CONSTANTS, p_cell->RATES, N_VGetArrayPointer_Serial(states_vec), p_cell->ALGEBRAIC);
        if( cvode_retval == CV_SUCCESS ){
          iout++;
          tnext += p_param->dt;
        }
        else{
          fprintf(stderr, "CVode error at sample_ID %d and concentration %.2lf at rank %d\n", sample_id, concs[idx], mympi::rank);
          break;
        }

        // executed at each period of BCL
        if( iout % (int)( p_cell->CONSTANTS[i_Stim_Period] / p_param->dt ) == 0 ) {
          // executed at the end of AP pacing
          if( is_eligible_AP && pace >= p_param->num_pace1-249) {
            temp_result->qnet_curr = p_cell->STATES[qnet];
            temp_result->INaL_auc_prev = p_cell->STATES[INaL_AUC];
            temp_result->ICaL_auc_prev = p_cell->STATES[ICaL_AUC];
            temp_result->vm_dia = p_cell->STATES[V];
            temp_result->ca_dia = p_cell->STATES[cai];
            fprintf( fp_dvmdt_250, "%d %.15f %.2f %.2f %.2f\n", pace, temp_result->dvmdt_repol, steepest_time, vm_repol30, vm_repol90 );
            // replace result with steeper repolarization AP
            if( temp_result->dvmdt_repol > cipa_result->dvmdt_repol ) {
              if(mympi::rank == 0 && steepest_pace != 0)printf("Pace %d dvmdt_repol %.15lf replaced by pace %d dvmdt_repol %.15lf\n", steepest_pace, cipa_result->dvmdt_repol,pace, temp_result->dvmdt_repol);
              steepest_pace = pace;
              memcpy(cipa_result, temp_result, sizeof(cipa_t));
              //if(mympi::rank == 0)printf("Pace %d replaced by pace %d\n", pace-1,pace);
              //if(mympi::rank == 0)printf("Pace %d dvmdt_repol %lf replaced by pace %d dvmdt_repol %lf\n", pace-1, cipa_result->dvmdt_repol,pace, temp_result->dvmdt_repol);
              //if(mympi::rank == 0)printf("Pace %d: Qnet_curr=%lf Qnet_prev=%lf Delta=%lf DVmDT_Repol=%lf\n", pace, cipa_result->qnet_curr, cipa_result->qnet_prev, cipa_result->qnet_curr - cipa_result->qnet_prev, cipa_result->dvmdt_repol);
              //if(mympi::rank == 0)printf("Pace %d: Vm_Dia: %lf Ca_Dia: %lf\n", pace, cipa_result->vm_dia, cipa_result->ca_dia);
            }
          }

          pace++;

          // executed at the beginning of AP pacing
          if( pace >= p_param->num_pace1-250 ) {
            //if(mympi::rank == 0)printf("Time new cycle: %lf with pace: %d\n", tcurr, pace);
            temp_result->qnet_prev = p_cell->STATES[qnet];
            temp_result->INaL_auc_prev = p_cell->STATES[INaL_AUC];
            temp_result->ICaL_auc_prev = p_cell->STATES[ICaL_AUC];
            temp_result->vm_valley = p_cell->STATES[V];
            temp_result->dvmdt_max = p_cell->STATES[V];
            temp_result->dvmdt_repol = -1.;
            temp_result->ca_valley = p_cell->STATES[cai];
            is_eligible_AP = false;
            is_ca_peak_found = false;
            is_ca_amp90_inc = true;
            is_ca_amp50_inc = true;
          }
        }

        if( pace == p_param->num_pace1-251 ){
          if( ca_prev > p_cell->STATES[cai] && !is_ca_peak_found ){
            temp_result->ca_peak = ca_prev;
            ca_amp50 = temp_result->ca_peak - (0.5 * (temp_result->ca_peak - temp_result->ca_valley));
            ca_amp90 = temp_result->ca_peak - (0.9 * (temp_result->ca_peak - temp_result->ca_valley));
            is_ca_peak_found = true;
            is_cad50_found = false;
            is_cad90_found = false;
            //if(group_id == 0) printf("Ca_Peak: %lf Ca50: %lf Ca90: %lf\n", temp_result->ca_peak, ca_amp50, ca_amp90 );
          }
          if( !is_ca_amp50_inc && is_ca_peak_found && ca_amp50 > p_cell->STATES[cai] && !is_cad50_found ){
            t1_cad50 = tcurr;
            is_ca_amp50_inc = true;
            //if(mympi::rank == 0) printf("T2_CAD50: %lf T1_CAD50: %lf CAD50: %lf\n", tcurr, t1_cad50, temp_result->cad50 );
          }
          if( !is_ca_amp90_inc && is_ca_peak_found && ca_amp90 > p_cell->STATES[cai] && !is_cad90_found ){
            t1_cad90 = tcurr;
            is_ca_amp90_inc = true;
            //if(mympi::rank == 0) printf("T2_CAD90: %lf T1_CAD90: %lf CAD90: %lf\n", tcurr, t1_cad90, temp_result->cad90 );
          }
          if( !is_ca_peak_found ) ca_prev = p_cell->STATES[cai];
        }

        // entering the last 250 paces
        if( p_param->num_pace1 > 250 && p_param->num_pace1-250 <= pace && pace < p_param->num_pace1 ){
          // get maximum dvmdt
          if( temp_result->dvmdt_max < p_cell->RATES[V] ) temp_result->dvmdt_max = p_cell->RATES[V];
          // get the peak Vm 6 secs after depolarization (when Na channel just closed after bursting)
          if( iout % (int)(((p_cell->CONSTANTS[i_Stim_Period]*pace)+(p_cell->CONSTANTS[i_Stim_Start]+6.)) / p_param->dt) == 0) {
            temp_result->vm_peak = p_cell->STATES[V];
            t_depol = tcurr;
            // only finding the vm_repol30, vm_repol_50 and vm_repol90 if the peak Vm more than 0 mV
            if( temp_result->vm_peak > 0. ){
              vm_repol30 = temp_result->vm_peak - (0.3 * (temp_result->vm_peak - temp_result->vm_valley));
              vm_repol50 = temp_result->vm_peak - (0.5 * (temp_result->vm_peak - temp_result->vm_valley));
              vm_repol90 = temp_result->vm_peak - (0.9 * (temp_result->vm_peak - temp_result->vm_valley));
              is_eligible_AP = true;
              //if(mympi::rank == 0) printf("Peak Vm and vm_valley and Vm Repol 30%% and 90%% at pace %d concentration %lf: %lf %lf %lf %lf\n", pace, concs[idx], p_cell->STATES[V], temp_result->vm_valley, vm_repol30, vm_repol90);
            }
          }
          // these operations will be executed if it's eligible AP and done after the beginning of repolarization
          else if( is_eligible_AP && tcurr > (p_cell->CONSTANTS[i_Stim_Period]*pace)+(p_cell->CONSTANTS[i_Stim_Start]+6.) ){
            // get the steepest dvmdt_repol around vm_repol30 and vm_repol90 if AP shape is eligible
            if( vm_repol30 >= p_cell->STATES[V] && p_cell->STATES[V] >= vm_repol90
                && temp_result->dvmdt_repol < p_cell->RATES[V] ){
                temp_result->dvmdt_repol = p_cell->RATES[V];
                steepest_time = tcurr;
            }
            // get the APD50
            if( vm_repol50 > p_cell->STATES[V] && p_cell->STATES[V] > vm_repol50-0.1 ){
              temp_result->apd50 = tcurr - t_depol;
              //if(mympi::rank == 0) printf("REPOL50: %lf DEPOL: %lf APD50: %lf vm_repol50=%lf curr_v=%lf\n", tcurr, t_depol, temp_result->apd50, vm_repol50, p_cell->STATES[V] );
            }
            // get the APD90
            if( vm_repol90 > p_cell->STATES[V] && p_cell->STATES[V] > vm_repol90-0.2 ){
              temp_result->apd90 = tcurr - t_depol;
              //if(mympi::rank == 0) printf("REPOL90: %lf DEPOL: %lf APD90: %lf vm_repol90=%lf curr_v=%lf\n", tcurr, t_depol, temp_result->apd90, vm_repol90, p_cell->STATES[V] );
            }
            // get the peak calcium and amplitude 50% and 90%
            if( ca_prev > p_cell->STATES[cai] && !is_ca_peak_found ){
              temp_result->ca_peak = ca_prev;
              ca_amp50 = temp_result->ca_peak - (0.5 * (temp_result->ca_peak - temp_result->ca_valley));
              ca_amp90 = temp_result->ca_peak - (0.9 * (temp_result->ca_peak - temp_result->ca_valley));
              is_ca_peak_found = true;
              is_cad90_found = false;
              is_cad50_found = false;
              //if(group_id == 0) printf("Ca_Peak: %lf Ca50: %lf Ca90: %lf\n", temp_result->ca_peak, ca_amp50, ca_amp90 );
            } 

            if( !is_ca_amp50_inc && is_ca_peak_found && ca_amp50 > p_cell->STATES[cai] ){
              t1_cad50 = tcurr;
              is_ca_amp50_inc = true;
            }
            if( !is_ca_amp90_inc && is_ca_peak_found && ca_amp90 > p_cell->STATES[cai] ){
              t1_cad90 = tcurr;
              is_ca_amp90_inc = true;
            }
          }

          if( is_ca_amp50_inc && is_ca_peak_found && ca_amp50 < p_cell->STATES[cai] ){
            temp_result->cad50 = p_cell->CONSTANTS[i_Stim_Period] - (tcurr - t1_cad50);
            is_ca_amp50_inc = false;
            is_cad50_found = true;
            //if(mympi::rank == 0) printf("T2_CAD50: %lf T1_CAD50: %lf CAD50: %lf pace: %d concs: %lf\n", tcurr, t1_cad50, temp_result->cad50, pace, concs[idx] );
          }
          if( is_ca_amp90_inc && is_ca_peak_found && ca_amp90 < p_cell->STATES[cai] ){
            temp_result->cad90 =  p_cell->CONSTANTS[i_Stim_Period] - (tcurr - t1_cad90);
            is_ca_amp90_inc = false;
            is_cad90_found = true;
            //if(mympi::rank == 0) printf("T2_CAD90: %lf T1_CAD90: %lf CAD90: %lf pace: %d concs: %lf\n", tcurr, t1_cad90, temp_result->cad90, pace, concs[idx] );
          }

          if( !is_ca_peak_found ) ca_prev = p_cell->STATES[cai];

          if( tcurr > p_cell->CONSTANTS[i_Stim_Period] * (p_param->num_pace1-1) ){
            // print graph
            if(iout % print_freq == 0 && p_param->is_print_graph == 1 ){
              fprintf( fp_vm, "%lf %lf\n", tcurr, p_cell->STATES[V] );
              fprintf( fp_ca, "%lf %lf\n", tcurr, p_cell->STATES[cai] );
              fprintf( fp_ires, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
                       tcurr, p_cell->ALGEBRAIC[INa], p_cell->ALGEBRAIC[INaL], p_cell->ALGEBRAIC[ICaL],
                       p_cell->ALGEBRAIC[Ito], p_cell->ALGEBRAIC[IKr], p_cell->ALGEBRAIC[IKs],
                       p_cell->ALGEBRAIC[IK1], p_cell->ALGEBRAIC[INaCa_i], p_cell->ALGEBRAIC[INaK] );
            }
          }
        }// end of last 250 paces operations


        // last cycle
/*
        if( tcurr > p_cell->CONSTANTS[i_Stim_Period] * (p_param->num_pace1-1) ){
          // print graph
          if(iout % print_freq == 0 && p_param->is_print_graph == 1 ){
            fprintf( fp_vm, "%lf %lf\n", tcurr, p_cell->STATES[V] );
            fprintf( fp_ca, "%lf %lf\n", tcurr, p_cell->STATES[cai] );
            fprintf( fp_ires, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
                     tcurr, p_cell->ALGEBRAIC[INa], p_cell->ALGEBRAIC[INaL], p_cell->ALGEBRAIC[ICaL],
                     p_cell->ALGEBRAIC[Ito], p_cell->ALGEBRAIC[IKr], p_cell->ALGEBRAIC[IKs],
                     p_cell->ALGEBRAIC[IK1], p_cell->ALGEBRAIC[INaCa_i], p_cell->ALGEBRAIC[INaK] );
          }
        } // end of last cycle
*/
      } // end simulation loop

      // print last result
      if( cvode_retval == CV_SUCCESS ){
        fprintf( fp_ap_profile, "%d %lf %lf %lf %lf %lf %lf %d\n",
               sample_id,
               cipa_result->dvmdt_max,
               cipa_result->vm_peak,
               cipa_result->vm_dia,
               cipa_result->apd90,
               cipa_result->apd50,
               cipa_result->apd90-cipa_result->apd50,
               steepest_pace);
        fprintf( fp_ca_profile, "%d %lf %lf %lf %lf %lf\n",
               sample_id,
               cipa_result->ca_peak,
               cipa_result->ca_dia,
               cipa_result->cad90,
               cipa_result->cad50,
               cipa_result->cad90-cipa_result->cad50 );
        fprintf( fp_qnet, "%d %lf\n", sample_id, (cipa_result->qnet_curr-cipa_result->qnet_prev)/1000.0 );
        if((int)ceil(concs[idx]) == 0){
          INaL_auc_control = cipa_result->INaL_auc_curr - cipa_result->INaL_auc_prev;
          ICaL_auc_control = cipa_result->ICaL_auc_curr - cipa_result->ICaL_auc_prev;
        }
        else{
          INaL_auc_drug = cipa_result->INaL_auc_curr - cipa_result->INaL_auc_prev;
          ICaL_auc_drug = cipa_result->ICaL_auc_curr - cipa_result->ICaL_auc_prev;
          fprintf( fp_qinward, "%d %lf\n", sample_id, ( (INaL_auc_drug/INaL_auc_control) + (ICaL_auc_drug/ICaL_auc_control) ) * 0.5 );
        }
      }
      else{
        fprintf( fp_ap_profile, "%d ERR ERR ERR ERR ERR ERR\n", sample_id);
        fprintf( fp_ca_profile, "%d ERR ERR ERR ERR ERR\n", sample_id);
        fprintf( fp_qnet, "%d ERR\n", sample_id);
        fprintf( fp_qinward, "%d ERR\n", sample_id);
      }
      printf( "Steepest AP pace for sample_ID %d at concentration %.2lf: %d\n", sample_id, concs[idx], steepest_pace );

      if( p_param->is_print_graph == 1 ){
        fprintf( fp_vm, "---------------------------------------end of sample %d\n", sample_id );
        fprintf( fp_ca, "---------------------------------------end of sample %d\n", sample_id );
        fprintf( fp_ires, "---------------------------------------end of sample %d\n", sample_id );
      }
     fprintf( fp_dvmdt_250, "---------------------------------------end of sample %d\n", sample_id );

      if( p_param->is_print_graph == 1 ){
        fclose( fp_vm );
        fclose( fp_ca );
        fclose( fp_ires );
      }
      fclose( fp_ap_profile );
      fclose( fp_ca_profile );
      fclose( fp_qinward );
      fclose( fp_qnet );
      fclose( fp_dvmdt_250 );

      N_VDestroy(states_vec);
      CVodeFree(&cvode_mem);
      delete temp_result;
      delete cipa_result;

    }  //end concs loop
  } // end sample loop
  double end = MPI_Wtime();

  // master rank will make a performance log and clean all used memory
  if( mympi::rank == 0 ){
    char buff_time[15] = {'\0'};
    FILE *fp_perf_log;
    struct tm* tm_info;
    time_t timer;
    
    timer = time(NULL);
    tm_info = localtime(&timer);
    strftime( buff_time, 15, "%Y%m%d%H%M%S", tm_info );
    snprintf( buffer, sizeof(buffer), "result/performance_%s.log", buff_time );
    fp_perf_log = fopen( buffer, "w" );
    fprintf( fp_perf_log, "Iteration CPU Sample Pacing Time(Minutes)\n" );
    fprintf( fp_perf_log, "%lld %d %d %d %lf\n", (int)(sample_size*concs_size*p_param->num_pace1*p_cell->CONSTANTS[i_Stim_Period]/p_param->dt), mympi::size, sample_size, p_param->num_pace1, (end-begin)/60. );
    printf("Computation Time: %lf minutes\n", (end-begin)/60.);

    fclose( fp_perf_log );
  }
  for( idx = 0; idx < sample_size; idx++ ) delete[] hills[idx];
  delete[] hills;
  delete[] concs;
  delete p_cell;

}
