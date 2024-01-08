#include "patches/Ohara_Rudy_2011.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <cvode/cvode.h>
#include <cvode/cvode_dense.h>
#include <nvector/nvector_serial.h>

double glob_vm = 0.0;


int rhs_fn( realtype t, N_Vector y, N_Vector ydot, void *user_data );

int rhs_fn( realtype t, N_Vector y, N_Vector ydot, void *user_data )
{
  NV_Ith_S(y, V) = glob_vm;
  patch_clamp *data = (patch_clamp*)user_data;
  data->computeRates( t,
                      data->CONSTANTS,
                      N_VGetArrayPointer_Serial(ydot),
                      N_VGetArrayPointer_Serial(y),
                      data->ALGEBRAIC );
  return 0;
}


int main( int argc, char **argv )
{
  if(argc < 2){
    printf("Please provide the membrane potential input file!!\n");
    return 1;
  }

  // I/O variables
  char buff[255];
  char hill_sample[255];
  char *pch;
  FILE *fp_curr;
  FILE *fp_hill;
  FILE *fp_vm;

  // Patch clamp variables
  patch_clamp *p_patch;
  double bcl,dt;

  // Drug related variables
  bool is_drug;
  double hill[14];
  double conc;

  // SUNDIALs variables
  int retval;
  double tnext, tcurr;
  void *cvode_main;
  N_Vector states_vec;

  // extra variables
  int idx;


  is_drug = false;
  if(argc >= 3){
    fp_hill = fopen(argv[2], "r");
    if(fp_hill != NULL) {is_drug = true;}
    else
    {  
      printf("File %s not found!! Applying the normal Patch Clamp!\n", argv[2]);
    }
  }

  if(is_drug == true){
    /* Skip the header in the sample file */
    fgets(hill_sample, sizeof(hill_sample), fp_hill);
    fgets(hill_sample, sizeof(hill_sample), fp_hill);
    pch = strtok(hill_sample, ",");
    if(pch == NULL){
      printf("Incorrect drug sample file. Program aborted!!\n");
      exit(0);
    }
    else{
      idx = 0;
      /* begin splitting sample loop */
      do
      {
        hill[idx] = strtod(pch, NULL);
        printf("hill[%d] = %lf\n", idx, hill[idx]);
        pch = strtok(NULL, ",");
        idx++;
      }
      while( pch != NULL );
      /* end splitting sample loop */
    }

    conc = strtod(argv[3], NULL);
    printf("Concentration: %lf\n", conc);
    if((int)conc == 0){
      //printf("Wrong concentration value!! Applying the normal Patch Clamp!\n");
      is_drug == false;
    }
  }

  
  p_patch = new Ohara_Rudy_2011();
  if( is_drug == true ) p_patch->initConsts(0, conc, hill, false);
  else p_patch->initConsts();
  dt = 0.001;
  tcurr = 0.0;

  // Create CVODE solver
  cvode_main = CVodeCreate( CV_BDF,CV_NEWTON );
  // Give p_cell as CVode User Data
  CVodeSetUserData( cvode_main, p_patch );
  // Create the states vector based on the STATES array
  states_vec = N_VMake_Serial( p_patch->states_size, p_patch->STATES );
  // Initalize CVODE solver
  CVodeInit( cvode_main, rhs_fn, 0.0, states_vec );
  // Set up the future time (must be more than 0.0)
  tnext = dt;
  // Set up the maximum step size (dt)
  CVodeSetMaxStep( cvode_main, dt );
  // Set up the linear solver type
  CVDense( cvode_main, p_patch->states_size );
  // Set up the numerical error tolerances
  CVodeSStolerances( cvode_main, 1.0e-7, 1.0e-7 );

  fp_curr = fopen("curr_ORd_Cab.plt", "w");
  fp_vm = fopen(argv[1], "r");

  fgets( buff, 255, fp_vm );

  fprintf(fp_curr, "%s %s\n", "TIME", "ICab");

  while(fgets( buff, 255, fp_vm ) != NULL){
    pch = strtok(buff, ", ");
    pch = strtok(NULL, ", ");

    // Grab the Vm from input file
    glob_vm = strtod(pch, NULL);

    retval = CVode( cvode_main, tnext, states_vec, &tcurr, CV_NORMAL  );
    if( retval == CV_SUCCESS ){
      tnext += dt;
    }
    else{
      break;
    }

    fprintf(fp_curr, "%lf %lf\n", tcurr, p_patch->ALGEBRAIC[ICab]);

  }

  fclose(fp_curr);
  fclose(fp_vm);


  delete p_patch;
  printf("FINISHED!\n");

  return 0;
}
