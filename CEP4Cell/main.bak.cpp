#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>

#include "commons/helper.hpp"
#include "cellmodels/enum_courtemanche_ramirez_nattel_1998.hpp"
#include "cellmodels/courtemanche_ramirez_nattel_1998.hpp"

#define T0 RCONST(0.0)
#define T1 RCONST(0.1)
#define DT RCONST(0.1)

static int rhs_fn( realtype t, N_Vector y, N_Vector ydot, void *user_data );
static int AP_sim();

int main( int argc, char **arg )
{
  int err_code;
  param_t *p_param;


  p_param = (param_t*)malloc( sizeof( param_t ) );
  err_code = 0;

  set_default_params(p_param);
  AP_sim();
  return err_code;
}

static int rhs_fn( realtype t, N_Vector y, N_Vector ydot, void *user_data )
{
  CellML *data = (CellML*)user_data;
  data->computeRates( t, 
                      data->CONSTANTS,
                      N_VGetArrayPointer_Serial(ydot),
                      N_VGetArrayPointer_Serial(y),
                      data->ALGEBRAIC ); 
  return 0; 
}

static int AP_sim()
{
  CellML *p_cell;

  int retval;
  double t_out, t;
  void *cvode_mem;
  N_Vector states_vec;
  SUNMatrix states_mat;
  SUNLinearSolver linear_solver;

 // Create CVODE solver
  cvode_mem = CVodeCreate( CV_BDF );
  // Give p_cell as CVode User Data
  p_cell = new courtemanche_ramirez_nattel_1998();
  p_cell->initConsts();
  CVodeSetUserData( cvode_mem, p_cell );
  // Create the states vector based on the STATES array
  states_vec = N_VMake_Serial( p_cell->states_size, p_cell->STATES );
  // Initalize CVODE solver
  CVodeInit( cvode_mem, rhs_fn, T0, states_vec );
  // Set up the maximum step size (dt)
  CVodeSetMaxStep( cvode_mem, 0.1 );
  // Set up the linear solver
  states_mat = SUNDenseMatrix( p_cell->states_size, p_cell->states_size );
  linear_solver = SUNLinSol_Dense( states_vec, states_mat );
  CVodeSetLinearSolver( cvode_mem, linear_solver, states_mat );
  // Set up the tolerances
  CVodeSStolerances( cvode_mem, 1.0e-7, 1.0e-7 );

  // Set up the starting time (must be more than 0.0)
  t_out = T1;
  int iout = 0;
  FILE *fptr = fopen( "test.plt", "w" );
  while(1){
    retval = CVode( cvode_mem, t_out, states_vec, &t, CV_NORMAL  );
    fprintf( fptr, "%lf %lf %lf\n", 
             t, 
             p_cell->STATES[V], 
             p_cell->ALGEBRAIC[i_Kr] );
    if( retval == CV_SUCCESS ){
      iout++;
      t_out += 0.1;
    }
    if (iout >= 10000) break;
  }
  fclose( fptr );

  // Memory Cleanup
  N_VDestroy(states_vec);
  SUNLinSolFree(linear_solver);
  SUNMatDestroy(states_mat);
  CVodeFree(&cvode_mem);
  delete p_cell;

  return 0;
}
