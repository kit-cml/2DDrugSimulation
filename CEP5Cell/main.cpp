#include "modules/ep.hpp"
#include "modules/drugs_hill.hpp"
#include "modules/population_ep.hpp"

#include <cstdio>
#include <cstdlib>

int main( int argc, char **argv )
{
  int err_code;
  param_t *p_param;

  p_param = (param_t*)malloc( sizeof( param_t ) );
  err_code = 0;

  set_default_values(p_param);

  edison_assign_params_single(argc, argv, p_param);

  /* check the modules folder to look the implementation */
  if( p_param->simulation_mode == 0 ){
    err_code = EP( argc, argv, p_param );
  }
#if defined ORUDY2011_STATIC
  else if( p_param->simulation_mode == 4 ){
    err_code = drugs_hill( argc, argv, p_param, false);
  }
  else if( p_param->simulation_mode == 8 ){
    err_code = drugs_hill( argc, argv, p_param, true);
  }
#endif
#if defined ORUDY2017_DYNAMIC
  else if( p_param->simulation_mode == 5 ){
    err_code = drugs_herg( argc, argv, p_param);
  }
#endif
  else{
    printf("Selection unknown!!\n");
  }

  return err_code;
}
