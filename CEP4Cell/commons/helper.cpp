#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "helper.hpp"

void set_default_values( param_t *p_param )
{
  p_param->celltype = 0;
  p_param->simulation_mode = 0;
  p_param->t_write_vtk = 0.0;
  p_param->dt = 0.1;
  p_param->dt_write = 1.0;
  p_param->stim_amp = 1.0;
  p_param->stim_dur = 0.5;
  p_param->num_pace1 = 2;
  p_param->bcl_decrement = 50.0;
  p_param->bcl_end = 2000000.0;
  p_param->bcl_init = 1000.0;
  p_param->mutation_scaling = 0.0;
  p_param->v_apd90 = -40.0;
  p_param->is_print_vm = 1;
  p_param->is_s1s2 = 0;
  p_param->is_using_output = 0;
  sprintf( p_param->variant, "%s", "ori" );
  sprintf( p_param->drug_name, "%s", "bepridil" );
  p_param->cmax = 33;
  sprintf( p_param->dosages, "%s", "1.0" );
}

double get_cmax( char *drug_name )
{
  double cmax;
  if(strcasecmp(drug_name, "quinidine") == 0){
    cmax = 3237.;
  }
  else if(strcasecmp(drug_name, "bepridil") == 0){
    cmax = 33.;
  }
  else if(strcasecmp(drug_name, "dofetilide") == 0){
    cmax = 2.;
  }
  else if(strcasecmp(drug_name, "sotalol") == 0){
    cmax = 14690.;
  }
  else if(strcasecmp(drug_name, "chlorpromazine") == 0){
    cmax = 38.;
  }
  else if(strcasecmp(drug_name, "cisapride") == 0){
    cmax = 2.6;
  }
  else if(strcasecmp(drug_name, "terfenadine") == 0){
    cmax = 4.;
  }
  else if(strcasecmp(drug_name, "ondansetron") == 0){
    cmax = 139.;
  }
  else if(strcasecmp(drug_name, "diltiazem") == 0){
    cmax = 122.;
  }
  else if(strcasecmp(drug_name, "mexiletine") == 0){
    cmax = 4129.;
  }
  else if(strcasecmp(drug_name, "ranolazine") == 0){
    cmax = 1948.2;
  }
  else if(strcasecmp(drug_name, "verapamil") == 0){
    cmax = 81.;
  }
  else if(strcasecmp(drug_name, "amiodarone") == 0){
    cmax = 0.7;
  }
  else if(strcasecmp(drug_name, "kitmoxi") == 0){
    cmax = 5979.;
  }
  else if(strcasecmp(drug_name, "kitranolazine") == 0){
    cmax = 1948.2;
  }
  else if(strcasecmp(drug_name, "kitverapamil") == 0){
    cmax = 81.;
  }
  else if(strcasecmp(drug_name, "control") == 0){
    cmax = 1.;
  }
  else{
    fprintf(stdout, "Drug name %s is not listed! Cmax set to 0!!\n", drug_name);
    cmax = 0.;
  }
  return cmax;
}


void assign_params_single(int argc, char **argv, param_t *p_param)
{

  for (int k = 2; k <= argc; k += 2) {
    if (strcmp(argv[k - 1], "-ionmodel") == 0) {
      strcpy(p_param->ionmodel, argv[k]);
    }
    else if (strcmp(argv[k - 1], "-twrite") == 0) {
      p_param->t_write_vtk = strtod(argv[k], NULL);
    }
    else if (strcmp(argv[k - 1], "-dt") == 0) {
      p_param->dt = strtod(argv[k], NULL);
    }
    else if (strcmp(argv[k - 1], "-stim_dur") == 0) {
      p_param->stim_dur = strtod(argv[k], NULL);
    }
    else if (strcmp(argv[k - 1], "-stim_amp") == 0) {
      p_param->stim_amp = strtod(argv[k], NULL);
    }
    else if (strcmp(argv[k - 1], "-num_pace1") == 0) {
      p_param->num_pace1 = strtol(argv[k], NULL, 10);
    }
    else if (strcmp(argv[k - 1], "-bcl_init") == 0) {
      p_param->bcl_init = strtod(argv[k], NULL);
    }
    else if (strcmp(argv[k - 1], "-bcl_end") == 0) {
      p_param->bcl_end = strtod(argv[k], NULL);
    }
    else if (strcmp(argv[k - 1], "-bcl_decrement") == 0) {
      p_param->bcl_decrement = strtod(argv[k], NULL);
    }
    else if (strcmp(argv[k - 1], "-v_apd90") == 0) {
      p_param->v_apd90 = strtod(argv[k], NULL);
    }
    else if (strcmp(argv[k - 1], "-mutation_scaling") == 0) {
      p_param->mutation_scaling = strtod(argv[k], NULL);
    }
    else if (strcmp(argv[k - 1], "-variant") == 0) {
      strcpy(p_param->variant, argv[k]);
    }
    else if (strcmp(argv[k - 1], "-simulation_mode") == 0) {
      p_param->simulation_mode = strtol( argv[k], NULL, 10 );
    }
    else if (strcmp(argv[k - 1], "-scaled_gate") == 0) {
      p_param->simulation_mode = strtol( argv[k], NULL, 10 );
    }
   

   else {
      fprintf( stderr, "%s: %s!!!\n %s.\n", "Unknown argument", argv[k - 1],
               "Check manual for details" );
    }
  }

}

void edison_assign_params_single(int argc, char *args[], param_t *p_param)
{
  char buffer[100];
  char key[100];
  char value[100];
  char file_name[255];
  FILE *fp_inputdeck;

  for (int i = 1; i < argc; i += 2) {
    if (!strcmp(args[i], "-input_deck"))
      strcpy(file_name, args[i + 1]);
    else if (!strcmp(args[i], "-hill_file"))
      strcpy(p_param->hill_file, args[i + 1]);
    else if (!strcmp(args[i], "-herg_file"))
      strcpy(p_param->herg_file, args[i + 1]);
  }

  if(strlen(p_param->hill_file) != 0){
    printf( "Hill File: %s\n", p_param->hill_file );
  }
  if(strlen(p_param->herg_file) != 0){
    printf( "hERG File: %s\n", p_param->herg_file );
  }
 
  fp_inputdeck = fopen( file_name, "r");
  if(fp_inputdeck == NULL){
    fprintf(stderr, "Cannot open file %s!!!\n", file_name);
    exit(0);
  }

  // read input_deck line by line
  // and store each line to the buffer
  while ( fgets( buffer, 100, fp_inputdeck ) != NULL ) {
    // parse the buffer and store it to key and value
    sscanf( buffer, "%s %*s %s", key, value );
    if (strcasecmp(key, "Celltype") == 0) {
      p_param->celltype = strtod( value, NULL );
      printf( "%s -- %lf\n", "Celltype", p_param->celltype);
    } 
    else if (strcasecmp(key, "Variant") == 0) {
      strcpy( p_param->variant, value );
      printf( "%s -- %s\n", "Variant", p_param->variant);
    }
    else if (strcasecmp(key, "Simulation_Mode") == 0) {
      p_param->simulation_mode = strtol( value, NULL, 10 );
      printf( "%s -- %d\n", "Simulation_Mode", p_param->simulation_mode);
    } 
    else if (strcasecmp(key, "Write_Time") == 0) {
      p_param->t_write_vtk = strtod( value, NULL );
      printf( "%s -- %lf\n", "Write_Time", p_param->t_write_vtk);
    }
    else if (strcasecmp(key, "Number_Of_Pacing") == 0) {
      p_param->num_pace1 = strtol( value, NULL, 10 );
      printf( "%s -- %d\n", "Number_of_Pacing", p_param->num_pace1);
    }
    else if (strcasecmp(key, "Basic_Cycle_Length") == 0) {
      p_param->bcl_init = strtod( value, NULL );
      printf( "%s -- %lf\n", "Basic_Cycle_Length", p_param->bcl_init);
    }
    else if (strcasecmp(key, "Time_Step") == 0) {
      p_param->dt = strtod( value, NULL );
      printf( "%s -- %lf\n", "Time_Step", p_param->dt);
    }
    else if (strcasecmp(key, "Write_Step") == 0) {
      p_param->dt_write = strtod( value, NULL );
      printf( "%s -- %lf\n", "Writing_Step", p_param->dt_write);
    }
    else if (strcasecmp(key, "Stim_Amplitude") == 0) {
      p_param->stim_amp = strtod( value, NULL );
      printf( "%s -- %lf\n", "Stim_Amplitude", p_param->stim_amp);
    }
    else if (strcasecmp(key, "Stim_Duration") == 0) {
      p_param->stim_dur = strtod( value, NULL );
      printf( "%s -- %lf\n", "Stim_Duration", p_param->stim_dur);
    }
    else if (strcasecmp(key, "BCL_Step") == 0) {
      p_param->bcl_decrement = strtod( value, NULL );
      printf( "%s -- %lf\n", "BCL_Step", p_param->bcl_decrement);
    }
    else if (strcasecmp(key, "BCL_End") == 0) {
      p_param->bcl_end = strtod( value, NULL );
      printf( "%s -- %lf\n", "BCL_End", p_param->bcl_end);
    }
    else if (strcasecmp(key, "APD_90_Percent") == 0) {
      p_param->v_apd90 = strtod( value, NULL );
      printf( "%s -- %lf\n", "APD_90_Percent", p_param->v_apd90);
    }
    else if (strcasecmp(key, "Scaled_Gate") == 0) {
      p_param->scaled_gate = strtol( value, NULL, 10 );
      printf( "%s -- %d\n", "Scaled_Gate", p_param->scaled_gate);
    }
    else if (strcasecmp(key, "Drug_Name") == 0) {
      strcpy( p_param->drug_name, value );
      printf( "%s -- %s\n", "Drug_Name", p_param->drug_name);
    }
    else if (strcasecmp(key, "Cmax") == 0) {
      p_param->cmax = strtod(value, NULL);
      printf( "%s -- %lf\n", "Cmax", p_param->cmax);
    }
    else if (strcasecmp(key, "Dosages") == 0) {
      strcpy( p_param->dosages, "0," );
      strcat( p_param->dosages, value );
      printf( "%s -- %s\n", "Dosages", p_param->dosages);
    }
    else if (strcasecmp(key, "GKs_Mult") == 0) {
      p_param->gks_mult =  strtod(value, NULL) ;
      printf( "%s -- %lf\n", "GKs_Mult", p_param->gks_mult);
    }
    else if (strcasecmp(key, "GCaL_Mult") == 0) {
      p_param->gcal_mult =  strtod(value, NULL) ;
      printf( "%s -- %lf\n", "GCaL_Mult", p_param->gcal_mult);
    }
    else if (strcasecmp(key, "GNa_Mult") == 0) {
      p_param->gna_mult =  strtod(value, NULL) ;
      printf( "%s -- %lf\n", "GNa_Mult", p_param->gna_mult);
    }
    else if (strcasecmp(key, "Is_Print_Vm") == 0) {
      p_param->is_print_vm = strtol( value, NULL, 10 );
      printf( "%s -- %d\n", "Is_Print_Vmcheck", p_param->is_print_vm );
    }
    else if (strcasecmp(key, "Is_Using_Output") == 0) {
      p_param->is_using_output = strtol( value, NULL, 10 );
      printf( "%s -- %d\n", "Is_Using_Output", p_param->is_using_output );
    }


  }
  fclose( fp_inputdeck );
}

int save_last_state(CellML *p_cell, double t_laststate)
{
  FILE *fp_laststate;
  char buff[20];
  int i;

  fp_laststate = fopen("output.dat", "w");
  fprintf(fp_laststate, "%lf\n", t_laststate);
  for(i = 0; i < p_cell->states_size; i++)
  {
    fprintf(fp_laststate, "%lf\n", p_cell->STATES[i]);
  }
  fclose(fp_laststate);
  //printf("Write initial state from the last %lf ms simulation...\n", t_laststate);

  return 0;
}

int load_last_state(CellML *p_cell)
{
  FILE *fp_laststate;
  double t_laststate;
  char buff[20];
  int i;

  fp_laststate = fopen("output.dat", "r");
  if(fp_laststate != NULL){
    i = 0;
    fgets(buff, sizeof(buff), fp_laststate);
    t_laststate = strtod(buff, NULL);
    printf("Get initial state from the last %lf ms simulation...\n", t_laststate);
    while(fgets(buff, sizeof(buff), fp_laststate) != NULL)
    {
      p_cell->STATES[i++] = strtod(buff, NULL);
    }
    fclose(fp_laststate);
  }
  else{
    printf("Cannot find output.dat in this path\n");
  }

  return 0;
}




CellML* init_cell()
{
  CellML *p_cell;
  char model_ID[255];

  p_cell = NULL;
#if defined(CRN1998)
  p_cell = new courtemanche_ramirez_nattel_1998();
  strcpy( model_ID, "crn1998" );
#elif defined(HHUXLEY1952)
  p_cell = new hodgkin_huxley_squid_axon_model_1952();
  strcpy( model_ID, "HHuxley1952" );
#elif defined(OHARA_RUDY2011)
  p_cell = new Ohara_Rudy_2011();
  strcpy( model_ID, "ORudy2011" );
#elif defined(ORUDY_CIPA2017)
  p_cell = new ohara_rudy_cipa_v1_2017();
  strcpy( model_ID, "ORudy_CiPA2017" );
#elif defined(TN2004ENDO)
  p_cell = new tentusscher_noble_noble_panfilov_2004_a();
  strcpy( model_ID, "tn2004endo" );
#elif defined(TN2004EPI)
  p_cell = new tentusscher_noble_noble_panfilov_2004_b();
  strcpy( model_ID, "tn2004epi" );
#elif defined(TN2004M)
  p_cell = new tentusscher_noble_noble_panfilov_2004_c();
  strcpy( model_ID, "tn2004m" );
#elif defined(TN2006ENDO)
  p_cell = new ten_tusscher_model_2006_IK1Ko_endo_units();
  strcpy( model_ID, "tn2006endo" );
#elif defined(TN2006EPI)
  p_cell = new ten_tusscher_model_2006_IK1Ko_epi_units();
  strcpy( model_ID, "tn2006epi" );
#elif defined(TN2006M)
  p_cell = new ten_tusscher_model_2006_IK1Ko_M_units();
  strcpy( model_ID, "tn2006M" );
#endif


  if (p_cell != NULL) {
    printf( "Using cell model %s\n", model_ID );
  }
  else {
    fprintf( stderr, "Please make sure you use correct model ID.\n" );
    fprintf( stderr, "%s\n",
             "If you have difficulties, contact Marcell or Prof. Lim." );
    exit(EXIT_FAILURE);
  }

  return p_cell;
}


CellML *init_mutation(const char variant_ID[], double scale, CellML *p_cell)
{
  if( p_cell == NULL ){
    fprintf( stderr, "Make sure the p_cell object is already initialized.\n" );
    exit(EXIT_FAILURE);
  }

  if (strcasecmp(variant_ID, "L532P") == 0) {
#ifdef CRN1998
      p_cell->mutation = new hERG("L532P");
#endif
  }
  else if (strcasecmp(variant_ID, "N588K") == 0) {
#ifdef CRN1998
      p_cell->mutation = new hERG("N588K");
#endif
  }
  else if (strcasecmp(variant_ID, "S140G") == 0) {
#ifdef CRN1998
      p_cell->mutation = new S140G(scale);
#endif
  }
  else if (strcasecmp(variant_ID, "Fibrosis") == 0) {
#ifdef CRN1998
      p_cell->mutation = new Fibrosis();
#endif
  }
  else if (strcasecmp(variant_ID, "A1656D") == 0) {
#if defined TN2004ENDO || defined TN2004M || defined TN2004EPI
      p_cell->mutation = new A1656D("A1656D");
#endif
  }
  else if (strcasecmp(variant_ID, "A1656D_WT") == 0) {
#if defined TN2004ENDO || defined TN2004M || defined TN2004EPI
      p_cell->mutation = new A1656D("WT");
#endif
  }
  else if (strcasecmp(variant_ID, "A1656D_MEX") == 0) {
#if defined TN2004ENDO || defined TN2004M || defined TN2004EPI
      p_cell->mutation = new A1656D("MEX");
#endif
  }
  else if (strcasecmp(variant_ID, "A1656D_FLE") == 0) {
#if defined TN2004ENDO || defined TN2004M || defined TN2004EPI
      p_cell->mutation = new A1656D("FLE");
#endif
  }
  else if (strcasecmp(variant_ID, "A1656D_RAN") == 0) {
#if defined TN2004ENDO || defined TN2004M || defined TN2004EPI
      p_cell->mutation = new A1656D("RAN");
#endif
  }
  else if (strcasecmp(variant_ID, "G229D_WT") == 0) {
#if defined OHARA_RUDY2011 || defined ORUDY_CIPA2017
      p_cell->mutation = new G229D("WT");
#endif
  }
  else if (strcasecmp(variant_ID, "G229D") == 0) {
#if defined OHARA_RUDY2011 || defined ORUDY_CIPA2017
      p_cell->mutation = new G229D("G229D");
#endif
  }



  if (p_cell->mutation != NULL) {
    printf( "Apply mutation %s\n", variant_ID );
  }
  else{
    fprintf( stderr, "Please make sure the mutation ID %s is registered.\n", variant_ID );
    exit(EXIT_FAILURE);
  }

  return p_cell;
}
