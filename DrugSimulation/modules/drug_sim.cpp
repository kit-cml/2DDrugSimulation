#include "drug_sim.hpp"

#ifdef TOMEK_2019
#include "../cellmodels/Tomek_model.hpp"
#elif ORUDY_2017
#include "../cellmodels/ohara_rudy_cipa_v1_2017.hpp"
#else
#include "../cellmodels/Ohara_Rudy_2011.hpp"
#endif

#include "glob_funct.hpp"

#include <cmath>
#include <cstdlib>

bool do_drug_sim(const double conc, row_data ic50, row_data2 herg_pars,
const param_t* p_param, const unsigned short sample_id, const unsigned short group_id, qinward_t *p_qin)
{
	bool is_ead = false;
	unsigned short idx;
	
	// to normalize the small values
	// so it can compressed the file size.
	// Information of the scaling should be
	// described in the header result.
	static const int CALCIUM_SCALING = 1000000;
	static const int CURRENT_SCALING = 1000;

	// cell object pointer
	Cellmodel* p_cell;

	// qnet_ap/inet_ap values
	double inet_ap, qnet_ap, inet4_ap, qnet4_ap, inet_cl, qnet_cl, inet4_cl, qnet4_cl;
	double inal_auc_ap, ical_auc_ap,inal_auc_cl, ical_auc_cl, qinward_cl;

	// variables for I/O
	char buffer[255];
	FILE* fp_output;
	FILE* fp_qni;
	FILE* fp_vmdebug;
	
	FILE* fp_time_series;
	FILE* fp_time_series_all;
	
	// simulation parameters
#ifdef DEBUG_MODE
	bool is_print_graph = true;
	bool is_dutta = false;
	const char *drug_name = "bepridil";
	const double bcl = 2000.;
	const double inet_vm_threshold = -88.0;
	const unsigned short pace_max = 1000;
	const unsigned short celltype = 0.;
	const unsigned short last_drug_check_pace = 250;
	const unsigned int print_freq = (1./dt) * dtw;
	unsigned short pace_count = 0;
	unsigned short pace_steepest = 0;
#else
	bool is_print_graph = p_param->is_print_graph;
	bool is_dutta = p_param->is_dutta;
	const char *drug_name = p_param->drug_name;
	const double bcl = p_param->bcl;
	const double inet_vm_threshold = p_param->inet_vm_threshold;
	const unsigned short pace_max = p_param->pace_max;
	const unsigned short celltype = p_param->celltype;
	const unsigned short last_drug_check_pace = p_param->last_drug_check_pace;
	unsigned short pace_count = 0;
	unsigned short pace_steepest = 0;
#endif

	// variables to store features
	// temp_result is the result of features in 1 pace,
	// will be interchanged during the simulation.
	// cipa_result is the final result of the simulation.
	cipa_t cipa_result, temp_result;
	// eligible AP shape means the Vm_peak > 0.
	bool is_eligible_AP;
	// Vm value at 30% repol, 50% repol, and 90% repol, respectively.
	double vm_repol30, vm_repol50, vm_repol90;

	// // these values are from the supplementary materials of ORd2011
	double dt = 0.005;
	double max_time_step = 1.0;
	double time_point = 25.0;
	double tcurr = 0.0;
	double tmax = pace_max*bcl;
	double t_peak_capture = 0.0;
	double dt_set;

#ifdef TOMEK_2019
	printf("Using Tomek cell model\n");
	p_cell = new Tomek_model();
	p_cell->initConsts((double)celltype, conc, ic50.data);
#elif ORUDY_2017
	printf("Using ORd2017-dyn cell model\n");
	p_cell = new ohara_rudy_cipa_v1_2017();
	//p_cell->initConsts();
	p_cell->initConsts((double)celltype, conc, ic50.data, herg_pars.data);
	dt = 0.01; // maximum dt for RK4 method.
#else
	printf("Using O'Hara Rudy cell model\n");
	p_cell = new Ohara_Rudy_2011();
	p_cell->initConsts((double)celltype, conc, ic50.data, is_dutta);
#endif
	p_cell->CONSTANTS[BCL] = bcl;

	FILE *fp_states;
	if( p_param->is_using_output > 0 ){
		fp_states = fopen("output_orudy.dat", "r");
		if( fp_states != NULL ){
			mpi_printf(0, "Using initial condition from steady state!\n");
			int idx = 0;
			while(fgets( buffer, sizeof(buffer), fp_states) != NULL){
				p_cell->STATES[idx++] = strtod(buffer, NULL);
			}
		}
		else{
			mpi_printf(0, "No initial file found! Skipped using initial file!\n");
		}
	}
	
	snprintf(buffer, sizeof(buffer), "result/%.4lf/%s_%.4lf_vmdebug_smp%hu.plt", 
			  conc, drug_name, conc, sample_id );
	fp_vmdebug = fopen( buffer, "w" );
	snprintf(buffer, sizeof(buffer), "result/%.4lf/%s_%.4lf_qni_proc%hu.plt", 
			  conc, drug_name, conc, mympi::rank );
	fp_qni = fopen( buffer, "a" );
	
	snprintf(buffer, sizeof(buffer), "result/%.4lf/%s_%.4lf_time_series_result_smp%hu.plt", 
			  conc, drug_name, conc, sample_id );
	fp_time_series = fopen( buffer, "w" );
/*
	snprintf(buffer, sizeof(buffer), "result/%.4lf/%s_%.4lf_time_series_result_all_smp%hu.plt", 
			  conc, drug_name, conc, sample_id );
	fp_time_series_all = fopen( buffer, "w" );	
*/
	if(group_id == 0){
		fprintf( fp_qni, "%s,%s,%s,%s,%s,%s,%s\n","Sample_ID", "Qnet_AP", "Qnet4_AP", "Qinward_CL","Qnet_CL", "Qnet4_CL", "Qinward_CL");
	}

	fprintf(fp_vmdebug, "%s,%s,%s,%s,%s,%s,%s\n", "Pace","T_Peak", "Vmpeak","Vm_repol30","Vm_repol50","Vm_repol90","dVmdt_repol");
	fprintf(fp_time_series,"%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n",
			"Time","Vm","dVm/dt","Cai(x1.000.000)(milliM->picoM)",
			"INa(x1.000)(microA->picoA)","INaL(x1.000)(microA->picoA)","ICaL(x1.000)(microA->picoA)",
			"IKs(x1.000)(microA->picoA)","IKr(x1.000)(microA->picoA)","IK1(x1.000)(microA->picoA)",
			"Ito(x1.000)(microA->picoA)");
/*
	fprintf(fp_time_series_all,"%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n",
			"Time","Vm","dVm/dt","Cai(x1.000.000)(milliM->picoM)",
			"INa(x1.000)(microA->picoA)","INaL(x1.000)(microA->picoA)","ICaL(x1.000)(microA->picoA)",
			"IKs(x1.000)(microA->picoA)","IKr(x1.000)(microA->picoA)","IK1(x1.000)(microA->picoA)",
			"Ito(x1.000)(microA->picoA)");
*/	
	inet_ap = 0.;
	qnet_ap = 0.;
	inet4_ap = 0.;
	qnet4_ap = 0.;
	inal_auc_ap = 0.;
	ical_auc_ap = 0.;
	qinward_cl = 0.;
	inet_cl = 0.;
	qnet_cl = 0.;
	inet4_cl = 0.;
	qnet4_cl = 0.;
	inal_auc_cl = 0.;
	ical_auc_cl = 0.;
	qinward_cl = 0.;
	pace_count = 0;
	
	unsigned int icount = 0;
	unsigned int imax = (unsigned int)((pace_max * bcl)/dt);
	double dtw = 2.;
	const unsigned int print_freq = (1./dt) * dtw;
	
	cipa_result.init( p_cell->STATES[V]);
	temp_result.init( p_cell->STATES[V]);	

#ifdef ORUDY_2017
	while(icount < imax){
		// compute ODE at tcurr
		p_cell->computeRates(tcurr,
					p_cell->CONSTANTS,
					p_cell->RATES,
					p_cell->STATES,
					p_cell->ALGEBRAIC);
		
		// new pace code
		if(icount % (int)(bcl/dt) == 0 && icount > 0){
			if( is_eligible_AP && pace_count >= pace_max-last_drug_check_pace) {
				temp_result.qnet_ap = qnet_ap;
				temp_result.qnet4_ap = qnet4_ap;
				temp_result.inal_auc_ap = inal_auc_ap;
				temp_result.ical_auc_ap = ical_auc_ap;
				temp_result.qnet_cl = qnet_cl;
				temp_result.qnet4_cl = qnet4_cl;
				temp_result.inal_auc_cl = inal_auc_cl;
				temp_result.ical_auc_cl = ical_auc_cl;
				fprintf(fp_vmdebug, "%hu,%.2lf,%.2lf,%.2lf,%.2lf,%.2lf,%.2lf\n", pace_count,t_peak_capture,temp_result.vm_peak,vm_repol30,vm_repol50,vm_repol90,temp_result.dvmdt_repol);
				// replace result with steeper repolarization AP or first pace from the last 250 paces
				if( temp_result.dvmdt_repol > cipa_result.dvmdt_repol ) {
					pace_steepest = pace_count;
					cipa_result = temp_result;
				}
			};
			inet_ap = 0.;
			qnet_ap = 0.;
			inet4_ap = 0.;
			qnet4_ap = 0.;
			inal_auc_ap = 0.;
			ical_auc_ap = 0.;
			inet_cl = 0.;
			qnet_cl = 0.;
			inet4_cl = 0.;
			qnet4_cl = 0.;
			inal_auc_cl = 0.;
			ical_auc_cl = 0.;
			t_peak_capture = 0.;
			temp_result.init( p_cell->STATES[V]);	
			pace_count++;
			is_eligible_AP = false;
		}
		
		// Compute the solution.
		// Analytical method not available in ORd2017-dyn,
		// hence in the meantime, using RK4 method (10x faster than Euler).
		//p_cell->solveAnalytical(dt);
		p_cell->solveRK4(tcurr,dt);
		
		// begin the last 250 pace operations
		if(pace_count >= pace_max-last_drug_check_pace){
			// Find peak vm around 2 msecs and  40 msecs after stimulation
			// and when the sodium current reach 0
			if( tcurr > ((p_cell->CONSTANTS[BCL]*pace_count)+(p_cell->CONSTANTS[stim_start]+2)) && 
				tcurr < ((p_cell->CONSTANTS[BCL]*pace_count)+(p_cell->CONSTANTS[stim_start]+10)) && 
				abs(p_cell->ALGEBRAIC[INa]) < 1){
				if( p_cell->STATES[V] > temp_result.vm_peak ){
					temp_result.vm_peak = p_cell->STATES[V];
					if(temp_result.vm_peak > 0){
						vm_repol30 = temp_result.vm_peak - (0.3 * (temp_result.vm_peak - temp_result.vm_valley));
						vm_repol50 = temp_result.vm_peak - (0.5 * (temp_result.vm_peak - temp_result.vm_valley));
						vm_repol90 = temp_result.vm_peak - (0.9 * (temp_result.vm_peak - temp_result.vm_valley));
						is_eligible_AP = true;
						t_peak_capture = tcurr;
					}
					else is_eligible_AP = false;
				}
			}
			else if( tcurr > ((p_cell->CONSTANTS[BCL]*pace_count)+(p_cell->CONSTANTS[stim_start]+10)) && is_eligible_AP ){
				if( p_cell->RATES[V] > temp_result.dvmdt_repol &&
					p_cell->STATES[V] <= vm_repol30 &&
					p_cell->STATES[V] >= vm_repol90 ){
					temp_result.dvmdt_repol = p_cell->RATES[V];
				}
				
			}
			
			if(is_eligible_AP && p_cell->STATES[V] > vm_repol90){
				// inet_ap/qnet_ap under APD.
				inet_ap = (p_cell->ALGEBRAIC[INaL]+p_cell->ALGEBRAIC[ICaL]+p_cell->ALGEBRAIC[Ito]+p_cell->ALGEBRAIC[IKr]+p_cell->ALGEBRAIC[IKs]+p_cell->ALGEBRAIC[IK1]);
				inet4_ap = (p_cell->ALGEBRAIC[INaL]+p_cell->ALGEBRAIC[ICaL]+p_cell->ALGEBRAIC[IKr]+p_cell->ALGEBRAIC[INa]);
				qnet_ap += (inet_ap * dt)/1000.;
				qnet4_ap += (inet4_ap * dt)/1000.;
				inal_auc_ap += (p_cell->ALGEBRAIC[INaL]*dt);
				ical_auc_ap += (p_cell->ALGEBRAIC[ICaL]*dt);
			}
			// inet_ap/qnet_ap under Cycle Length
			inet_cl = (p_cell->ALGEBRAIC[INaL]+p_cell->ALGEBRAIC[ICaL]+p_cell->ALGEBRAIC[Ito]+p_cell->ALGEBRAIC[IKr]+p_cell->ALGEBRAIC[IKs]+p_cell->ALGEBRAIC[IK1]);
			inet4_cl = (p_cell->ALGEBRAIC[INaL]+p_cell->ALGEBRAIC[ICaL]+p_cell->ALGEBRAIC[IKr]+p_cell->ALGEBRAIC[INa]);
			qnet_cl += (inet_cl * dt)/1000.;
			qnet4_cl += (inet4_cl * dt)/1000.;
			inal_auc_cl += (p_cell->ALGEBRAIC[INaL]*dt);
			ical_auc_cl += (p_cell->ALGEBRAIC[ICaL]*dt);
					
			if(icount % print_freq == 0){
				snprintf( buffer, sizeof(buffer), "%.2lf,%.2lf,%.0lf,%.0lf,%.0lf,%.0lf,%0.lf,%.0lf,%.0lf,%.0lf",
						p_cell->STATES[V], p_cell->RATES[V], p_cell->STATES[cai]*CALCIUM_SCALING,
						p_cell->ALGEBRAIC[INa]*CURRENT_SCALING, p_cell->ALGEBRAIC[INaL]*CURRENT_SCALING, 
						p_cell->ALGEBRAIC[ICaL]*CURRENT_SCALING, p_cell->ALGEBRAIC[Ito]*CURRENT_SCALING,
						p_cell->ALGEBRAIC[IKr]*CURRENT_SCALING, p_cell->ALGEBRAIC[IKs]*CURRENT_SCALING, 
						p_cell->ALGEBRAIC[IK1]*CURRENT_SCALING);
				temp_result.time_series_data.insert( std::pair<double, string> (tcurr, string(buffer)) );
			}

		}// end the last 250 pace operations

/*
		if( pace_count >= 0 && pace_count <= 10 && icount % print_freq == 0){
		fprintf( fp_time_series_all, "%.4lf,%.2lf,%.2lf,%.0lf,%.0lf,%.0lf,%.0lf,%0.lf,%.0lf,%.0lf,%.0lf\n",
						tcurr,p_cell->STATES[V], p_cell->RATES[V], p_cell->STATES[cai]*CALCIUM_SCALING,
						p_cell->ALGEBRAIC[INa]*CURRENT_SCALING, p_cell->ALGEBRAIC[INaL]*CURRENT_SCALING, 
						p_cell->ALGEBRAIC[ICaL]*CURRENT_SCALING, p_cell->ALGEBRAIC[Ito]*CURRENT_SCALING,
						p_cell->ALGEBRAIC[IKr]*CURRENT_SCALING, p_cell->ALGEBRAIC[IKs]*CURRENT_SCALING, 
						p_cell->ALGEBRAIC[IK1]*CURRENT_SCALING);
				temp_result.time_series_data.insert( std::pair<double, string> (tcurr, string(buffer)) );
		}
*/
		tcurr += dt;
		icount++;
	}
#else
	while(tcurr < tmax)
	{
		// compute ODE at tcurr
		p_cell->computeRates(tcurr,
					p_cell->CONSTANTS,
					p_cell->RATES,
					p_cell->STATES,
					p_cell->ALGEBRAIC);
		dt_set = Ohara_Rudy_2011::set_time_step(tcurr,
				time_point,
				max_time_step,
				p_cell->CONSTANTS,
				p_cell->RATES,
				p_cell->STATES,
				p_cell->ALGEBRAIC);
		// compute accepted timestep
		if (floor((tcurr + dt_set) / bcl) == floor(tcurr / bcl)) {
		  dt = dt_set;
		}
		// new cycle length code.
		// this is the place for comparing current and previous paces.
		// if the AP shape is eligible and dvmdt_repol is bigger,
		// that AP shape become the resultant pace.
		// also for re-initializing stuffs.
		else {
			dt = (floor(tcurr / bcl) + 1) * bcl - tcurr;
			if( is_eligible_AP && pace_count >= pace_max-last_drug_check_pace) {
				temp_result.qnet_ap = qnet_ap;
				temp_result.qnet4_ap = qnet4_ap;
				temp_result.inal_auc_ap = inal_auc_ap;
				temp_result.ical_auc_ap = ical_auc_ap;
				temp_result.qnet_cl = qnet_cl;
				temp_result.qnet4_cl = qnet4_cl;
				temp_result.inal_auc_cl = inal_auc_cl;
				temp_result.ical_auc_cl = ical_auc_cl;
				fprintf(fp_vmdebug, "%hu,%.2lf,%.2lf,%.2lf,%.2lf,%.2lf,%.2lf\n", pace_count,t_peak_capture,temp_result.vm_peak,vm_repol30,vm_repol50,vm_repol90,temp_result.dvmdt_repol);
				// replace result with steeper repolarization AP or first pace from the last 250 paces
				if( temp_result.dvmdt_repol > cipa_result.dvmdt_repol ) {
					pace_steepest = pace_count;
					cipa_result = temp_result;
				}
			};
			inet_ap = 0.;
			qnet_ap = 0.;
			inet4_ap = 0.;
			qnet4_ap = 0.;
			inal_auc_ap = 0.;
			ical_auc_ap = 0.;
			inet_cl = 0.;
			qnet_cl = 0.;
			inet4_cl = 0.;
			qnet4_cl = 0.;
			inal_auc_cl = 0.;
			ical_auc_cl = 0.;
			t_peak_capture = 0.;
			temp_result.init( p_cell->STATES[V]);	
			pace_count++;
			is_eligible_AP = false;
		}
		//Compute the analytical solution
		p_cell->solveAnalytical(dt);
		//printf("%lf\t%lf\n",tcurr,p_cell->STATES[V]);
		
		// begin the last 250 pace operations
		if(pace_count >= pace_max-last_drug_check_pace){
			// Find peak vm around 2 msecs and  40 msecs after stimulation
			// and when the sodium current reach 0
			if( tcurr > ((p_cell->CONSTANTS[BCL]*pace_count)+(p_cell->CONSTANTS[stim_start]+2)) && 
				tcurr < ((p_cell->CONSTANTS[BCL]*pace_count)+(p_cell->CONSTANTS[stim_start]+10)) && 
				abs(p_cell->ALGEBRAIC[INa]) < 1){
				if( p_cell->STATES[V] > temp_result.vm_peak ){
					temp_result.vm_peak = p_cell->STATES[V];
					if(temp_result.vm_peak > 0){
						vm_repol30 = temp_result.vm_peak - (0.3 * (temp_result.vm_peak - temp_result.vm_valley));
						vm_repol50 = temp_result.vm_peak - (0.5 * (temp_result.vm_peak - temp_result.vm_valley));
						vm_repol90 = temp_result.vm_peak - (0.9 * (temp_result.vm_peak - temp_result.vm_valley));
						is_eligible_AP = true;
						t_peak_capture = tcurr;
					}
					else is_eligible_AP = false;
				}
			}
			else if( tcurr > ((p_cell->CONSTANTS[BCL]*pace_count)+(p_cell->CONSTANTS[stim_start]+10)) && is_eligible_AP ){
				if( p_cell->RATES[V] > temp_result.dvmdt_repol &&
					p_cell->STATES[V] <= vm_repol30 &&
					p_cell->STATES[V] >= vm_repol90 ){
					temp_result.dvmdt_repol = p_cell->RATES[V];
				}
				
			}
			
			if(is_eligible_AP && p_cell->STATES[V] > vm_repol90){
				// inet_ap/qnet_ap under APD.
				inet_ap = (p_cell->ALGEBRAIC[INaL]+p_cell->ALGEBRAIC[ICaL]+p_cell->ALGEBRAIC[Ito]+p_cell->ALGEBRAIC[IKr]+p_cell->ALGEBRAIC[IKs]+p_cell->ALGEBRAIC[IK1]);
				inet4_ap = (p_cell->ALGEBRAIC[INaL]+p_cell->ALGEBRAIC[ICaL]+p_cell->ALGEBRAIC[IKr]+p_cell->ALGEBRAIC[INa]);
				qnet_ap += (inet_ap * dt)/1000.;
				qnet4_ap += (inet4_ap * dt)/1000.;
				inal_auc_ap += (p_cell->ALGEBRAIC[INaL]*dt);
				ical_auc_ap += (p_cell->ALGEBRAIC[ICaL]*dt);
			}
			// inet_ap/qnet_ap under Cycle Length
			inet_cl = (p_cell->ALGEBRAIC[INaL]+p_cell->ALGEBRAIC[ICaL]+p_cell->ALGEBRAIC[Ito]+p_cell->ALGEBRAIC[IKr]+p_cell->ALGEBRAIC[IKs]+p_cell->ALGEBRAIC[IK1]);
			inet4_cl = (p_cell->ALGEBRAIC[INaL]+p_cell->ALGEBRAIC[ICaL]+p_cell->ALGEBRAIC[IKr]+p_cell->ALGEBRAIC[INa]);
			qnet_cl += (inet_cl * dt)/1000.;
			qnet4_cl += (inet4_cl * dt)/1000.;
			inal_auc_cl += (p_cell->ALGEBRAIC[INaL]*dt);
			ical_auc_cl += (p_cell->ALGEBRAIC[ICaL]*dt);
			
			// save temporary result		
			snprintf( buffer, sizeof(buffer), "%.2lf,%.2lf,%.0lf,%.0lf,%.0lf,%.0lf,%0.lf,%.0lf,%.0lf,%.0lf",
					p_cell->STATES[V], p_cell->RATES[V], p_cell->STATES[cai]*CALCIUM_SCALING,
					p_cell->ALGEBRAIC[INa]*CURRENT_SCALING, p_cell->ALGEBRAIC[INaL]*CURRENT_SCALING, 
					p_cell->ALGEBRAIC[ICaL]*CURRENT_SCALING, p_cell->ALGEBRAIC[Ito]*CURRENT_SCALING,
					p_cell->ALGEBRAIC[IKr]*CURRENT_SCALING, p_cell->ALGEBRAIC[IKs]*CURRENT_SCALING, 
					p_cell->ALGEBRAIC[IK1]*CURRENT_SCALING);
			temp_result.time_series_data.insert( std::pair<double, string> (tcurr, string(buffer)) );
		} // end the last 250 pace operations
		
		tcurr += dt;
	}
#endif
	printf("time series size outer loop: %d\n", cipa_result.time_series_data.size());

	if( cipa_result.dvmdt_repol > 0 ) is_ead = true;
	
	// print result from the map to the file.
	// first_itrmap is used to normalize the time unit.
	std::multimap<double, string>::iterator first_itrmap_str = cipa_result.time_series_data.begin();
	for(std::multimap<double, string>::iterator itrmap = cipa_result.time_series_data.begin(); itrmap != cipa_result.time_series_data.end() ; itrmap++ ){
		fprintf(fp_time_series, "%.4lf,%s\n", itrmap->first-first_itrmap_str->first, (itrmap->second).c_str());
	}

	fprintf(fp_vmdebug,"Steepest pace: %hu\n", pace_steepest);
	
	// for qnet_ap and qinward_cl output
    if( (int)ceil(conc) == 0 ) {
      p_qin->inal_auc_control = cipa_result.inal_auc_cl;
      p_qin->ical_auc_control = cipa_result.ical_auc_cl;
      fprintf( fp_qni, "%hu,%lf,%lf,%lf,%lf,%lf,%lf\n", sample_id, cipa_result.qnet_ap, cipa_result.qnet4_ap, 0.0, cipa_result.qnet_cl, cipa_result.qnet4_cl, 0.0);
    }
    else{
      p_qin->inal_auc_drug = cipa_result.inal_auc_cl;
      p_qin->ical_auc_drug = cipa_result.ical_auc_cl;
      qinward_cl =  ( (p_qin->inal_auc_drug/p_qin->inal_auc_control) + (p_qin->ical_auc_drug/p_qin->ical_auc_control) ) * 0.5;
      fprintf( fp_qni, "%hu,%lf,%lf,%lf,%lf,%lf,%lf\n", sample_id, cipa_result.qnet_ap, cipa_result.qnet4_ap, qinward_cl, cipa_result.qnet_cl, cipa_result.qnet4_cl, qinward_cl );
    }


	// clean the memories
	//fclose(fp_time_series_all);
	fclose(fp_time_series);
	fclose(fp_vmdebug);
	fclose(fp_qni);

	delete p_cell;


	return is_ead;
}