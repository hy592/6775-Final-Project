#include "AHC_low.h"
//#include "iostream"
//#include <fstream>

//#define DEBUG
#ifdef DEBUG
using std::cout;
using std::endl;
float x1_debug[150][65];
float x2_debug[150][65];
#endif

AHC::AHC(){
#pragma HLS INLINE

// Partition local variables to improve bandwidth
#pragma HLS ARRAY_PARTITION variable=J complete dim=1
#pragma HLS ARRAY_PARTITION variable=x complete
#pragma HLS ARRAY_PARTITION variable=MVM_out complete
#pragma HLS ARRAY_PARTITION variable=error_var complete

	// Initialize the AHC solver
	this->dt   = 0.01;	// time_step
	this->r    = 0.98;
	this->beta = 1;
	this->gamma = 1/(256*0.01);
	this->coupling_str = 100;
	this->mu = 1;
	this->num_time_steps = 256;
	this->target_a_baseline = 0;
	this->target_a = target_a_baseline;

	// Initialize the spins, J matrix and MVM output
	// initialize_vectors:
	// for(int i=0;i<N;i++){
	// 	this->error_var[i] = 1.0;
	// 	// this->x[i] = x_init[i];
	// }

    // initialize_matrix:
	// for(int j=0;j<N;j++){
	// #pragma HLS PIPELINE
	// 	for(int i=0;i<N;i++){
    //         this->J[i][j] = J_init[i][j];
    //     }
    // }
}

// setSpins
// This function update the lastSpins based on the current spin values
// This is used to determine the sign of the spins in the MVM
void AHC::setSpins(){
	#pragma HLS INLINE off
	setSpins_loop:
	for(int i=0; i<N; i++){
		#pragma HLS PIPELINE
		if(this->x[i] > 0){
			this->lastSpins[i] = 1;
		}
		else if(this->x[i] < 0){
			this->lastSpins[i] = -1;
		}
		else{
			this->lastSpins[i] = 0;
		}
	}
	//cout << endl;
}

// Matrix vector product
void AHC::Mat_Vec_Mal()
{
	#pragma HLS INLINE off
	// Matrix vector product
	// MVM = (J).dot(np.sign(x))
	MVM_init:
	for (int row = 0; row < N; row++) {
		#pragma HLS PIPELINE
		this->MVM_out[row] = 0.0;
	}

	// using column method MVM = \sum_j J[:][j] * x[j]
	MVM_outer:
		// for each element in x
	for(int j = 0; j < N; j++){
		#pragma HLS PIPELINE II=1
		MVM_inner:
		for(int i=0; i<N; i++){
			#pragma HLS UNROLL
			// multiply with each element on i th column of J
			if(this->lastSpins[j] == 1){
				this->MVM_out[i] += this->J[i][j];
			}
			else if(this->lastSpins[j] == -1){
				this->MVM_out[i] += -(this->J[i][j]);
			}
			else {
				// this->MVM_out[i] += 0;
				continue;
			}
		}
	}
}

// Calculates the Ising energy
void AHC::IsingEnergy(){
	#pragma HLS INLINE off
	data_type_e energy = 0.0;
	IsingEnergy_loop: 
	for(int i=0; i<N; i++){
		#pragma HLS PIPELINE
		data_type_e temp;
		if(this->lastSpins[i]==1){
			energy += -((this->MVM_out[i]))/2;
		}
		else if(this->lastSpins[i]==-1){
			energy += (this->MVM_out[i])/2;
		}
	}

	if(energy < this->bestEnergy){
		this->bestEnergy = energy;
		for(int k=0; k<N; k++){
			this->bestSpins[k] = this->lastSpins[k];
		}
	}
}

// Update the spins and error vectors
void AHC::update(){
	#pragma HLS INLINE off
	// #pragma HLS LATENCY min=4 max=10
	// data_type_x xx[N];
	// data_type_e de[N];

	// #pragma HLS ARRAY_PARTITION variable=xx dim=0 complete
	// #pragma HLS ARRAY_PARTITION variable=de dim=0 complete
	
	//int spinFlips = 0;
	//cout << "x update IS ";

	update_spin_and_error:
	for(int i=0; i<N; i++){
		#pragma HLS PIPELINE
		// Evolve Spin Variables
        // xx = x**2
		this->xx[i] = this->x[i]*this->x[i];
        // x += (x*((r-1) - mu*xx)) * time_step
		this->x[i] += (this->x[i]*((this->r-1) - this->mu*this->xx[i])) * this->dt;
        // x += coupling_str*(MVM*error_var) * time_step
		this->x[i] += this->coupling_str * (this->MVM_out[i]*this->error_var[i]) * this->dt;

		// Evolve Error Variables
		// error_var += -(beta*((xx) - target_a)*error_var) * time_step
		this->error_var[i] += -(this->beta*((this->xx[i]) - this->target_a)*this->error_var[i]) * this->dt;
	}

	// if modulate_parameter=="target_a":
	// 	#Modulate Target Amplitude
	// 	# delta_a = coupling_str*np.mean((J.dot(sig))*sig*etc_flag)
	// 	delta_a = gamma*time_step
	// 	target_a  = target_a_baseline + delta_a
	this->target_a = this->target_a_baseline + this->gamma * this->dt;
	// elif modulate_parameter=="coupling_str":
	// 	coupling_str += time_step * gamma
	// this->coupling_str += this->dt * this->gamma;
}

void AHC::reset(){
	#pragma HLS INLINE off
	// Reset MVM
	reset_MVM:
	for(int i = 0; i < N; i++){
		#pragma HLS PIPELINE
		this->MVM_out[i] = 0.0;
	}
}

void AHC::updateX(
	data_type_x x_init[N],
	data_type_e error_var_init[N]
	// ap_fixed<MAX_WIDTH, 2> coupling_strength_new, ap_fixed<MAX_WIDTH,2> new_mu
){
	initialize_vectors:for(int i=0; i<N; i++){
		// x = np.sqrt(0.001) * np.random.normal(0,1,(N, num_anneals))
    	// error_var = np.abs(np.sqrt(0.001)* np.random.normal(0,1,(N, num_anneals)))
		// if (x_init[N-i-1]>=0) {
		// 	this->error_var[i] = x_init[N-i-1];
		// }
		// else{
		// 	this->error_var[i] = -x_init[N-i-1];
		// }
		this->error_var[i] = error_var_init[i];
		// this->error_var[i] = 0.015;
		this->x[i] = x_init[i];
		// this->lastSpins[i] = x_init[i];
		// this->de[i] = 0;
		// this->coupling_str = coupling_strength_new;
		// this->mu = new_mu;
	}
	this->target_a = target_a_baseline;
}

void AHC::updateJ(data_type_J J_init[N][N]){
	this->bestEnergy = 10;	// reset bestEnergy
	
	// reset bestSpins
	initialize_matrix_X:
	for (int i=0; i<N; i++){
		this->bestSpins[i] = 0;
	}

	initialize_matrix_J:
	for(int i=0; i<N; i++){
		for(int j=0;j<N;j++){
            this->J[i][j] = J_init[i][j];
        }
    }
}

void AHC::ahc_solver(
	data_type_x x_init[N],
	data_type_e error_var_init[N]
){
	updateX(x_init, error_var_init);

	setSpins();	// initialize vectors
	// Mat_Vec_Mal();

	Mat_Vec_Mal();
	TIME_STEP_LOOP:
	for(int time_step=0; time_step < num_time_steps; time_step++){
		update();
		setSpins();
		Mat_Vec_Mal();
		IsingEnergy();
	}
}

data_type_e AHC::bestEnergySpins(
	spin_sign bestSpins[N]
){
	for(int i=0; i<N; i++){
		bestSpins[i] = this->bestSpins[i];
	}

	return this->bestEnergy;
}
//----------------------------------------------------------
// Top function
//----------------------------------------------------------

// void dut(hls::stream<bit32_t> &strm_in, hls::stream<bit32_t> &strm_out) {
// 	#pragma HLS INTERFACE axis port=strm_in
// 	#pragma HLS INTERFACE axis port=strm_out

// 	data_type_J J_in[N][N];
// 	data_type_x x_in[N];
// 	data_type_x bestSpins[N];

// 	bit32_t input_1;
// 	bit32_t spin_output;

// 	// Read J matrix from the input stream
// 	for (int i = 0; i < N; i++) {
// 		for (int j = 0; j < N; j++) {
// 			input_1 = strm_in.read();
// 			J_in[i][j] = input_1(MAX_WIDTH-1,0);
// 		}
// 	}

// 	// Read x array from the input stream
// 	for (int i = 0; i < N; i++) {
// 		input_1 = strm_in.read();
// 		x_in[i] = input_1(MAX_WIDTH-1,0);
// 	}

// 	// Run the AHC solver
// 	static AHC ahc_instance;
// 	ahc_instance.updateJ(J_in);  // Assuming this function is adapted to accept J[N][N]

// 	ahc_instance.ahc_solver(x_in);

// 	// Get the best spins and write to the output stream
// 	//ahc_instance.bestEnergySpins(bestSpins);  // Assuming this function is adapted to write to bestSpins[N]
// 	for (int i = 0; i < N; i++) {
// 		spin_output = bestSpins[i];
// 		strm_out.write(spin_output);
// 	}
// }

#if USE_6775_FPGA
void dut(hls::stream<bit32_t> &strm_in, hls::stream<bit32_t> &strm_out) {
	data_type_J J_in[N][N];
	data_type_x x_in[N];
	spin_sign bestSpinsOut[N];
	data_type_e bestEnergy;

	bit32_t input_l;
	bit32_t output_energy;
	bit32_t output_spin;

	static AHC ahc_instance;
	static data_type_e error_var_init[N] = {
		// load the error_var_init array with random values
		0.05242187853736398,0.008994496034482512,0.05480038244046159,0.014535002581538762,0.030323180132189084,0.02071990675188693,0.01213663603773126,0.00960427400309387,0.017547715916464597,0.024983051674625344,0.03859390568213095,0.034408063058388484,0.0636682426201058,0.0026453633012165404,0.004943585304561,0.007960975245106408,0.005659684029419596,0.02215511338451815,0.001726779419402328,0.025851092401670953,0.026165299804272325,0.02988735373792587,0.006131210329798413,0.009222319281435049,0.015215984098860743,0.03222454678140047,0.04533535815500804,0.02579067020782812,0.002089247779940635,0.01987630401193418,0.020922910140258277,0.07817702324671731,0.020167845558607785,0.04794447537786903,0.02262055299674045,0.04876876714536372,0.06158700291924554,0.02247247608951715,0.029853129896046538,0.031169784715797493,0.023949098756832352,0.015571072434411177,0.037654580840756374,0.008904025901565473,0.00916569097574489,0.043537247439120214,0.02526344484950409,0.03788821512993865,0.011778126741224991,0.05141053224292928,0.001778414510275391,0.0011715480839723674,0.012272890148973526,0.010841586182997884,0.030748680876433288,0.0504170902169587,0.0066111367253571885,0.04569360064802776,0.014188751421206104,0.001515377842436009,0.012149846868276077,0.0027996573416023776,0.053971628955546094,0.05138659062773086,0.015673517622038453
	};

	// read J matrix
	for (int i = 0; i < N; i++) {
		#pragma HLS pipeline off
		for (int j = 0; j < N; j++) {	
			#pragma HLS pipeline
			input_l = strm_in.read();
			data_type_J J_receive;
			J_receive(MAX_WIDTH-1,0) = input_l(MAX_WIDTH-1,0);
			J_in[i][j] = J_receive;
			// cout << "J " << J_receive << endl;
		}
	}
	ahc_instance.updateJ(J_in);

	// run num_anneals sets of X
	for (int x_iter=0; x_iter<num_anneals; x_iter++){
		#pragma HLS pipeline off
		// read x_init
		for (int i = 0; i < N; i++) {
			#pragma HLS pipeline
			input_l = strm_in.read();
			data_type_x X_receive;
			X_receive(MAX_WIDTH-1,0) = input_l(MAX_WIDTH-1,0);
			x_in[i] = X_receive;
			// cout << "x " << X_receive << endl;
		}
		ahc_instance.ahc_solver(x_in, error_var_init);
	}


	// return the best energy
	bestEnergy = ahc_instance.bestEnergySpins(bestSpinsOut);
	output_energy(MAX_WIDTH-1,0) = bestEnergy(MAX_WIDTH-1,0);
	strm_out.write(output_energy);

	// write out the result
	for (int i = 0; i < 8; i++) {
		// pack 8 spins into one transmit
		output_spin(1,0)   = reinterpret_cast<bit2_t&>(bestSpinsOut[8*i  ]);
		output_spin(3,2)   = reinterpret_cast<bit2_t&>(bestSpinsOut[8*i+1]);
		output_spin(5,4)   = reinterpret_cast<bit2_t&>(bestSpinsOut[8*i+2]);
		output_spin(7,6)   = reinterpret_cast<bit2_t&>(bestSpinsOut[8*i+3]);
		output_spin(9,8)   = reinterpret_cast<bit2_t&>(bestSpinsOut[8*i+4]);
		output_spin(11,10) = reinterpret_cast<bit2_t&>(bestSpinsOut[8*i+5]);
		output_spin(13,12) = reinterpret_cast<bit2_t&>(bestSpinsOut[8*i+6]);
		output_spin(15,14) = reinterpret_cast<bit2_t&>(bestSpinsOut[8*i+7]);
		strm_out.write(output_spin);
	}
	output_spin(1,0)   = reinterpret_cast<bit2_t&>(bestSpinsOut[64]);
	strm_out.write(output_spin);
}
#else
void dut(
	data_type_J J_matrix[N][N], 
	data_type_x x_init[N], 
	spin_sign bestSpinsOut[N], 
	data_type_e bestEnergyOut
){
	#pragma HLS INTERFACE s_axilite port=return // export the control interface on the top function
	#pragma HLS INTERFACE m_axi port=J_matrix offset=slave bundle=gmem
	#pragma HLS INTERFACE m_axi port=x_init offset=slave bundle=gmem
	#pragma HLS INTERFACE m_axi port=bestSpinsOut offset=slave bundle=gmem
	#pragma HLS INTERFACE m_axi port=bestEnergyOut offset=slave bundle=gmem

	// Run the AHC solver
	static AHC ahc_instance;
	static data_type_e error_var_init[N] = {
		// load the error_var_init array with random values
		0.05242187853736398,0.008994496034482512,0.05480038244046159,0.014535002581538762,0.030323180132189084,0.02071990675188693,0.01213663603773126,0.00960427400309387,0.017547715916464597,0.024983051674625344,0.03859390568213095,0.034408063058388484,0.0636682426201058,0.0026453633012165404,0.004943585304561,0.007960975245106408,0.005659684029419596,0.02215511338451815,0.001726779419402328,0.025851092401670953,0.026165299804272325,0.02988735373792587,0.006131210329798413,0.009222319281435049,0.015215984098860743,0.03222454678140047,0.04533535815500804,0.02579067020782812,0.002089247779940635,0.01987630401193418,0.020922910140258277,0.07817702324671731,0.020167845558607785,0.04794447537786903,0.02262055299674045,0.04876876714536372,0.06158700291924554,0.02247247608951715,0.029853129896046538,0.031169784715797493,0.023949098756832352,0.015571072434411177,0.037654580840756374,0.008904025901565473,0.00916569097574489,0.043537247439120214,0.02526344484950409,0.03788821512993865,0.011778126741224991,0.05141053224292928,0.001778414510275391,0.0011715480839723674,0.012272890148973526,0.010841586182997884,0.030748680876433288,0.0504170902169587,0.0066111367253571885,0.04569360064802776,0.014188751421206104,0.001515377842436009,0.012149846868276077,0.0027996573416023776,0.053971628955546094,0.05138659062773086,0.015673517622038453
	};

	ahc_instance.updateJ(J_matrix);  // Assuming this function is adapted to accept J[N][N]
	ahc_instance.ahc_solver(x_init, error_var_init);

	// Get the best spins and write to the output stream
	bestEnergyOut = ahc_instance.bestEnergySpins(bestSpinsOut);  // Assuming this function is adapted to write to bestSpins[N]
}
#endif