#include "AHC_low.h"
//#include "iostream"
//#include <fstream>

//Define output array

//#define DEBUG

#ifdef DEBUG
using std::cout;
using std::endl;
float x1_debug[150][65];
float x2_debug[150][65];
#endif


#define LAST
//using std::cout;
//using std::endl;

AHC::AHC(){
// Partition local variables to improve bandwidth
#pragma HLS ARRAY_PARTITION variable=this->J complete dim=1
#pragma HLS ARRAY_PARTITION variable=this->x complete
#pragma HLS ARRAY_PARTITION variable=this->MVM_out complete
#pragma HLS ARRAY_PARTITION variable=this->e complete

	// Initialize the AHC solver
	this->dt   = 0.01;
	this->r    = 0.98;
	this->beta = 0.1;
	this->coupling_strength = 0.020;
	this->mu = 1.0;
	this->num_time_steps = 200;
	this->target_a_baseline = 0.2;
	this->target_a = target_a_baseline;

	// Initialize the spins, J matrix and MVM output
	// initialize_vectors:
	// for(int i=0;i<N;i++){
	// 	this->e[i] = 1.0;
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
}

// Matrix vector product
void AHC::matmul()
{
	#pragma HLS INLINE off
	// Matrix vector product
	// MVM = (J).dot(np.sign(x))
	for (int row = 0; row < N; row++) {
		#pragma HLS PIPELINE
		this->MVM_out[row] = 0.0;
	}

	// using column method MVM = \sum_j J[:][j] * x[j]
	MVM_outer:
	for(int i=0; i<N; i++){
		// for each element in x
		#pragma HLS PIPELINE
		MVM_inner:
		for(int row = 0; row < N; row++){
			// multiply with each element on i th column of J
			if(this->x[i] == 1){
				this->MVM_out[row] += this->J[row][i];
			}
			else if(this->x[i] == -1){
				this->MVM_out[row] -= this->J[row][i];
			}
			// else{
			// 	this->MVM_out[row] += 0;
			// }
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
			temp = -((this->MVM_out[i]) >> 1);
		}
		else if(this->lastSpins[i]==-1){
			temp = (this->MVM_out[i] >> 1);
		}
		energy += temp;
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
	
	update_spin_and_error:
	for(int i=0; i<N; i++){
		#pragma HLS PIPELINE

		data_type_x xx;
		data_type_e de;
		// dt   = (1 >> 7) + (1 >> 9); // 0.01
		// r    = 1 - (1 >> 6) - (1 >> 8); // 0.98
		// beta = (1 >> 4) + (1 >> 5) + (1 >> 7); // 0.1
		// coupling_strength = (1 >> 6) + (1 >> 8) + (1 >> 11); // 0.020
		// mu = 1.0;
		// num_time_steps = 200;
		// target_a_baseline = (1 >> 3) + (1 >> 4) + (1 >> 6); // 0.2
		// target_a = target_a_baseline;
		
		// Update spin vector
		// this->x[i] += dt * (coupling_strength * this->MVM_out[i])*this->e[i];
		// this->xx[i] = (this->x[i] >> 4);
		// this->x[i] += -dt * this->x[i] * ((data_type_x(0.02)) + mu*this->xx[i]);
		// this->de[i] = dt*(-beta * this->e[i] * (this->xx[i] - target_a));
        // this->e[i] += de[i];

		data_type_x this_x_tmp;
		// Update spin vector
		// this->x[i] += dt * (coupling_strength * this->MVM_out[i])*this->e[i];
		this_x_tmp = (this->MVM_out[i] >> 6) + (this->MVM_out[i] >> 8) + (this->MVM_out[i] >> 11);
		this->x[i] += ((this_x_tmp >> 7) + (this_x_tmp >> 9)) * (this->e[i]);

		// xx[i] = (this->x[i] >> 4);
		xx = (this->x[i] >> 4);

		// this->x[i] += -dt * this->x[i] * ((data_type_x(0.02)) + mu*this->xx[i]);
		data_type_x this_x_tmp_2;
		this_x_tmp_2 = (data_type_x(0.02)) + mu*xx;
		this->x[i] += -((this_x_tmp_2 >> 7) + (this_x_tmp_2 >> 9)) * this->x[i];

		// this->de[i] = dt*(-beta * this->e[i] * (this->xx[i] - target_a));
		data_type_x this_tmp_de;
		// this_tmp_de = -(((xx[i] - target_a) >> 4) + ((xx[i] - target_a) >> 5) + ((xx[i] - target_a) >> 7));
		this_tmp_de = -(((xx - target_a) >> 4) + ((xx - target_a) >> 5) + ((xx - target_a) >> 7));

		de = ((this_tmp_de >> 7) + (this_tmp_de >> 9)) * (this->e[i]);
		this->e[i] += de;
	}
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
	data_type_x x_init[N]
	// ap_fixed<MAX_WIDTH, 2> coupling_strength_new, ap_fixed<MAX_WIDTH,2> new_mu
){
	initialize_vectors:for(int i=0; i<N; i++){
		this->e[i] = 1.0;
		this->x[i] = x_init[i];
		// this->lastSpins[i] = x_init[i];
		// this->de[i] = 0;
		// this->coupling_strength = coupling_strength_new;
		// this->mu = new_mu;
	}
}

void AHC::updateJ(data_type_J J_init[N][N]){
	this->bestEnergy = 10;	// reset bestEnergy
	
	// reset bestSpins
	for (int i=0; i<N; i++){
		this->bestSpins[i] = 0;
	}

	initialize_matrix:
	for(int j=0;j<N;j++){
	#pragma HLS PIPELINE
		for(int i=0; i<N; i++){
            this->J[i][j] = J_init[i][j];
        }
    }
}

void AHC::ahc_solver(
	data_type_x x_init[N]
){
	updateX(x_init);

	setSpins();	// initialize vectors
	matmul();

	TIME_STEP_LOOP:
	for(int time_step=0; time_step < 200; time_step++){
		update();
		setSpins();
		matmul();
		setSpins();
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

// void ahc_top(
// 	data_type_J J_matrix[N][N], 
// 	data_type_x x_init[N], 
// 	spin_sign bestSpinsOut[N]
// ){
// 	// Partition the first dim
// 	#pragma HLS ARRAY_PARTITION variable=J_matrix dim=1 complete
// 	#pragma HLS ARRAY_PARTITION variable=x_init dim=1 complete
// 	#pragma HLS ARRAY_PARTITION variable=bestSpinsOut dim=1 complete

// 	// #pragma HLS INTERFACE bram port=bestSpinsOut
// 	static AHC ahc_instance(J_matrix);
// 	ahc_instance.ahc_solver(x_init, bestSpinsOut);
// }

//----------------------------------------------------------
// Top function
//----------------------------------------------------------

void dut(hls::stream<bit32_t> &strm_in, hls::stream<bit32_t> &strm_out) {
	data_type_J J_in[N][N];
	data_type_x x_in[N];
	spin_sign bestSpinsOut[N];
	
	// Partition the first dim
	#pragma HLS ARRAY_PARTITION variable=J_in dim=1 complete
	#pragma HLS ARRAY_PARTITION variable=x_in dim=1 complete
	#pragma HLS ARRAY_PARTITION variable=bestSpinsOut dim=1 complete

	bit32_t input_l;
	bit32_t output;

	// read J matrix
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {	
			input_l = strm_in.read();
			J_in[i][j] = input_l(15,0);
		}
	}

	// read x_init
	for (int i = 0; i < N; i++) {
		input_l = strm_in.read();
		x_in[i] = input_l(15,0);
	}

	// call ahc
	// ahc_top(J_in, x_in, bestSpinsOut);
	static AHC ahc_instance;
	ahc_instance.updateJ(J_in);
	ahc_instance.ahc_solver(x_in);

	// write out the result
	for (int i = 0; i < N; i++) {
		output(15,0) = bestSpinsOut[i];
		strm_out.write(output);
	}
}