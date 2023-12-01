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

AHC::AHC(data_type_x x_init[N], data_type_J J_init[N][N]){
// Partition local variables to improve bandwidth
// #pragma HLS ARRAY_PARTITION variable=J complete dim=0
// #pragma HLS ARRAY_PARTITION variable=x complete
#pragma HLS ARRAY_PARTITION variable=this->MVM_out complete
#pragma HLS ARRAY_PARTITION variable=this->e complete

	// Initialize the AHC solver
	dt   = 0.01;
	r    = 0.98;
	beta = 0.1;
	coupling_strength = 0.020;
	mu = 1.0;
	num_time_steps = 200;
	target_a_baseline = 0.2;
	target_a = target_a_baseline;

	// Initialize the spins, J matrix and MVM output
	initialize_vectors:
	for(int i=0;i<N;i++){
		e[i] = 1.0;
		this->x[i] = x_init[i];
		MVM_out[i] = 0.0;
		this->lastSpins[i] = x_init[i];
	}

    initialize_matrix:
	for(int j=0;j<N;j++){
	#pragma HLS PIPELINE
		for(int i=0;i<N;i++){
            this->J[i][j] = J_init[i][j];
        }
    }
}

// returns the square of the input val
data_type_x square_single(data_type_x val) {
	data_type_x tmp_xx;
	// use LUT to perform mul, instead of using DSP
	// #pragma HLS BIND_OP variable=tmp_xx op=mul impl=fabric
	tmp_xx = val * val;
	return(tmp_xx);
}

// setSpins
// This function update the lastSpins based on the current spin values
// This is used to determine the sign of the spins in the MVM
void AHC::setSpins(){
	#pragma HLS INLINE
	setSpins_loop:
	for(int i =0; i < N; i++){
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

void matmul_inner(data_type_J J[N][N], data_type_x x[N], data_type_x MVM_out[N]) {
	// #pragma HLS INLINE
	// Matrix vector product
	// MVM = (J).dot(np.sign(x))
	data_type_J J_tmp[N][N];
	#pragma HLS ARRAY_PARTITION variable=J_tmp complete dim=2
	#pragma HLS ARRAY_PARTITION variable=J complete dim=2
	#pragma HLS ARRAY_PARTITION variable=x complete
	for(int i = 0; i < N; ++i) {
		#pragma HLS pipeline
		for (int j = 0; j < N; ++j) {
			J_tmp[i][j] = J[i][j];
		}
	}

	for (int i = 0; i < N; ++i) {
		#pragma HLS PIPELINE
		MVM_out[i] = 0.0;
	}
	MVM_outer:
	for(int i = 0; i < N; i++){
		#pragma HLS PIPELINE
		data_type_x tmp = 0.0;
		MVM_inner:
		for(int j = 0; j < N; j++){
			if(x[j] == 1){
				tmp += J_tmp[i][j];
			}
			else if(x[j] == -1){
				tmp -= (J_tmp[i][j]);
				//this->lastSpins[i] = -1;
			}
			else{
				tmp += 0;
				// this->lastSpins[i] = 0;
			}
		}
		MVM_out[i] += tmp;
	}
}

// Matrix vector product
void AHC::matmul()
{
	// #pragma HLS INLINE
	// // Matrix vector product
	// // MVM = (J).dot(np.sign(x))
	// for (int i = 0; i < N; ++i) {
	// 	#pragma HLS PIPELINE
	// 	this->MVM_out[i] = 0.0;
	// }
	// MVM_outer:
	// for(int i = 0; i < N; i++){
	// 	#pragma HLS PIPELINE
	// 	data_type_x tmp = 0.0;
	// 	MVM_inner:
	// 	for(int j = 0; j < N; j++){
	// 		if(this->x[j] == 1){
	// 			tmp += this->J[i][j];
	// 		}
	// 		else if(this->x[j] == -1){
	// 			tmp -= (this->J[i][j]);
	// 			//this->lastSpins[i] = -1;
	// 		}
	// 		else{
	// 			tmp += 0;
	// 			// this->lastSpins[i] = 0;
	// 		}
	// 	}
	// 	this->MVM_out[i] += tmp;
	// }
	matmul_inner(this->J, this->x, this->MVM_out);
}



// Calculates the Ising energy
data_type_e AHC::IsingEnergy(){
	data_type_e energy = 0.0;
	IsingEnergy_loop: 
	for(int i = 0; i < N; i++){
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
	return energy;
}

// Update the spins and error vectors
void AHC::update(){
	#pragma HLS INLINE
	// #pragma HLS LATENCY min=4 max=10

	update_spin_and_error:
	for(int i=0;i<N;i++){
		#pragma HLS PIPELINE
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

		this->xx[i] = (this->x[i] >> 4);

		// this->x[i] += -dt * this->x[i] * ((data_type_x(0.02)) + mu*this->xx[i]);
		data_type_x this_x_tmp_2;
		this_x_tmp_2 = (data_type_x(0.02)) + mu*this->xx[i];
		this->x[i] += -((this_x_tmp_2 >> 7) + (this_x_tmp_2 >> 9)) * this->x[i];

		// this->de[i] = dt*(-beta * this->e[i] * (this->xx[i] - target_a));
		data_type_x this_tmp_de;
		this_tmp_de = -(((this->xx[i] - target_a) >> 4) + ((this->xx[i] - target_a) >> 5) + ((this->xx[i] - target_a) >> 7));
		this->de[i] = ((this_tmp_de >> 7) + (this_tmp_de >> 9)) * (this->e[i]);
		this->e[i] += de[i];
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

void AHC::updateSpins(data_type_x x_init[N], ap_fixed<MAX_WIDTH, 2> coupling_strength_new, ap_fixed<MAX_WIDTH,2> new_mu){
	initialize_vectors:for(int i=0;i<N;i++){
		this->e[i] = 1.0;
		this->x[i] = x_init[i];
		this->MVM_out[i] = 0.0;
		this->lastSpins[i] = x_init[i];
		this->de[i] = 0;
		this->coupling_strength = coupling_strength_new;
		this->mu = new_mu;
	}
}

void AHC::ahc_solver(){
	#pragma HLS ARRAY_PARTITION variable=xx dim=0 complete
	#pragma HLS ARRAY_PARTITION variable=dx dim=0 complete
	#pragma HLS ARRAY_PARTITION variable=dx2 dim=0 complete
	#pragma HLS ARRAY_PARTITION variable=de dim=0 complete

	data_type_e energy = 0.0;

	setSpins();	// initialize vectors
	matmul();

	TIME_STEP_LOOP:
	for(int time_step=0; time_step < 200; time_step++){
		update();
		setSpins();
		matmul();
		setSpins();
		energy = IsingEnergy();
		if(energy < this->bestEnergy){
			bestEnergy = energy;
			for(int k = 0; k < N; k++){
				bestSpins[k] = this->lastSpins[k];
			}
		}
	}
}

void ahc_top(
	data_type_J J_matrix[N][N], 
	data_type_x x_init[N], 
	spin_sign bestSpinsOut[N]
){
	// Partition the first dim
	// #pragma HLS ARRAY_PARTITION variable=J_matrix dim=1 complete
	// #pragma HLS ARRAY_PARTITION variable=x_init dim=1 complete

	#pragma HLS INTERFACE bram port=bestSpinsOut
	static AHC ahc_instance(x_init, J_matrix);
	ahc_instance.ahc_solver();

    for (int i = 0; i < N; i++) {
        bestSpinsOut[i] = ahc_instance.bestSpins[i];
    }
}

// add a new dut function
