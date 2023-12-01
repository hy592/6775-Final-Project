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

AHC::AHC(data_type_J J_init[N][N]){
// Partition local variables to improve bandwidth
#pragma HLS ARRAY_PARTITION variable=this->J complete
#pragma HLS ARRAY_PARTITION variable=this->x complete
#pragma HLS ARRAY_PARTITION variable=this->MVM complete
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
	initialize_MVM_vectors:
	for(int i=0;i<N;i++){
		e[i] = 1.0;
		// this->x[i] = x_init[i];
		MVM[i] = 0.0;
		// this->lastSpins[i] = x_init[i];
	}

    initialize_J_matrix:
	for(int j=0;j<N;j++){
	#pragma HLS PIPELINE
		J_init:for(int i=0;i<N;i++){
            this->J[i][j] = J_init[i][j];
        }
    }
}

// returns the square of the input val
// data_type_x square_single(data_type_x val) {
// 	data_type_x tmp_xx;
// 	// use LUT to perform mul, instead of using DSP
// 	#pragma HLS BIND_OP variable=tmp_xx op=mul impl=fabric
// 	tmp_xx = val * val;
// 	return(tmp_xx);
// }

// setSpins
// This function update the lastSpins based on the current spin values
// This is used to determine the sign of the spins in the MVM
void AHC::setSpins(
	data_type_x x_in[N],
	spin_sign lastSpins_out[N]
){
	// input x
	// output lastSpins
	#pragma HLS inline off
	#pragma HLS pipeline
	setSpins_loop:
	for(int i =0; i < N; i++){
		if(this->x[i] > 0){
			lastSpins_out[i] = 1;
		}
		else if(this->x[i] < 0){
			lastSpins_out[i] = -1;
		}
		else{
			lastSpins_out[i] = 0;
		}
	}
}

// Matrix vector product
void AHC::matmul(
	data_type_x x_in[N],
	data_type_x MVM_new[N]
){
	#pragma HLS inline off
	// Matrix vector product
	// MVM = (J).dot(np.sign(x))
	data_type_x MVM_temp[N];
	MVM_init:for(int i=0; i<N; i++) {
		#pragma HLS unroll
		MVM_temp[i]=0;
	}

	MVM_inner:for(int j=0; j<N; j++){
		#pragma HLS pipeline rewind
		#pragma HLS unroll factor=8
		MVM_outer:for(int i=0; i<N; i++){
			if(x_in[j]==1){
				MVM_temp[i] += this->J[i][j];
			}
			else if(x_in[j]==-1){
				MVM_temp[i] -= (this->J[i][j]);
				//this->lastSpins[i] = -1;
			}
			else{
				MVM_temp[i] += 0;
				// this->lastSpins[i] = 0;
			}
		}
	}
	update_MVM:for(int i=0; i<N; i++){
		MVM_new[i] = MVM_temp[i];
	}
}

// Calculates the Ising energy
void AHC::IsingEnergy(
	data_type_x MVM_in[N],
	spin_sign lastSpins_in[N]
){
	#pragma HLS inline off
	data_type_e energy = 0.0;
	// calculate the new energy
	IsingEnergy_loop: for(int i = 0; i < N; i++){
		#pragma HLS pipeline rewind
		#pragma HLS unroll factor=8
		data_type_e temp;
		if(lastSpins_in[i]==1){
			temp = -((MVM_in[i]) >> 1);
		}
		else if(lastSpins_in[i]==-1){
			temp = (MVM_in[i] >> 1);
		}
		energy += temp;
	}

	// update bestEnergy and bestSpins
	if(energy < this->bestEnergy){
		this->bestEnergy = energy;
		bestSpins:for(int k = 0; k < N; k++){
			this->bestSpins[k] = lastSpins_in[k];
		}
	}
}

// Update the spins and error vectors
void AHC::update(
	data_type_x x_old[N],
	data_type_x MVM_in[N],
	data_type_x x_new[N]
){
	#pragma HLS inline off
	// #pragma HLS LATENCY min=4 max=10

	update_spin_and_error:
	for(int i=0;i<N;i++){
		#pragma HLS pipeline rewind
		#pragma HLS unroll factor=8
		data_type_x this_x_tmp;
		data_type_x xx, de;

		// Update spin vector
		// this->x[i] += dt * (coupling_strength * this->MVM_out[i])*this->e[i];
		this_x_tmp = (MVM_in[i] >> 6) + (MVM_in[i] >> 8) + (MVM_in[i] >> 11);
		x_old[i] += ((this_x_tmp >> 7) + (this_x_tmp >> 9)) * (this->e[i]);

		xx = (x_old[i] >> 4);

		// this->x[i] += -dt * this->x[i] * ((data_type_x(0.02)) + mu*this->xx[i]);
		data_type_x this_x_tmp_2;
		this_x_tmp_2 = (data_type_x(0.02)) + mu*xx;
		x_old[i] += -((this_x_tmp_2 >> 7) + (this_x_tmp_2 >> 9)) * x_old[i];

		// this->de[i] = dt*(-beta * this->e[i] * (this->xx[i] - target_a));
		data_type_x this_tmp_de;
		this_tmp_de = -(((xx - target_a) >> 4) + ((xx - target_a) >> 5) + ((xx - target_a) >> 7));
		de = ((this_tmp_de >> 7) + (this_tmp_de >> 9)) * (this->e[i]);
		this->e[i] += de;

		x_new[i] = x_old[i];
	}
}

// void AHC::reset(){
// 	#pragma HLS INLINE off
// 	#pragma HLS PIPELINE
// 	// Reset MVM
// 	reset_MVM:for(int i=0;i<N;i++){
// 		this->MVM_out[i] = 0.0;
// 	}
// }

// void AHC::updateSpins(data_type_x x_init[N], ap_fixed<MAX_WIDTH, 2> coupling_strength_new, ap_fixed<MAX_WIDTH,2> new_mu){
// 	initialize_vectors:for(int i=0;i<N;i++){
// 		this->e[i] = 1.0;
// 		this->x[i] = x_init[i];
// 		this->MVM_out[i] = 0.0;
// 		this->lastSpins[i] = x_init[i];
// 		this->de[i] = 0;
// 		this->coupling_strength = coupling_strength_new;
// 		this->mu = new_mu;
// 	}
// }

void Split (
	data_type_x in[N], 
	data_type_x out1[N], 
	data_type_x out2[N]
) {
	// Duplicated data
	L1:for(int i=1;i<N;i++) {
		out1[i] = in[i]; 
		out2[i] = in[i]; 
	}
}

void AHC::ahc_solver(
	data_type_x x_init[N]
){
	data_type_e energy = 0.0;
	data_type_x x_old[N];
	data_type_x x_new[N];
	data_type_x MVM_new[N];
	
	// setSpins();	// initialize vectors
	solver_init:for(int i=0; i<N; i++) {
		this->e[i] = 1.0;
		x_old[i] = x_init[i];
	}

	matmul(x_old, MVM_new);

	TIME_STEP_LOOP:
	for(int time_step=0; time_step < 200; time_step++){
		data_type_x x_new_1[N], x_new_2[N];

		update(x_old, MVM_new, x_new);
		Split(x_new, x_new_1, x_new_2);

		#pragma HLS dataflow
		// in parallel
		spin_sign lastSpins[N];
		setSpins(x_new_1, lastSpins);
		matmul(x_new_2, MVM_new);

		// setSpins();
		// merge
		IsingEnergy(MVM_new, lastSpins);
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
	static AHC ahc_instance(J_matrix);
	ahc_instance.ahc_solver(x_init);

    bestSpinsOut:for (int i = 0; i < N; i++) {
        bestSpinsOut[i] = ahc_instance.bestSpins[i];
    }
}
