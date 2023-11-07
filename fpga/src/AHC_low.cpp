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
// #pragma HLS PIPELINE
#pragma HLS ARRAY_PARTITION variable=this->x dim=0 complete
#pragma HLS ARRAY_PARTITION variable=this->J dim=0 complete
#pragma HLS ARRAY_PARTITION variable=this->MVM_out dim=0 complete
#pragma HLS ARRAY_PARTITION variable=this->e dim=0 complete

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
	for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            this->J[i][j] = J_init[i][j];
        }
    }
}

// returns the square of the input val
data_type_x square_single(data_type_x val) {
	data_type_x tmp_xx;
	#pragma HLS BIND_OP variable=tmp_xx op=mul impl=fabric
	tmp_xx = val * val;
	return(tmp_xx);
}

// Squares the input vector
void AHC::square(){
	#pragma HLS INLINE off
	#pragma HLS PIPELINE
	// Element wise square
	square_loop:for(int i=0;i<N;i++){
		#pragma HLS UNROLL FACTOR=8
		this->xx[i]= square_single(this->x[i]);
	}
}

// setSpins
// This function update the lastSpins based on the current spin values
// This is used to determine the sign of the spins in the MVM
void AHC::setSpins(){
	setSpins_loop:
	for(int i =0; i < N; i++){
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
	#pragma HLS PIPELINE
	#pragma HLS LATENCY min=4 max=5
	// Matrix vector product
	// MVM = (J).dot(np.sign(x))
	MVM_outer:for(int i=0;i<N;i++){
		this->MVM_out[i]=0.0;
		MVM_inner:for(int j=0;j<N;j++){
			#pragma HLS UNROLL factor = 8
			if(this->x[j]==1){
				this->MVM_out[i] += this->J[i][j];
			}
			else if(this->x[j]==-1){
				this->MVM_out[i] -= (this->J[i][j]);
				//this->lastSpins[i] = -1;
			}
			else{
				this->MVM_out[i] += 0;
				// this->lastSpins[i] = 0;
			}
		}
	}
}

// Calculates the Ising energy
data_type_e AHC::IsingEnergy(){
#pragma HLS PIPELINE
	data_type_e energy = 0.0;
	IsingEnergy_loop:
	for(int i = 0; i < N; i++){
		if(this->lastSpins[i]==1){
			data_type_e temp = -(this->MVM_out[i]) >> 1;
			energy += temp;
		}
		else if(this->lastSpins[i]==-1){
			data_type_e temp = this->MVM_out[i] >> 1;
			energy += temp;
		}
	}
	return energy;
}

// Update the spins and error vectors
void AHC::update(){
	#pragma HLS INLINE off
	#pragma HLS PIPELINE
	#pragma HLS LATENCY min=4 max=10

	update_spin_and_error:
	for(int i=0;i<N;i++){
		#pragma HLS UNROLL factor=8
		// Update spin vector
		this->x[i] += dt * (coupling_strength * this->MVM_out[i])*this->e[i];
		this->xx[i] = this->x[i] >> 4;
		this->x[i] += -dt * this->x[i] * ((data_type_x(0.02)) + mu*this->xx[i]);
		this->de[i] = dt*(-beta * this->e[i] * (this->xx[i] - target_a));
        this->e[i] += de[i];
	}
}

void AHC::reset(){
	#pragma HLS INLINE off
	#pragma HLS PIPELINE
	// Reset MVM
	reset_MVM:for(int i=0;i<N;i++){
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

	setSpins();
	matmul();
	iterations:
	for(int time_step=0;time_step<this->num_time_steps;time_step++){
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
	#pragma HLS INTERFACE bram port=bestSpinsOut
	static AHC ahc_instance(x_init, J_matrix);
	ahc_instance.ahc_solver();

    for (int i = 0; i < N; i++) {
        bestSpinsOut[i] = ahc_instance.bestSpins[i];
    }
}
