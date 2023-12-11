#include "AHC_low.h"
#include "iostream"
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
using std::cout;
using std::endl;




AHC::AHC(){
// Partition local variables to improve bandwidth
#pragma HLS ARRAY_PARTITION variable=this->J complete dim=1
#pragma HLS ARRAY_PARTITION variable=this->x complete
#pragma HLS ARRAY_PARTITION variable=this->MVM_out complete
#pragma HLS ARRAY_PARTITION variable=this->e complete

	// Initialize the AHC solver
	this->dt   = 0.01;
	this->r    = 0.98;
	this->beta = 1;
	this->coupling_strength = 0.20;
	this->mu = 1.0;
	this->num_time_steps = 1;
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
	cout << endl;
	cout << "INIT SPINS ";

	#pragma HLS INLINE off
	setSpins_loop:
	for(int i=0; i<N; i++){
		#pragma HLS PIPELINE
		if(this->x[i] > 0){
			this->lastSpins[i] = 1;
			cout << "1, ";
		}
		else if(this->x[i] < 0){
			this->lastSpins[i] = -1;
			cout << "-1, ";

		}
		else{
			this->lastSpins[i] = 0;
			cout << "0, ";
		}
	}
	cout << endl;
}

// Matrix vector product
void AHC::matmul()
{
	#pragma HLS INLINE off
	// Matrix vector product
	// MVM = (J).dot(np.sign(x))
	for (int row = 0; row < N; row++) {
		cout << "MVM[row] is " << MVM_out[row] << " ";
		#pragma HLS PIPELINE
		this->MVM_out[row] = 0.0;
		cout << endl;
	}

	// using column method MVM = \sum_j J[:][j] * x[j]
	MVM_outer:
	for(int i=0; i<N; i++){
		// for each element in x

		#pragma HLS PIPELINE
		MVM_inner:
		for(int j = 0; j < N; j++){
			// multiply with each element on i th column of J

			if(this->x[j] == 1){
				this->MVM_out[i] += this->J[i][j];
			}
			else if(this->x[j] == -1){
				this->MVM_out[i] += -(this->J[i][j]);
			}

			else {
				this->MVM_out[i] += 0;

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
			energy += -((this->MVM_out[i])) >> 1;
		}
		else if(this->lastSpins[i]==-1){
			energy += (this->MVM_out[i]) >> 1;
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
	
	update_spin_and_error:
	for(int i=0; i<N; i++){
		#pragma HLS PIPELINE

		
		// Update spin vector
		this->x[i] += dt * (coupling_strength * this->MVM_out[i])*this->e[i];
		xx[i] = (this->x[i] * this->x[i]);
		this->x[i] += -dt * this->x[i] * ((data_type_x(0.02)) + mu*xx);
		this->de[i] = dt*(-(this->beta) * this->e[i] * (xx[i] - target_a));
        this->e[i] += de[i];
		
		/////////////////////////// option 2 using shift ////////////////////////////
		// data_type_x this_x_tmp;
		// // Update spin vector
		// // this->x[i] += dt * (coupling_strength * this->MVM_out[i])*this->e[i];
		// this_x_tmp = (this->MVM_out[i] >> 6) + (this->MVM_out[i] >> 8) + (this->MVM_out[i] >> 11);
		// this->x[i] += ((this_x_tmp >> 7) + (this_x_tmp >> 9)) * (this->e[i]);

		// // xx[i] = (this->x[i] >> 4);
		// xx = (this->x[i] >> 4);

		// // this->x[i] += -dt * this->x[i] * ((data_type_x(0.02)) + mu*this->xx[i]);
		// data_type_x this_x_tmp_2;
		// this_x_tmp_2 = (data_type_x(0.02)) + mu*xx;
		// this->x[i] += -((this_x_tmp_2 >> 7) + (this_x_tmp_2 >> 9)) * this->x[i];

		// // this->de[i] = dt*(-beta * this->e[i] * (this->xx[i] - target_a));
		// data_type_x this_tmp_de;
		// // this_tmp_de = -(((xx[i] - target_a) >> 4) + ((xx[i] - target_a) >> 5) + ((xx[i] - target_a) >> 7));
		// this_tmp_de = -(((xx - target_a) >> 4) + ((xx - target_a) >> 5) + ((xx - target_a) >> 7));

		// de = ((this_tmp_de >> 7) + (this_tmp_de >> 9)) * (this->e[i]);
		// this->e[i] += de;
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

	cout << " XINIT ";
	initialize_vectors:for(int i=0; i<N; i++){
		this->e[i] = 1.0;
		this->x[i] = x_init[i];
		cout << this->x[i] << ", ";
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
	for(int time_step=0; time_step < num_time_steps; time_step++){
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
	data_type_e bestEnergy;
	
	// Partition the first dim
	#pragma HLS ARRAY_PARTITION variable=J_in dim=1 complete
	#pragma HLS ARRAY_PARTITION variable=x_in dim=1 complete
	#pragma HLS ARRAY_PARTITION variable=bestSpinsOut dim=1 complete

	bit32_t input_l;
	bit32_t output;

	static AHC ahc_instance;

	// read J matrix
	for (int i = 0; i < N; i++) {
		#pragma HLS pipeline off
		for (int j = 0; j < N; j++) {	
			#pragma HLS pipeline
			input_l = strm_in.read();
			J_in[i][j] = input_l(MAX_WIDTH-1,0);
		}
	}
	ahc_instance.updateJ(J_in);

	// run 100 sets of X
	for (int x_iter=0; x_iter<20; x_iter++){
		#pragma HLS pipeline off
		// read x_init
		for (int i = 0; i < N; i++) {
			#pragma HLS pipeline
			input_l = strm_in.read();
			x_in[i] = input_l(MAX_WIDTH-1,0);
		}
		ahc_instance.ahc_solver(x_in);
	}
	
	// return the best energy
	bestEnergy = ahc_instance.bestEnergySpins(bestSpinsOut);
	output(MAX_WIDTH-1,0) = bestEnergy;
	strm_out.write(output);

	// write out the result
	for (int i = 0; i < N; i++) {
		output(MAX_WIDTH-1,0) = bestSpinsOut[i];
		strm_out.write(output);
	}
}
