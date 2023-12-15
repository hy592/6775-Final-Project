#include "AHC_low.h"
#include "iostream"
//#include <fstream>

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
	this->coupling_strength = 0.02;
	this->mu = 1;
	this->num_time_steps = 200;
	this->target_a_baseline = 0.2;
	this->target_a = target_a_baseline;
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
	for (int row = 0; row < N; row++) {
		#pragma HLS PIPELINE
		this->MVM_out[row] = 0.0;
	}

	// using column method MVM = \sum_j J[:][j] * x[j]
	MVM_outer:
	// for each element in x
	for(int j = 0; j < N; j++){
		#pragma HLS PIPELINE
		MVM_inner:
		for(int i=0; i<N; i++){
			// multiply with each element on i th column of J
			if(this->x[j] > 0){
				this->MVM_out[i] += this->J[i][j];
			}
			else if(this->x[j] < 0){
				this->MVM_out[i] += -(this->J[i][j]);
			}
			else {
				this->MVM_out[i] += 0;
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
	// #pragma HLS ARRAY_PARTITION variable=xx dim=0 complete
	// #pragma HLS ARRAY_PARTITION variable=de dim=0 complete

	data_type_x prevX;
	update_spin_and_error:
	for(int i=0; i<N; i++){
		#pragma HLS PIPELINE
		prevX = x[i];
		
		// Update spin vector
		this->x[i] += dt * (coupling_strength * this->MVM_out[i])*this->e[i];
		//if(prevX != x[i])
		//	spinFlips +=1;
		xx[i] = (this->x[i] << 2);
		//cout << -dt * this->x[i] * ((data_type_x(0.02)) + mu*xx[i]) << " ";
		this->x[i] += -dt * this->x[i] * ((data_type_x(0.02)) + mu*xx[i]); //0.2 is r-1 (1.2-1)==0.2
		this->de[i] = dt*(-(this->beta) * this->e[i] * (xx[i] - target_a));
		//cout << de[i] << " ";
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
	//cout << endl;

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
	Mat_Vec_Mal();

	TIME_STEP_LOOP:
	for(int time_step=0; time_step < num_time_steps; time_step++){
		update();
		IsingEnergy();
		setSpins();
		Mat_Vec_Mal();
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
	bit32_t output_energy;
	bit32_t output_spin;

	static AHC ahc_instance;

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

	// run 20 sets of X
	for (int x_iter=0; x_iter<20; x_iter++){
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
		ahc_instance.ahc_solver(x_in);
	}
	
	// return the best energy
	bestEnergy = ahc_instance.bestEnergySpins(bestSpinsOut);
	output_energy(MAX_WIDTH-1,0) = bestEnergy(MAX_WIDTH-1,0);
	strm_out.write(output_energy);

	// write out the result
	for (int i = 0; i < 8; i++) {
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
