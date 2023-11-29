#include "AHC_new.h"

AHC::AHC(const data_type_x x_init[N], const data_type_J J_init[N][N]){
// Partition local variables to improve bandwidth
#pragma HLS ARRAY_PARTITION variable=this->J dim=2 complete
#pragma HLS ARRAY_PARTITION variable=this->x dim=1 complete
#pragma HLS ARRAY_PARTITION variable=this->MVM_out dim=1 complete
#pragma HLS ARRAY_PARTITION variable=this->e dim=1 complete
#pragma HLS ARRAY_PARTITION variable=this->bestSpins dim=1 complete
#pragma HLS ARRAY_PARTITION variable=this->lastSpins dim=1 complete
	// Initialize the AHC solver
	// dt   = 0.01;
	// r    = 0.98;
	// beta = 0.1;
	// coupling_strength = 0.020;
	// mu = 1.0;
	// num_time_steps = 200;
	// target_a_baseline = 0.2;
	// target_a = target_a_baseline;
	target_a = 0.2;

	// Initialize the spins, J matrix and MVM output
	Initialize:for(int i=0;i<N;i++){
		#pragma HLS pipeline off
		#pragma HLS unroll factor=8
		// #pragma HLS pipeline 
		this->e[i] = 1.0;
		this->x[i] = x_init[i];
		this->MVM_out[i] = 0.0;
		this->lastSpins[i] = x_init[i];
	}

	Initialize_J:for(int j=0;j<N;j++){
		// #pragma HLS pipeline 
		Initialize_J_inn:for(int i=0;i<N;i++){
			#pragma HLS pipeline off
			#pragma HLS unroll factor=8
            this->J[i][j] = J_init[i][j];
        }
    }
}

void AHC::ahc_solver(){
	#pragma HLS ARRAY_PARTITION variable=xx dim=0 complete
	#pragma HLS ARRAY_PARTITION variable=dx dim=0 complete
	#pragma HLS ARRAY_PARTITION variable=dx2 dim=0 complete
	#pragma HLS ARRAY_PARTITION variable=de dim=0 complete

	Initialize_ahc_solver:for(int i=0; i<N; i++){
		#pragma HLS pipeline off
		// #pragma HLS pipeline 
		#pragma HLS unroll factor=8
		setSpins(i);	// initialize vectors
		matmul(i);
	}

	// ahc_solver_time_step:for(int time_step=0; time_step < this->num_time_steps; time_step++){
	ahc_solver_time_step:for(int time_step=0; time_step < num_time_step; time_step++){
		// #pragma HLS pipeline off
		// #pragma HLS pipeline
		// #pragma HLS latency min=4 max=5
		
		data_type_e energy = 0.0;
		ahc_solver_update:for(int i=0; i<N; i++){
			#pragma HLS pipeline
            update(i);
            setSpins(i);
		}
		ahc_solver_matmul:for(int i=0; i<N; i++){
			#pragma HLS pipeline
			// #pragma HLS LATENCY min=4 max=5
            matmul(i);
            energy += IsingEnergy(i);
		}
		// ahc_solver_energy:for(int i=0; i<N; i++){
		// 	#pragma HLS pipeline
        //     // setSpins(i);
        // }

        // find bestEnergy for all iteration
		if(energy < this->bestEnergy){
			this->bestEnergy = energy;
			ahc_solver_best_energy:for(int k = 0; k < N; k++){
				#pragma HLS pipeline
				this->bestSpins[k] = this->lastSpins[k];
			}
		}
	}
}

// Update the spins and error vectors
void AHC::update(int i){
	#pragma HLS INLINE
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
	// this_x_tmp_2 = (data_type_x(0.02)) + mu*this->xx[i];
	this_x_tmp_2 = (data_type_x(0.02)) + this->xx[i];
	this->x[i] += -((this_x_tmp_2 >> 7) + (this_x_tmp_2 >> 9)) * this->x[i];

	// this->de[i] = dt*(-beta * this->e[i] * (this->xx[i] - target_a));
	data_type_x this_tmp_de;
	this_tmp_de = -(((this->xx[i] - target_a) >> 4) + ((this->xx[i] - target_a) >> 5) + ((this->xx[i] - target_a) >> 7));
	this->de[i] = ((this_tmp_de >> 7) + (this_tmp_de >> 9)) * (this->e[i]);
	this->e[i] += de[i];
}

// setSpins
// This function update the lastSpins based on the current spin values
// This is used to determine the sign of the spins in the MVM
void AHC::setSpins(int i){
	#pragma HLS INLINE
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

// Matrix vector product
void AHC::matmul(int i){
	#pragma HLS INLINE
	this->MVM_out[i]=0.0;
	matmul:for(int j=0;j<N;j++){
		#pragma HLS unroll factor=8
		data_type_x temp;
		if(this->x[j]==1){
			temp = this->J[i][j];
		}
		else if(this->x[j]==-1){
			temp = -(this->J[i][j]);
			//this->lastSpins[i] = -1;
		}
		else{
			temp = 0;
			// this->lastSpins[i] = 0;
		}
		this->MVM_out[i] += temp;
	}
}

// Calculates the Ising energy
data_type_e AHC::IsingEnergy(int i){
	#pragma HLS INLINE
	data_type_e temp;
	if(this->lastSpins[i]==1){
		temp = -((this->MVM_out[i]) >> 1);
	}
	else if(this->lastSpins[i]==-1){
		temp = (this->MVM_out[i] >> 1);
	}
	return temp;
}

void AHC::updateSpins(data_type_x x_init[N])
{
	#ifdef DEBUG
		cout << "SPINS INIT" << endl;
	#endif
	initialize_vectors:for(int i=0;i<N;i++){
		this->e[i] = 1.0;
		this->x[i] = x_init[i];
		this->MVM_out[i] = 0.0;
		this->lastSpins[i] = x_init[i];
		// this->de[i] = 0;
		// this->coupling_strength = coupling_strength_new;
		// this->mu = new_mu;
		#ifdef DEBUG
			cout << this->x[i] << " ";
		#endif
	}
	#ifdef DEBUG
	cout << endl;
	#endif
}

void ahc_top(
	const data_type_J J_matrix[N][N], 
	const data_type_x x_init[N], 
	spin_sign bestSpinsOut[N]
){
	// Partition the first dim
	// #pragma HLS ARRAY_PARTITION variable=J_matrix dim=1 complete
	// #pragma HLS ARRAY_PARTITION variable=x_init dim=1 complete

	#pragma HLS INTERFACE bram port=J_matrix
	#pragma HLS INTERFACE bram port=x_init
	#pragma HLS INTERFACE bram port=bestSpinsOut

	static AHC ahc_instance(x_init, J_matrix);
	ahc_instance.ahc_solver();

	#pragma HLS pipeline off
    top_chann:for (int i = 0; i < N; i++) {
		// #pragma HLS pipeline 
		#pragma HLS unroll factor=8
        bestSpinsOut[i] = ahc_instance.bestSpins[i];
    }
}
