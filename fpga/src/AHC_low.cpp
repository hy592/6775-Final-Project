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
#pragma HLS PIPELINE
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
            //cout<< this->J[i][j] << " ";
        }
        //cout<< endl;
    }

	/*
	// Initialize the spins, J matrix and MVM output
	#ifdef DEBUG
		cout<< "INIT SPINS HLS SIDE" << std::endl;
			for(int i = 0; i < N; i++){
				cout<< this->x[i]<< " ";
			}
	cout << endl;
	#endif
	*/
}

// returns the square of the input val
data_type_x square_single(data_type_x val) {
	data_type_x tmp_xx;
	#pragma HLS BIND_OP variable=tmp_xx op=mul impl=fabric
	//cout << "VAL " << val <<endl;
	tmp_xx = val * val;
	//cout << "SQVAL " << tmp_xx <<endl;
	return(tmp_xx);
}

/*
data_type_x square_single(data_type_x val) {
	data_type_x tmp_xx;
	#pragma HLS BIND_OP variable=tmp_xx op=mul impl=fabric
	//cout << "VAL " << val <<endl;
	tmp_xx = val * val;
	//cout << "SQVAL " << tmp_xx <<endl;
	return(tmp_xx);
}
*/

// Squares the input vector
void AHC::square(){
	#pragma HLS INLINE off
	#pragma HLS PIPELINE
	// Element wise square
	square_loop_1:for(int i=0;i<N/2;i++){
		#pragma HLS UNROLL FACTOR=8
		//tmp_xx = x[i][j] * x[i][j];
		this->xx[i]= square_single(this->x[i]);
	}
	square_loop_2:for(int i=N/2;i<N;i++){
		#pragma HLS UNROLL FACTOR=8
		//tmp_xx = x[i][j] * x[i][j];
		this->xx[i] = square_single(this->x[i]);
	}
}

// setSpins
// This function update the lastSpins based on the current spin values
// This is used to determine the sign of the spins in the MVM
void AHC::setSpins(){
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

#ifdef LAST
// Matrix vector product
void AHC::matmul()
{
	#pragma HLS INLINE off
	#pragma HLS PIPELINE
	#pragma HLS LATENCY min=4 max=5
	// Matrix vector product
	MVM_outer:for(int i=0;i<N;i++){
		this->MVM_out[i]=0.0;
		MVM_inner:for(int j=0;j<N;j++){
			#pragma HLS UNROLL factor = 8
			if(this->lastSpins[j]==1){
				//this->MVM_out[i] += data_type_x(-(this->J[i][j]));
				//cout << this->J[i][j] << " ";
				this->MVM_out[i] += this->J[i][j];
				//this->lastSpins[j] = 1;
				//(-(this->J[i][j]));
			}
			else if(this->lastSpins[j]==-1){
				//cout <<"LESS THAN ZERO     ";
				//cout << this->x[i]<< endl;
				//cout <<  -(this->J[i][j]) << " ";
				this->MVM_out[i] -= (this->J[i][j]);
				//this->lastSpins[i] = -1;
			}
			else{
				// cout << "hit a 0   ";
				// cout << this->x[i] << endl;
				this->MVM_out[i] += 0;
				// this->lastSpins[i] = 0;

			}
		}
	}
	/*
	#ifdef DEBUG
		cout << "MVM OUTPUT" << std::endl;
			for(int k = 0; k < N; k++){
				std::cout << this-> MVM_out[k] << " ";
		}
		cout <<endl;
		cout << "SPINS SIGN" << std::endl;
		for(int z = 0; z < N; z++){
			std::cout << this-> lastSpins[z] << " ";
		}
		cout <<endl;
		cout << "POS 48 " << MVM_out[48] << endl;
	#endif
	*/
}
#endif

#ifdef curr
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
				//this->MVM_out[i] += data_type_x(-(this->J[i][j]));
				//cout << this->J[i][j] << " ";
				this->MVM_out[i] += this->J[i][j];

				//this->lastSpins[j] = 1;
				//(-(this->J[i][j]));
			}
			else if(this->x[j]==-1){
				//cout <<"LESS THAN ZERO     ";
				//cout << this->x[i]<< endl;
				//cout <<  -(this->J[i][j]) << " ";
				this->MVM_out[i] -= (this->J[i][j]);
				//this->lastSpins[i] = -1;
			}
			else{
				//cout << "hit a 0   ";
					//	cout << this->x[i] << endl;
				this->MVM_out[i] += 0;
				// this->lastSpins[i] = 0;
			}
		}
	}
}
#endif


// Calculates the Ising energy
data_type_e AHC::IsingEnergy(){
#pragma HLS PIPELINE
	data_type_e energy = 0.0;
	for(int i = 0; i < N; i++){
		if(this->lastSpins[i]==1){
			data_type_e temp = -(this->MVM_out[i])/2;
			energy += temp;
			/*
			#ifdef DEBUG
				cout << temp<< endl;
				cout << "1, ";
			#endif
			*/
		}
		else if(this->lastSpins[i]==-1){
			data_type_e temp = this->MVM_out[i]/2;
			energy += temp;
			/*
			#ifdef DEBUG
				cout << temp<< endl;
				cout << "-1, ";
			#endif
			*/
		}
		else{
			//cout << "hit a 0   ";
			//cout << this->x[i] << endl;
			//cout << "0, ";
		}
		//energy += -1/2*(this->MVM_out[i]*this->x[i]);
	}
	//cout << endl;
	return energy;
}


/*
void AHC::update(){
#pragma HLS INLINE off
#pragma HLS PIPELINE

#pragma HLS LATENCY min=4 max=10
	update_spin:
	for(int i=0;i<N;i++){
		#pragma HLS UNROLL factor=8
		// Update spin vector
		this->dx[i] = dt * this->x[i] * ((r-1) - mu*xx[i]);
		this->dx2[i] = dt * (coupling_strength * this->MVM_out[i] * e[i]);
		this->de[i] = (-beta * this->e[i] * (this->xx[i] - target_a));
	}

	update_spin_and_error:
	for(int i=0;i<N;i++){
		#pragma HLS UNROLL factor=8
		// Update spin vector
		//Update error vector

		//de[i]

			this->e[i] +=	dt * this->de[i];
			this->x[i] += this->dx[i] + this->dx2[i];

		//e[i] += de[i];
	}
}
*/

// Update the spins and error vectors
void AHC::update(){
	#pragma HLS INLINE off
	#pragma HLS PIPELINE
	#pragma HLS LATENCY min=4 max=10

	update_spin_and_error:for(int i=0;i<N;i++){
		#pragma HLS UNROLL factor=8
		// Update spin vector
		//cout << "x squared is " << this->xx[i] << endl;
		this->xx[i] = square_single(this->x[i]);

		this->x[i] += dt * (coupling_strength * this->MVM_out[i])*this->e[i];
		this->xx[i] = this->x[i] >> 4;
		//square_single(this->x[i]);

		this->x[i] += -dt * this->x[i] * ((data_type_x(0.02)) + mu*this->xx[i]);


		//this->dx[i] = -dt * this->x[i] * (0.02 + mu*this->xx[i]);

		//cout << "Dx1 " << dx[i] << endl;

		//cout << "Dx2 " << dx2[i] << endl;
		//this->x[i] += this->dx[i]+ this->dx2[i];
				// * e[i]);
		//this->xx[i] = this->x[i] >> 3;

		//std::cout << "XX " << this->xx[i] << endl;
		this->de[i] = dt*(-beta * this->e[i] * (this->xx[i] - target_a));
		//std::cout << "de " << this->de[i] << endl;
		// Update error vector
		//this->e[i] +=	dt * this->de[i];
		//this->x[i] += this->dx[i] + this->dx2[i];
        this->e[i] += de[i];
	}

	/*
	#ifndef DEBUG
		cout << "X1" << endl;

		for(int k=0; k < N; k++){
			cout << this->dx[k] << " ";
		}
		cout << endl;
		cout << "X2" << endl;

			for(int z=0; z < N; z++){
				cout << dx2[z] << " ";
			}
		cout << endl;
	#endif
	*/
}

void AHC::reset(){
	#pragma HLS INLINE off
	#pragma HLS PIPELINE
	// Reset MVM
	reset_MVM:for(int i=0;i<N;i++){
		this->MVM_out[i] = 0.0;
	}
}

#ifdef DEBUG
void AHC::writeDebug(int index){
	for(int i =0; i < N; i++){
		x1_debug[index][i] = this->dx[i];
		x2_debug[index][i] = this->dx2[i];

	}
}
#endif

void AHC::updateSpins(data_type_x x_init[N], ap_fixed<MAX_WIDTH, 2> coupling_strength_new, ap_fixed<MAX_WIDTH,2> new_mu){
	#ifdef DEBUG
		cout << "SPINS INIT" << endl;
	#endif
	initialize_vectors:for(int i=0;i<N;i++){
		this->e[i] = 1.0;
		this->x[i] = x_init[i];
		this->MVM_out[i] = 0.0;
		this->lastSpins[i] = x_init[i];
		this->de[i] = 0;
		this->coupling_strength = coupling_strength_new;
		this->mu = new_mu;
		#ifdef DEBUG
			cout << this->x[i] << " ";
		#endif
	}

	#ifdef DEBUG
	cout << endl;
	#endif

}

void AHC::ahc_solver(){
	#pragma HLS ARRAY_PARTITION variable=xx dim=0 complete
	#pragma HLS ARRAY_PARTITION variable=dx dim=0 complete
	#pragma HLS ARRAY_PARTITION variable=dx2 dim=0 complete
	#pragma HLS ARRAY_PARTITION variable=de dim=0 complete

	data_type_e energy = 0.0;
	//cout < "NUMMM STEPS " << this->num_time_steps << endl;

	setSpins();
	matmul();
	iterations:
	for(int time_step=0;time_step<this->num_time_steps;time_step++){

		//square();
		update();
		//writeDebug(time_step);
		setSpins();
		matmul();

		/*
		#ifdef DEBUG
			cout << "SPINS" << endl;
			for(int z = 0; z < N; z++){
				cout << this->x[z] << " ";
			}
			cout << endl;
		#endif
		*/

		setSpins();
		energy = IsingEnergy();
		//	cout << "ISING ENERGY " << energy << endl;
		if(energy < this->bestEnergy){
			bestEnergy = energy;
			for(int k = 0; k < N; k++){
				bestSpins[k] = this->lastSpins[k];
			}
		}
	}


	// #ifdef DEBUG
	// 	cout << "FINAL SPINS" << std::endl;
	// 	for(int i = 0; i < N; i++){
	// 		std::cout << this->bestSpins[i]<< " ,";
	// 	}
	// 	cout << endl;
	// 	cout << "best energy " << bestEnergy << endl;
	// #endif


	/*
	#ifdef DEBUG
		std::ofstream outFileX1;
		std::ofstream outFileX2;

		outFileX1.open("data_x1.txt");
		outFileX2.open("data_x2.txt");

		// Loop over each pair
		for(int pairIndex = 0; pairIndex < 150; pairIndex++) {
			// Open a new file for each pair

			// Loop over the arrays and write the data to the file
			for(int i = 0; i < N; i++) {
				outFileX1 << x1_debug[pairIndex][i] << " ";  // Separate values by space
				outFileX2 << x2_debug[pairIndex][i] << " ";  // Separate values by space

			}
			outFileX1 << "\n";
			outFileX2 << "\n";

		}
		outFileX1.close();  // Don't forget to close the file
		outFileX2.close();  // Don't forget to close the file
	#endif
	*/
}

void ahc_top(
	data_type_J J_matrix[N][N], 
	data_type_x x_init[N], 
	spin_sign bestSpinsOut[N]
){
	#pragma HLS INTERFACE bram port=bestSpinsOut
	// static AHC ahc_instance = AHC(x_init, J_matrix);
	static AHC ahc_instance(x_init, J_matrix);
	ahc_instance.ahc_solver();

    for (int i = 0; i < N; i++) {
        bestSpinsOut[i] = ahc_instance.bestSpins[i];
    }
}
