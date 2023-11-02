#include "ahc_header.h"

AHC::AHC(){
#pragma HLS PIPELINE

	dt   = 0.01;
	r    = 1.0;
	beta = 0.05;
	coupling_strength = 0.07;
	mu = 1.0;
	num_time_steps = 4;
	target_a_baseline = 0.2;
	target_a = target_a_baseline;

	initialize_vectors:for(int i=0;i<N;i++){
		e[i] = 1.0;
		x[i] = 0.1 * ((random_integer()/2147483647.0) - 0.5);
		MVM[i] = 0.0;
	}

#ifdef __DEBUG__

	printf("\n Using parameters:\n");
	printf("\n dt = %f", dt.to_float());
	printf("\n r = %f", r.to_float());
	printf("\n beta = %f", beta.to_float());
	printf("\n coupling strength = %f", coupling_strength.to_float());
	printf("\n mu = %f", mu.to_float());
	printf("\n target_a_baseline = %f", target_a_baseline.to_float());

	printf("\n num_time_steps = %d", num_time_steps);

	debug(-1);
#endif
}


#ifdef __DEBUG__

void AHC::debug(int time_step){
#ifndef __USE_FLOAT__
	printf("\n\n");
	printf("Time step = %d\n", time_step);
	printf("x = \n");
	print_x:for(int i=0;i<N;i++){
		printf("%f, ", x[i].to_float());
	}

	printf("\n");

	printf("e = \n");
	print_e:for(int i=0;i<N;i++){
		printf("%f, ", e[i].to_float());
	}

	printf("\n");

	printf("MVM = \n");
	print_MVM:for(int i=0;i<N;i++){
		printf("%f, ", MVM[i].to_float());
	}
#endif

#ifdef __USE_FLOAT__

	printf("x = \n");
	print_x:for(int i=0;i<N;i++){
		printf("%f, ", x[i]);
	}

	printf("\n");

	printf("e = \n");
	print_e:for(int i=0;i<N;i++){
		printf("%f, ", e[i]);
	}

	printf("\n");

	printf("MVM = \n");
	print_MVM:for(int i=0;i<N;i++){
		printf("%f, ", MVM[i]);
	}
#endif

	printf("\n");
}

#endif

void AHC::square(){
#pragma HLS INLINE off
#pragma HLS PIPELINE
	// Element wise square
	square_loop:for(int i=0;i<N;i++){
		xx[i] = x[i] * x[i];
	}
}

void AHC::matmul(){
#pragma HLS INLINE off
#pragma HLS PIPELINE
	// Matrix vector product
	MVM_outer:for(int j=0;j<N;j++){
		MVM_inner:for(int i=0;i<N;i++){
			#pragma HLS UNROLL
			MVM[i] += J[i][j] * sign_bit(x[j]);
			printf("%d\n", sign_bit(x[j]));
		}
	}
}

void AHC::update(){
#pragma HLS INLINE off
#pragma HLS PIPELINE

	update_spin_and_error:for(int i=0;i<N;i++){
		#pragma HLS UNROLL
		// Update spin vector
		dx[i] = dt * x[i] * ((r-1) - mu*xx[i]);
		dx2[i] = dt * (coupling_strength * MVM[i] * e[i]);
		x[i] += dx[i] + dx2[i];
		//Update error vector
		de[i] = dt * (-beta * e[i] * (xx[i] - target_a));
		e[i] += de[i];
	}

	// Include target amplitude modulation here
	//#########
	//#########
}

void AHC::reset(){
#pragma HLS INLINE off
#pragma HLS PIPELINE
	// Reset MVM
	reset_MVM:for(int i=0;i<N;i++){
		MVM[i] = 0.0;
	}
}


void AHC::ahc_solver(){
	#pragma HLS ARRAY_PARTITION variable=xx dim=1 complete
	#pragma HLS ARRAY_PARTITION variable=dx dim=1 complete
	#pragma HLS ARRAY_PARTITION variable=dx2 dim=1 complete
	#pragma HLS ARRAY_PARTITION variable=de dim=1 complete

	iterations:for(int time_step=0;time_step<num_time_steps;time_step++){
		#pragma HLS PIPELINE off
		square();
		matmul();
		update();

#ifdef __DEBUG__
		debug(time_step);
#endif

		reset();
	}
}

extern "C" {
void ahc(float J_matrix[N*N], 
	float out_vector[N], 
	float value, 
	short option){
	
	static AHC ahc_instance;
	#pragma HLS ARRAY_PARTITION variable=ahc_instance.J dim=0 complete
	#pragma HLS ARRAY_PARTITION variable=ahc_instance.x dim=1 complete
	#pragma HLS ARRAY_PARTITION variable=ahc_instance.MVM dim=1 complete
	#pragma HLS ARRAY_PARTITION variable=ahc_instance.e dim=1 complete


	if(option == 1){
		ahc_instance.ahc_solver();
		// Write-out
		write_out:for(int i=0;i<N;i++){
		#pragma HLS UNROLL
			out_vector[i] = ahc_instance.x[i].to_float();
		}
	}

	else if(option==2){
		debug_write_out:for(int i=0;i<N;i++){
		#pragma HLS UNROLL
			out_vector[i] = float(i);
		}
	}

	else if(option==3){
		write_J_outer: for(int i=0;i<N;i++){
			write_J_inner: for(int j=0;j<N;j++){
				ahc_instance.J[i][j] = J_matrix[i*N + j];
			}
		}

#ifdef __DEBUG__
		printf("\n J matrix (fixed point) = \n");
		for(int i=0;i<N;i++){
				for(int j=0;j<N;j++){
					printf("%f ", ahc_instance.J[i][j].to_float());
				}
				printf("\n");
		}
		printf("\n\n");
#endif
	}
}
}
