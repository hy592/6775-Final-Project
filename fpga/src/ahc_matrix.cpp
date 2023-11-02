#include "ahc_matrix_header.h"

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
        for(int j=0;j<N;j++){
            e[i][j] = 1.0;
            x[i][j] = 0.1 * ((random_integer()/2147483647.0) - 0.5);
            MMM[i][j] = 0.0;
        }
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
        for(int j=0;j<N;j++){
		    printf("%f, ", x[i][j].to_float());
        }
	}

	printf("\n");

	printf("e = \n");
	print_e:for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
		    printf("%f, ", e[i][j].to_float());
        }
	}

	printf("\n");

	printf("MMM = \n");
	print_MMM:for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
		    printf("%f, ", MMM[i][j].to_float());
        }
	}
#endif

#ifdef __USE_FLOAT__

	printf("x = \n");
	print_x:for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
		    printf("%f, ", x[i][j]);
        }
	}

	printf("\n");

	printf("e = \n");
	print_e:for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
		    printf("%f, ", e[i][j]);
        }
	}

	printf("\n");

	printf("MMM = \n");
	print_MMM:for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
		    printf("%f, ", MMM[i][j]);
        }
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
        for(int j=0;j<N;j++){
		    xx[i][j] = x[i][j] * x[i][j];
        }
	}
}

void AHC::matmul(){
#pragma HLS INLINE off
//#pragma HLS PIPELINE
	// Matrix vector product
	MMM_outer:for(int i=0;i<N;i++){
        #pragma HLS PIPELINE
		MMM_inner:for(int j=0;j<N;j++){
            MMM_k:for(int k=0;k<N;k++){
                MMM[i][j] += J[i][k] * sign_bit(x[k][j]);
            }
		}
	}
}

void AHC::update(){
#pragma HLS INLINE off
#pragma HLS PIPELINE

	update_spin_and_error:for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            // Update spin vector
            dx[i][j] = dt * x[i][j] * ((r-1) - mu*xx[i][j]);
            dx2[i][j] = dt * (coupling_strength * MMM[i][j] * e[i][j]);
            x[i][j] += dx[i][j] + dx2[i][j];
            //Update error vector
            de[i][j] = dt * (-beta * e[i][j] * (xx[i][j] - target_a));
            e[i][j] += de[i][j];
        }
	}

	// Include target amplitude modulation here
	//#########
	//#########
}

void AHC::reset(){
#pragma HLS INLINE off
#pragma HLS PIPELINE
	// Reset MMM
	reset_MVM:for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
		    MMM[i][j] = 0.0;
        }
	}
}


void AHC::ahc_solver(){
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
	#pragma HLS ARRAY_PARTITION variable=ahc_instance.x dim=0 complete
	#pragma HLS ARRAY_PARTITION variable=ahc_instance.MMM dim=0 complete
	#pragma HLS ARRAY_PARTITION variable=ahc_instance.e dim=0 complete
	#pragma HLS ARRAY_PARTITION variable=ahc_instance.xx dim=0 complete
	#pragma HLS ARRAY_PARTITION variable=ahc_instance.dx dim=0 complete
	#pragma HLS ARRAY_PARTITION variable=ahc_instance.dx2 dim=0 complete
	#pragma HLS ARRAY_PARTITION variable=ahc_instance.de dim=0 complete


	if(option == 1){
		ahc_instance.ahc_solver();
        // Replace with finding best anneal
        int s = N-1;
		// Write-out
		write_out:for(int i=0;i<N;i++){
		#pragma HLS UNROLL
			out_vector[i] = ahc_instance.x[s][i].to_float();
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
