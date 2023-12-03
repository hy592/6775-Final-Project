#include <stdio.h>
#include "AHC_low.h"

using namespace std;

int main(){
	printf("Enter Simulation\n");
	/*  spin_sign* ahc_top(
			data_type_J J_matrix[N][N], 
			data_type_x x_init[N], 
			spin_sign bestSpinsOut[N]
		);
	*/
	data_type_J J_matrix[N][N];
	// initialize J_matrix
	for (int i=0; i<N; i++){
		for (int j=0; j<N; j++){
			J_matrix[i][j] = 0;
		}
	}

	data_type_x x_init[N];
	// initialize x_init
	for (int i=0; i<N; i++){
		x_init[i] = 0;
	}

	spin_sign bestSpinsOut[N];
	// initialize bestSpinsOut
	for (int i=0; i<N; i++){
		bestSpinsOut[i] = 0;
	}

	data_type_e bestEnergy = 10;

	// ahc_top(J_matrix, x_init, bestSpinsOut);
	static AHC ahc_instance;
	ahc_instance.updateJ(J_matrix);
	ahc_instance.ahc_solver(x_init);
	
	bestEnergy = ahc_instance.bestEnergySpins(bestSpinsOut);
	
	// print ahc_instance_bestSpins
	for (int i = 0; i < N; i++)
	{
		std::cout << bestSpinsOut[i] << "  ";
	}

	printf("Finish Simulation\n");
	return 0;
}
