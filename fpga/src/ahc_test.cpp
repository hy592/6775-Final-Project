#include "ahc_header.h"
#include <stdio.h>
#include <fstream>

int main(int argc, char **argv){
	float J[N*N];

	ifstream f;

    if(argc != 2){
        printf("Please include weight matrix file!\n");
        exit(3);
    }

    char *file_name = argv[1];
    float temp;

    printf("\n%s\n", file_name);

    f.open(file_name,ios::in);

	if(!f) {
		printf("Could not open weight file!\n");
		exit(1);
	}

    for(int i=0;i<N;i++){
		for(int j=0;j<N;j++){
			f >> temp;
			J[i*N+j] = temp;
		}
	}
    f.close();

    for(int i=0;i<N;i++){
    		for(int j=0;j<N;j++){
    			printf("%f ", J[i*N+j]);
    		}
    		printf("\n");
    	}

//	 for(int i=0;i<N;i++){
//	 	for(int j=0;j<N;j++){
//	 		if(abs(j-i) == 1){
//	 			J[i][j] = -1;
//	 			J[j][i] = -1;
//	 		}
//	 		else{
//	 			J[i][j] = 0.0;
//	 		}
//	 	}
//	 }

//	for(int i=0;i<N;i++){
//			for(int j=0;j<N;j++){
//				printf("%f ", J[i][j]);
//			}
//			printf("\n");
//	}

	printf("\n");

	float out_vector[N];
	ahc(J, out_vector, 0.0, 3);
	//run(J, out_vector, 0.0, 2);
	ahc(J, out_vector, 0.0, 1);

	printf("\nFinal out vector = \n");

#ifndef __USE_FLOAT__
	for(int i=0;i<N;i++){
		printf("%f, ", out_vector[i]);
	}
#endif
	return 0;
}
