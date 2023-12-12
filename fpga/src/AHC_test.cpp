#include <cstdlib>  // Required for rand() function
#include <fstream>
#include <iostream>
#include <random>
#include <stdlib.h>

#include "AHC_low.h"
#include "ap_fixed.h"

#define USER "zs343/ece6775"
#define PROJECT_DIR_NAME "final_project"
#define DIR "/home/" USER "/" PROJECT_DIR_NAME "/sample_MIMO"
#define JMATRICES "/home/" USER "/" PROJECT_DIR_NAME "/fixed_pt_data" 
#define IDEALS DIR "/ideal_sols_16_nr16_16QAM" 

using namespace std;

const int matrix_size = 65;

void top_func_test( data_type_J J[N][N], data_type_x init_spins[20][N], int problemNum )
{
  spin_sign bestSpins[N];
  data_type_e bestEnergy;

  AHC ahc_instance;  // AHC instance
  ahc_instance.updateJ(J);    // load problem set
  // ahc_instance.ahc_solver();
  for ( int i = 0; i < 20; i++ ) {
    // ahc_instance.updateSpins( init_spins[i], 0.02,1.3 );
    ahc_instance.ahc_solver(init_spins[i]);
  }

  bestEnergy = ahc_instance.bestEnergySpins(bestSpins);
  // spin_sign spin_guesses[N];

  std::cout << "BEST ENERGY = " << bestEnergy << endl;
  
  std::cout << "Spin = " << std::endl;
  for (int i = 0; i < N; ++i)
  {
    std::cout << bestSpins[i] << std::endl;
  }

  /* ifstream ideal_solutions(IDEALS + "/DI_MIMO_mmse_sol_" + std::to_string(problemNum) + ".txt");

  spin_sign ideal_spins[N]; */
  
  

  // outSpins.close();  // Don't forget to close the file
}
int main()
{
  data_type_x x_init_arrays[20][matrix_size];

  data_type_J matrix[matrix_size][matrix_size];

  // Seed the random number generator
  srand( 1234 );

  data_type_x random_array[matrix_size] = { 3.58964059e-05, 1.63794645e-04, 1.48891121e-05, 4.44594756e-04, 8.65550405e-05, 4.03401915e-04, -3.62525296e-04, -3.60723653e-04, 3.07391289e-04, -1.02323163e-04, -3.34645803e-04, 4.27508580e-04, -1.52234140e-04, 2.50812103e-04, 2.25997985e-04, 3.83306091e-04, 1.23672207e-04, 2.50942434e-04, -1.51101658e-04, -2.30072108e-04, 3.95886218e-04, -7.19088101e-05, 4.64840047e-04, 1.63441498e-04, 1.21695720e-04, -3.85254027e-04, 4.49489259e-04, -5.00878665e-05, 7.83896144e-05, -9.18631972e-05, -2.62973020e-04, 4.03379521e-04, 7.36794867e-05, -4.97129673e-04, 1.17144914e-04, -1.73355098e-04, 2.70581023e-05, 3.85942099e-04, -1.42730240e-04, 4.08535151e-04, 1.23360116e-04, -4.84178757e-04, 4.29437234e-04, 1.90896918e-04, 4.97322850e-04, -3.27659492e-04, -3.62864250e-04, 4.32595463e-04, 1.96818161e-04, -4.33999827e-04, 2.55463053e-04, 2.53876188e-04, 4.23024536e-04, 2.11524759e-04, -3.75729038e-04, -4.80119866e-04, -4.73789013e-04, -4.71693512e-04, -2.53788932e-04, 3.60027949e-04, 3.88310643e-05, 5.28219787e-05, 3.42030892e-04, -3.75826685e-04, -2.20816321e-04 };

  // cout << "INIT SPINS -----------------" << endl;

  ifstream x_init_file( "init_spins_big.txt" );
  if ( !x_init_file.is_open() ) {
    cout << "Error opening x_init file!" << endl;
    return 1;
  }

  // Define a 2D array to store all 20 x_init arrays

  for ( int array_idx = 0; array_idx < 20; array_idx++ ) {
    for ( int i = 0; i < matrix_size; i++ ) {
      if ( !( x_init_file >> x_init_arrays[array_idx][i] ) ) {
        cout << "Error reading x_init file!" << endl;
        return 1;
      }
    }
  }

  x_init_file.close();

  // int numProblems = 10;

  int startNum = 0;

  //  Do something with the matrix values here...

  for ( int k = startNum; k < startNum + numProblems; k++ ) {
    // ifstream input_file1( JMATRICES + (string)"/DI_MIMO_J_Conv_11_29_16_2_" + std::to_string( k ) + (string)".txt" );
    ifstream input_file1((string)"fixed_pt_data" + (string)"/DI_MIMO_J_Conv_11_29_16_2_" + to_string(k) + (string)".txt");
    // cout << JMATRICES + (string)"/DI_MIMO_J_Conv_11_29_16_2_" + to_string( k ) + (string)".txt" << endl;
    cout << (string)"fixed_pt_data" + (string)"/DI_MIMO_J_Conv_11_29_16_2_" + to_string(k) + (string)".txt" << endl;

#ifdef DEBUG
    cout << "INIT SPINS" << endl;
    for ( int z = 0; z < matrix_size; z++ )
      cout << random_array[z] << " ";
    cout << endl;
#endif DEBUG

    for ( int i = 0; i < matrix_size; i++ ) {
      for ( int j = 0; j < matrix_size; j++ ) {
        int hex_value;
        input_file1 >> hex >> hex_value;
        matrix[i][j] = static_cast<float>( hex_value ) / 4096.0f;  // 12_2
        // printf("%X\n", hex_value);
        // cout << matrix[i][j] << endl;
      }
    }
    top_func_test( matrix, x_init_arrays, k );
    input_file1.close();
  }

  return 0;
}
