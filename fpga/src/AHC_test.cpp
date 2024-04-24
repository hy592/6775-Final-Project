#include <stdlib.h>

#include <cstdlib>  // Required for rand() function
#include <fstream>
#include <iostream>
#include <random>

//#include "AHC_low_6775.h"
#include "AHC_low.h"
//#include "AHC_low_old.h"
// #include "ap_fixed.h"

using namespace std;

std::string SNR;  // SNR value

// define a list of SNR values to solve
int SNR_list[] = {30}; // {0, 5, 10, 15, 20, 25, 30};
const int numSNR = sizeof(SNR_list)/sizeof(SNR_list[0]); // calculate the number of elements in the SNR_list
const int numProblems = 10; // 100;

const int N_t = 16;   // Number Of Transmit antennas

// define the testbench path
const std::string GitHub_Repo_Path = "/home/hy592/MIMO";  // Path where GitHub Repo is cloned
const std::string Git_Repo_Name = "/6775-Final_Project";                           // the GitHub Repo Name
const std::string PATH_tb_data = "tb_data";                         // Path to testbench data

// the following are used to set the path to the testbench data
// std::string fullPathTo_Tb_data = GitHub_Repo_Path + Git_Repo_Name + PATH_tb_data; 
std::string fullPathTo_Tb_data = PATH_tb_data; 
std::string fullPathTo_fixed_pt_data = fullPathTo_Tb_data + "/HLS_tb_data/fixed_pt_data";
// std::string fullPathTo_x_init = fullPathTo_Tb_data + "/init_spins.txt";
// std::string fullPathTo_x_init = fullPathTo_Tb_data + "/HLS_tb_data/init_spins_big.txt"; // old file
std::string fullPathTo_x_init = fullPathTo_Tb_data + "/HLS_tb_data/x_file.txt"; // old file
std::string fullPathTo_error_init = fullPathTo_Tb_data + "/HLS_tb_data/error_var_file.txt"; // old file

// std::string fullPathTo_IDEALS = fullPathTo_Tb_data + "/sample_MIMO/ideal_sols_16_nr16_16QAM";
std::string fullPathTo_JMATRICES = fullPathTo_Tb_data + "/instances/Nt16_Nr16_16QAM";
// for example, the path to the J matrix for problem 0 is:
// C:\Users\user\Documents\GitHub\MIMO_Ising\fpga\tb_data\sample_MIMO\J_matrix_16_nr16_16QAM_snr_30\DI_MIMO_J_Quad_0.txt

std::string fullPathTo_hls_solved_solutions = fullPathTo_Tb_data + "/HLS_solved_solutions/Nt16_Nr16_16QAM";

const int matrix_size = 65;

void top_func_test( 
  data_type_J J[N][N], 
  data_type_x init_spins[num_anneals][N], 
  data_type_e error_var_init[num_anneals][N],
  string SNR,
  int k // problem number 
)
{
  spin_sign bestSpins[N];
  data_type_e bestEnergy;

  AHC ahc_instance;           // AHC instance
  ahc_instance.updateJ( J );  // load problem set
  // ahc_instance.ahc_solver();
  for ( int i = 0; i < num_anneals; i++ ) {
    // ahc_instance.updateSpins( init_spins[i], 0.02,1.3 );
    ahc_instance.ahc_solver( init_spins[i], error_var_init[i] );
  }

  bestEnergy = ahc_instance.bestEnergySpins( bestSpins );
  // spin_sign spin_guesses[N];

  std::cout << "BEST ENERGY = " << bestEnergy << std::endl;

  // std::cout << "Spin = ";
  // write to fullPathTo_hls_solved_solutions/{SNR}/DI_MIMO_sol_{k}.txt
  std::ofstream outSpins_file;
  outSpins_file.open( fullPathTo_hls_solved_solutions + "/" + SNR + "/DI_MIMO_sol_" + std::to_string(k) + ".txt" );
  
  if(!outSpins_file.is_open()) {
    std::cout << "Error opening outSpins_file files!" << std::endl;
    return;
  }

  // write the best spins to a file
  for ( int i = 0; i < N; ++i ) {
    outSpins_file << bestSpins[i] << std::endl;
    // printf("%i, ", bestSpins[i]);
  }
  printf("\n");
  outSpins_file.close();
  // std::cout << std::endl;

  /* ifstream ideal_solutions(IDEALS + "/DI_MIMO_mmse_sol_" + std::to_string(problemNum) + ".txt");

  spin_sign ideal_spins[N]; */
  // ofstream best_spins_out("best_spins_" + std::to_string(problemNum) + ".txt");
  // best_spins_out << "[ ";
  // for (int i = 0; i < N; ++i)
  // {
  //   best_spins_out << bestSpins[i] << ", ";
  // }
  // best_spins_out << "]";

  // for very high SNR, calculate bestSpins[0:n]+bestSpins[n:2*n], where n=2*N_t. All elements should be 0
  // spin_sign sol_out[2*N_t];
  // std::cout << "bestSpins[0:n]+bestSpins[n:2*n] = ";
  // for (int i = 0; i < 2*N_t; i++)
  // {
  //   sol_out[i] = bestSpins[i] + bestSpins[i+2*N_t];
  //   std::cout << sol_out[i] << ", ";
  // }
  // std::cout << std::endl;

  // outSpins.close();  // Don't forget to close the file
}

int main()
{
  data_type_x x_init_arrays[num_anneals][matrix_size];
  data_type_e error_var_init_arrays[num_anneals][N];

  data_type_J J_matrix[matrix_size][matrix_size];

  // Seed the random number generator
  srand( 1234 );

  // random_array not used
  // data_type_x random_array[matrix_size] = { 3.58964059e-05, 1.63794645e-04, 1.48891121e-05, 4.44594756e-04, 8.65550405e-05, 4.03401915e-04, -3.62525296e-04, -3.60723653e-04, 3.07391289e-04, -1.02323163e-04, -3.34645803e-04, 4.27508580e-04, -1.52234140e-04, 2.50812103e-04, 2.25997985e-04, 3.83306091e-04, 1.23672207e-04, 2.50942434e-04, -1.51101658e-04, -2.30072108e-04, 3.95886218e-04, -7.19088101e-05, 4.64840047e-04, 1.63441498e-04, 1.21695720e-04, -3.85254027e-04, 4.49489259e-04, -5.00878665e-05, 7.83896144e-05, -9.18631972e-05, -2.62973020e-04, 4.03379521e-04, 7.36794867e-05, -4.97129673e-04, 1.17144914e-04, -1.73355098e-04, 2.70581023e-05, 3.85942099e-04, -1.42730240e-04, 4.08535151e-04, 1.23360116e-04, -4.84178757e-04, 4.29437234e-04, 1.90896918e-04, 4.97322850e-04, -3.27659492e-04, -3.62864250e-04, 4.32595463e-04, 1.96818161e-04, -4.33999827e-04, 2.55463053e-04, 2.53876188e-04, 4.23024536e-04, 2.11524759e-04, -3.75729038e-04, -4.80119866e-04, -4.73789013e-04, -4.71693512e-04, -2.53788932e-04, 3.60027949e-04, 3.88310643e-05, 5.28219787e-05, 3.42030892e-04, -3.75826685e-04, -2.20816321e-04 };

  // cout << "INIT SPINS -----------------" << std::endl;

  // load the x_init arrays
  ifstream x_init_file( fullPathTo_x_init );
  std::cout << fullPathTo_x_init << std::endl;
  if ( !x_init_file.is_open() ) {
    std::cout << "Error opening x_init file!" << std::endl;
    return 1;
  }

  // Define a 2D array to store all 20 x_init arrays
  for ( int array_idx = 0; array_idx < num_anneals; array_idx++ ) {
    for ( int i = 0; i < matrix_size; i++ ) {
      if ( !( x_init_file >> x_init_arrays[array_idx][i] ) ) {
        std::cout << "Error reading x_init file!" << std::endl;
        return 1;
      }
    }
  }
  x_init_file.close();

  // Load the J matrix values
  printf("There are %d SNR values to solve: [", numSNR);

  // print each SNR value
  for (int SNR_idx = 0; SNR_idx < numSNR; SNR_idx++){
    printf("%d, ", SNR_list[SNR_idx]);
  }
  printf("]\n");
  

  for (int SNR_idx = 0; SNR_idx < numSNR; SNR_idx++){
    SNR = to_string(SNR_list[SNR_idx]);
    //  Do something with the J_matrix values here...
    int startNum = 0;
    for ( int k = startNum; k < startNum + numProblems; k++ ) {
      printf("Solving SNR=%d, problem %d\n", SNR_list[SNR_idx], k);
      std::string J_Matrix_file = SNR+"/DI_MIMO_J_Quad_"+to_string(k)+".txt";
      ifstream input_file1( fullPathTo_JMATRICES + "/" + J_Matrix_file );
      if ( !input_file1.is_open() ) {
        std::cout << "Error opening J_Matrix file!" << std::endl;
        return 1;
      }

      // std::cout << fullPathTo_JMATRICES << "/" << J_Matrix_file << std::endl;
      for ( int i = 0; i < matrix_size; i++ ) {
        for ( int j = 0; j < matrix_size; j++ ) {
          float J_in;
          input_file1 >> J_in;
          J_matrix[i][j] = J_in;
        }
      }
      input_file1.close();
      top_func_test( J_matrix, x_init_arrays, error_var_init_arrays, SNR, k );
    }
  }
  return 0;
}
