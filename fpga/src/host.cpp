#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <math.h>
#include <assert.h>

#include <iostream>
#include <fstream>

#include "AHC_low.h"
#include "timer.h"

using namespace std;

const int matrix_size = N; // Size of the matrix, N

// define a list of SNR values to solve
int SNR_list[] = {0, 5, 10, 15, 20, 25, 30};
const int numSNR = sizeof(SNR_list)/sizeof(SNR_list[0]); // calculate the number of elements in the SNR_list
const int numProblems = 100;

const int N_t = 16;   // Number Of Transmit antennas

// define the testbench path
const std::string GitHub_Repo_Path = "C:\\Users\\user\\Documents\\GitHub";  // Path where GitHub Repo is cloned
const std::string Git_Repo_Name = "\\MIMO_ISING";                           // the GitHub Repo Name
const std::string PATH_tb_data = "\\tb_data";                         // Path to testbench data

// the following are used to set the path to the testbench data
std::string fullPathTo_Tb_data = GitHub_Repo_Path + Git_Repo_Name + PATH_tb_data; // MIMO_ISING/tb_data
std::string fullPathTo_fixed_pt_data = fullPathTo_Tb_data + "\\HLS_tb_data\\fixed_pt_data";
// std::string fullPathTo_x_init = fullPathTo_Tb_data + "\\init_spins.txt";
// std::string fullPathTo_x_init = fullPathTo_Tb_data + "\\HLS_tb_data\\init_spins_big.txt"; // old file
std::string fullPathTo_x_init = fullPathTo_Tb_data + "\\HLS_tb_data\\x_file.txt"; // old file
std::string fullPathTo_error_init = fullPathTo_Tb_data + "\\HLS_tb_data\\error_var_file.txt"; // old file

// std::string fullPathTo_IDEALS = fullPathTo_Tb_data + "\\sample_MIMO\\ideal_sols_16_nr16_16QAM";
std::string fullPathTo_JMATRICES = fullPathTo_Tb_data + "\\instances\\Nt16_Nr16_16QAM";
// for example, the path to the J matrix for problem 0 is:
// C:\Users\user\Documents\GitHub\MIMO_Ising\fpga\tb_data\sample_MIMO\J_matrix_16_nr16_16QAM_snr_30\DI_MIMO_J_Quad_0.txt

std::string fullPathTo_hls_solved_solutions = fullPathTo_Tb_data + "\\HLS_solved_solutions\\Nt16_Nr16_16QAM";


void top_func_test_6775(
  data_type_J J[N][N], 
  data_type_x init_spins[num_anneals][N], 
  data_type_e error_var_init[num_anneals][N],
  string SNR,
  int k // problem number 
){
  // to receive data from FPGA
  spin_sign bestSpins[N];
  data_type_e bestEnergy;

  bit_Width_t energy_fpga;
  bit2_t spin_fpga[matrix_size];

  // send data to FPGA
  data_type_J test_data_J;
  data_type_x test_data_X;

  //--------------------------------------------------------------------
  Timer timer("Ising on FPGA");
  bit32_t nbytes;
  bit32_t data_write;
  data_write = 0;
  //--------------------------------------------------------------------

  //--------------------------------------------------------------------
  // Send data to accelerator
  //--------------------------------------------------------------------
  timer.start();
  std::cout << "Start FPGA Send Data" << std::endl;

  // std::cout << "k=" << k << std::endl;
  // std::cout << "send J"<< std::endl;
  for (int i = 0; i < matrix_size; i++) {
    for (int j = 0; j < matrix_size; j++) {
      test_data_J = J_Matrix[i][j];
      data_write(MAX_WIDTH-1,0) = test_data_J(MAX_WIDTH-1,0);
      nbytes = write(fdw, (void *)&data_write, sizeof(data_write));
      assert(nbytes == sizeof(data_write));
    }
  }
  // std::cout << "send J finish" << k << std::endl;
  for (int i = 0; i < num_anneals; i++) {
    // std::cout << "  i=" << i << std::endl;
    // std::cout << "  send x" << std::endl;
    for (int j = 0; j < matrix_size; j++) {
      test_data_X = x_init_arrays[i][j];
      data_write(MAX_WIDTH-1,0) = test_data_X(MAX_WIDTH-1,0);
      nbytes = write(fdw, (void *)&data_write, sizeof(data_write));
      assert(nbytes == sizeof(data_write));  
      // test_data_X = data_write(15, 0);
      // std::cout << "test=" << data_write << std::endl;
    }
    // std::cout << "  send x finish" << std::endl;
  }
  std::cout << "Finish FPGA send Data" << std::endl;

  //--------------------------------------------------------------------
  // Receive data from accelerator
  //--------------------------------------------------------------------
  std::cout << "Receive Output" << std::endl;
  std::cout << "k=" << k << std::endl;
  bit32_t energy_received;
  data_type_e energy_result;

  nbytes = read(fdr, (void *)&(energy_received), sizeof(energy_received));
  assert(nbytes == sizeof(energy_received));
  energy_fpga[k] = energy_received(MAX_WIDTH-1,0);
  energy_result(MAX_WIDTH-1,0) = energy_received(MAX_WIDTH-1,0);
  std::cout << "BEST ENERGY = " << energy_result << std::endl;

  std::cout << "Spin = " << std::endl;
  bit32_t spins_received;
  spin_sign spin_result;
  for (int i = 0; i < 8; ++i) {
    nbytes = read(fdr, (void *)&(spins_received), sizeof(spins_received));
    assert(nbytes == sizeof(spins_received));
    for (int j=0; j<8; j++){
      bit2_t temp_value;
      temp_value = (spins_received >> (2 * j)) & 0b11;
      // spin_result = (spins_received >> (2 * j)) & 0b11;
      spin_result = reinterpret_cast<ap_int<2>&>(temp_value);
      spin_fpga[k][8*i+j] = spin_result;
      std::cout << spin_result << std::endl;
    }
  }
  nbytes = read(fdr, (void *)&spins_received, sizeof(spins_received));
  assert(nbytes == sizeof(spins_received));
  bit2_t temp_value;
  temp_value = spins_received(1,0);
  spin_result = reinterpret_cast<ap_int<2>&>(temp_value);
  spin_fpga[k][64] = spin_result;
  std::cout << spin_result << std::endl;
  
  std::cout << "End Receive Output" << std::endl;
  timer.stop();

  //--------------------------------------------------------------------
  // Write Result into File
  //--------------------------------------------------------------------

  // std::cout << "Spin = ";
  // write to fullPathTo_hls_solved_solutions\\{SNR}\\DI_MIMO_sol_{k}.txt
  std::ofstream outSpins_file;
  outSpins_file.open( fullPathTo_hls_solved_solutions + "\\" + SNR + "\\DI_MIMO_sol_" + std::to_string(k) + ".txt" );
  
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

  std::cout << "Finish \n" << std::endl;
}


int main(int argc, char **argv) {
  int fdr = open("/dev/xillybus_read_32", O_RDONLY);
  int fdw = open("/dev/xillybus_write_32", O_WRONLY);  

  // Check that the channels are correctly opened
  if ((fdr < 0) || (fdw < 0)) {
    fprintf(stderr, "Failed to open Xillybus device channels\n");
    return -1;
  }

  data_type_x x_init_arrays[num_anneals][matrix_size];
  data_type_J J_matrix[matrix_size][matrix_size];
  // data_type_e error_var_init_arrays[num_anneals][N];

  // load the x_init arrays
  ifstream x_init_file( fullPathTo_x_init );
  if ( !x_init_file.is_open() ) {
    cout << "Error opening x_init file!" << endl;
    return 1;
  }
  else{
    // Define a 2D array to store all num_anneals x_init arrays
    for ( int array_idx = 0; array_idx < num_anneals; array_idx++ ) {
      for ( int i = 0; i < matrix_size; i++ ) {
        if ( !( x_init_file >> x_init_arrays[array_idx][i] ) ) {
          cout << "Error reading x_init file!" << endl;
          return 1;
        }
      }
    }
    x_init_file.close();
  }

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
      std::string J_Matrix_file = SNR+"\\DI_MIMO_J_Quad_"+to_string(k)+".txt";
      ifstream input_file1( fullPathTo_JMATRICES + "\\" + J_Matrix_file );
      if ( !input_file1.is_open() ) {
        cout << "Error opening J_Matrix file!" << endl;
        return 1;
      }

      // cout << fullPathTo_JMATRICES << "\\" << J_Matrix_file << endl;
      for ( int i = 0; i < matrix_size; i++ ) {
        for ( int j = 0; j < matrix_size; j++ ) {
          float J_in;
          input_file1 >> J_in;
          J_matrix[i][j] = J_in;
        }
      }
      input_file1.close();
      top_func_test_6775( J_matrix, x_init_arrays, error_var_init_arrays, SNR, k );
    }
  }
  return 0;
}
