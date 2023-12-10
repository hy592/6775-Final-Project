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

const int matrix_size = N; // Size of the matrix, N

// void comp_ref(data_type_J J[N][N], data_type_x init_spins[100][N],
//               data_type_e energy_ref[numProblems], spin_sign spin_ref[numProblems][N], int problemNum)
// {
//   AHC ahc_instance;  // AHC instance
//   ahc_instance.updateJ(J);    // load problem set
//   // solver.ahc_solver();
//   for (int i = 0; i < 100; i++) {
//     // run the AHC solver
//     ahc_instance.ahc_solver(init_spins[i]);
//   }
//   energy_ref[problemNum] = ahc_instance.bestEnergySpins(spin_ref[problemNum]);
// }

int main(int argc, char **argv) {
  // Define the J matrix and the initial spins
  data_type_J matrix[numProblems][matrix_size][matrix_size];
  data_type_x x_init_arrays[100][matrix_size];

  int fdr = open("/dev/xillybus_read_32", O_RDONLY);
  int fdw = open("/dev/xillybus_write_32", O_WRONLY);  

  // read the initial spins from the file
  std::ifstream x_init_file("init_spins_big.txt");
  if (!x_init_file.is_open())
  {
    // error handling
    std::cout << "Error opening x_init file!" << std::endl;
    return 1;
  }

  // Define a 2D array to store all 20 x_init arrays
  for (int array_idx = 0; array_idx < 100; array_idx++){
    for (int i = 0; i < matrix_size; i++){
      // read the x_init values from the file
      if (!(x_init_file >> x_init_arrays[array_idx][i])){
        // error handling
        std::cout << "Error reading x_init file!" << std::endl;
        return 1;
      }
    }
  }
  x_init_file.close();

  for (int k = 0; k < numProblems; k++) {
    // std::string file_path = (string)"fixed_pt_data/" + (string)"DI_MIMO_J_Conv_11_29_16_2_" + std::to_string(k) + (string)".txt";
    // ifstream input_file1((string)"fixed_pt_data/" + (string)"DI_MIMO_J_Conv_11_29_16_2_" + std::to_string(k) + (string)".txt");
    // Use stringstream to concatenate integer and string
    std::ifstream input_file1("fixed_pt_data/DI_MIMO_J_Conv_11_29_16_2_0.txt");

    if (!input_file1.is_open()){
      std::cout << "Error opening file!" << std::endl;
      return 1;
    }

    for (int i = 0; i < matrix_size; i++){
      for (int j = 0; j < matrix_size; j++){
        int hex_value;
        input_file1 >> std::hex >> hex_value;

        matrix[k][i][j] = static_cast<float>(hex_value) / 4096.0f; // 12_2
      }
    }
    input_file1.close();
  }
  
  data_type_J test_data_J;
  data_type_x test_data_X;
  data_type_e energy_fpga[numProblems];
  spin_sign spin_fpga[numProblems][matrix_size];

  data_type_e energy_ref[numProblems];
  spin_sign spin_ref[numProblems][matrix_size];

  Timer timer("Ising on FPGA");
  bit32_t nbytes;
  timer.start();
  //--------------------------------------------------------------------
  // Send data to accelerator
  //--------------------------------------------------------------------
  std::cout << "Start FPGA Send Data" << std::endl;
  for (int k = 0; k < numProblems; ++k) {
    std::cout << "k=" << k << std::endl;
    std::cout << "send J"<< std::endl;
    for (int i = 0; i < matrix_size; i++) {
      for (int j = 0; j < matrix_size; j++) {
        test_data_J = matrix[k][i][j];
        nbytes = write(fdw, (void *)&test_data_J, sizeof(test_data_J));
        assert(nbytes == sizeof(test_data_J));
      }
    }
    std::cout << "send J finish" << k << std::endl;
    for (int i = 0; i < 20; ++i) {
      std::cout << "  i=" << i << std::endl;
      std::cout << "  send x" << std::endl;
      for (int j = 0; j < matrix_size; ++j) {
        test_data_X = x_init_arrays[i][j];
        nbytes = write(fdw, (void *)&test_data_X, sizeof(test_data_X));
        assert(nbytes == sizeof(test_data_X));  
      }
      std::cout << "  send x finish" << std::endl;
    }
  }
  std::cout << "Finish FPGA send Data" << std::endl;
  //--------------------------------------------------------------------
  // Receive data from accelerator
  //--------------------------------------------------------------------
  std::cout << "Receive Output" << std::endl;
  for (int k = 0; k < numProblems; ++k) {
    std::cout << "k=" << k << std::endl;

    nbytes = read(fdr, (void *)&(energy_fpga[k]), sizeof(energy_fpga[k]));
    assert(nbytes == sizeof(energy_fpga[k]));
    std::cout << "  Energy = " << energy_fpga[k] << std::endl;

    std::cout << "  Spin = " << std::endl;
    for (int i = 0; i < matrix_size; ++i) {
      nbytes = read(fdr, (void *)&(spin_fpga[k][i]), sizeof(spin_fpga[k][i]));
      assert(nbytes == sizeof(spin_fpga[k][i]));
      std::cout << spin_fpga[k][i] << std::endl;
    }
  }
  std::cout << "End Receive Output" << std::endl;
  timer.stop();
  // for (int k = 0; k < numProblems; ++k) {
  //   comp_ref(matrix[k], x_init_arrays, energy_ref, spin_ref, k);
  // }

  // float energy_error = 0;
  // int spin_error = 0;
  // for (int k = 0; k < numProblems; ++k) {
  //   energy_error += abs(energy_fpga - energy_ref);
  //   for (int i = 0; i < 100; ++i) {
  //     spin_error += abs(spin_fpga - spin_ref);
  //   }
  // }
  
  // std::cout << "Energy Testing size: " << numProblems << std::endl;
  // std::cout << "Accuracy: " << energy_error / numProblems << std::endl;

  // std::cout << "Spin Testing size: " << numProblems * 100 << std::endl;
  // std::cout << "Accuracy: " << spin_error / numProblems * 100 << std::endl;

  std::cout << "Finish \n" << std::endl;
  return 0; // Return 0 if everything executed properly
}