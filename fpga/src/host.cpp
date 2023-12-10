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

std::string data_path = "/root/GitHub_Linux/MIMO_Ising/data";

const int matrix_size = N; // Size of the matrix, N

void comp_ref(data_type_J J[N][N], data_type_x init_spins[100][N],
              data_type_e &energy_ref[numProblems], spin_sign &spin_ref[numProblems][N], int problemNum)
{
  AHC solver(J);
  // solver.ahc_solver();
  for (int i = 0; i < 100; i++) {
    // run the AHC solver
    solver.ahc_solver(init_spins[i]);
  }
  energy_ref[problemNum] = solver.bestEnergy;
  for (int i = 0; i < N; i++) {
    spin_ref[problemNum][i] = solver.bestSpins[i];
  }
}

int main(int argc, char **argv) {
  // Define the J matrix and the initial spins
  data_type_J matrix[numProblems][matrix_size][matrix_size];
  data_type_x x_init_arrays[100][matrix_size];
  
  // read the initial spins from the file
  ifstream x_init_file("init_spins_big.txt");
  if (!x_init_file.is_open())
  {
    // error handling
    cout << "Error opening x_init file!" << endl;
    return 1;
  }

  // Define a 2D array to store all 20 x_init arrays
  for (int array_idx = 0; array_idx < 100; array_idx++){
    for (int i = 0; i < matrix_size; i++){
      // read the x_init values from the file
      if (!(x_init_file >> x_init_arrays[array_idx][i])){
        // error handling
        cout << "Error reading x_init file!" << endl;
        return 1;
      }
    }
  }
  x_init_file.close();

  std::string file_directory_instances = data_path+"di_mimo\\instances\\Nt16_Nr16_16QAM\\30\\";
  for (int k = 0; k < numProblems; k++) {
    std::string file_name = "DI_MIMO_J_Conv_6_19_12_2_" + std::to_string(k) + ".txt";
    ifstream input_file1(file_directory_instances + file_name);

    if (!input_file1.is_open()){
      cout << "Error opening file!" << endl;
      return 1;
    }

    for (int i = 0; i < matrix_size; i++){
      for (int j = 0; j < matrix_size; j++){
        int hex_value;
        input_file1 >> hex >> hex_value;

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
  timer.start();
  //--------------------------------------------------------------------
  // Send data to accelerator
  //--------------------------------------------------------------------
  for (int k = 0; k < numProblems; ++k) {
    for (int i = 0; i < matrix_size; i++) {
      for (int j = 0; j < matrix_size; j++) {
        test_data_J = matrix[k][i][j];
        nbytes = write(fdw, (void *)&test_data_J, sizeof(test_data_J));
        assert(nbytes == sizeof(test_data_J));
      }
    }
    for (int i = 0; i < 100; ++i) {
      for (int j = 0; j < matrix_size; ++j) {
        test_data_X = x_init_arrays[i][j];
        nbytes = write(fdw, (void *)&test_data_X, sizeof(test_data_X));
        assert(nbytes == sizeof(test_data_X));  
      }
    }
  }
  //--------------------------------------------------------------------
  // Receive data from accelerator
  //--------------------------------------------------------------------
  for (int k = 0; k < numProblems; ++k) {
    for (int i = 0; i < matrix_size + 1; ++i) {
      if (i == 0) {
        nbytes = read(fdr, (void *)&(energy_fpga[k]), sizeof(energy_fpga[k]));
        assert(nbytes == sizeof(energy_fpga[k]));
      }
      nbytes = read(fdr, (void *)&(spin_fpga[k][i]), sizeof(spin_fpga[k][i]));
      assert(nbytes == sizeof(spin_fpga[k][i]));
    }
  }
  timer.stop();
  for (int k = 0; k < numProblems; ++k) {
    comp_ref(matrix[k], x_init_arrays, energy_ref[numProblems], spin_ref[numProblems][N], k);
  }

  float energy_error = 0;
  int spin_error = 0;
  for (int k = 0; k < numProblems; ++k) {
    energy_error += abs(energy_fpga - energy_ref);
    for (int i = 0; i < 100; ++i) {
      spin_error += abs(spin_fpga - spin_ref);
    }
  }
  
  std::cout << "Energy Testing size: " << numProblems << std::endl;
  std::cout << "Accuracy: " << energy_error / numProblems << std::endl;

  std::cout << "Spin Testing size: " << numProblems * 100 << std::endl;
  std::cout << "Accuracy: " << spin_error / numProblems * 100 << std::endl;
  return 0; // Return 0 if everything executed properly
}