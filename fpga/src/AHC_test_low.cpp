#include <cstdlib> // Required for rand() function
#include <iostream>
#include <fstream>
#include <random>

#include "AHC_low.h"

using namespace std;

// set the data_path where di_mimo instances are stored
// In my computer, it is C:\Users\user\Documents\GitHub\DesignProject\IsingMachine\data
//std::string data_path = "C:\\Users\\user\\Documents\\GitHub\\DesignProject\\IsingMachine\\data\\";
std::string data_path = "/root/GitHub_Linux/MIMO_Ising/data";

const int matrix_size = N; // Size of the matrix, N

/*
  This function is used to test the AHC solver
  It takes in the J matrix and the initial spins
  It then runs the AHC solver for 100 iterations
  and outputs the best spins to a file
*/
void top_func_test(data_type_J J[N][N], data_type_x init_spins[100][N], int problemNum)
{
  AHC solver(J);
  // solver.ahc_solver();
  for (int i = 0; i < 100; i++)
  {
    // initialize the temperature and the annealing rate
    // if (i < 100) {
    //   solver.updateSpins(init_spins[i], 0.02, 1.3);
    // }
    /*
    else if(i < 60){
      solver.updateSpins(init_spins[i], 0.02,1.2);
    }
    else if(i < 80){
      solver.updateSpins(init_spins[i], 0.02,1);
    }
    else if(i < 100){
      solver.updateSpins(init_spins[i], 0.024, 1);
    }
    */

    // run the AHC solver
    solver.ahc_solver(init_spins[i]);
  }
  // output filestream object
  std::ofstream outSpins; 
  // path to the directory where the solution file will be stored
  std::string file_directory_hls_sols = data_path+"di_mimo\\hls_sols\\";
  // name of the solution file
  std::string file_name = "DI_MIMO_J_Conv_snr_30_12_2_" + std::to_string(problemNum) + ".txt";

  // Print the best energy and spins
  std::cout << "BEST ENERGY" << solver.bestEnergy << endl;

  // Write the best energy to the file
  outSpins.open(file_directory_hls_sols + file_name);
  for (int i = 0; i < N; i++){
    outSpins << solver.bestSpins[i] << " "; // Separate values by space
  }
  outSpins.close(); // Don't forget to close the file
}


int main() {
  // Define the J matrix and the initial spins
  data_type_J matrix[matrix_size][matrix_size];
  data_type_x x_init_arrays[100][matrix_size];
  
  // test sign bit
  data_type_J zed = 0.0;
  if (zed > 0){
    cout << "ZERO IS bigger than zpointz" << endl;
  }
  else if (zed == 0)
    cout << "ZEROdotZero IS ZED" << endl;
  // else if (sign_bit(zed) == 1)
  //   cout << " SIGN BIT ONE" << endl;
  else if (zed < 0.0)
    cout << "zed less zpointz" << endl;
  else
    cout << "SIGN BIT ISSUE" << endl;

  // Seed the random number generator
  srand(1234);
  /*
    // Fill the random array with values of either 0 or 1
    for (int i = 0; i < matrix_size; i++) {
      random_array[i] = static_cast<data_type_x>(0.002 * (rand() / (double)RAND_MAX) - 1.0);
      cout << "INIT SPINS" <<endl;
      std::cout << random_array[i] << ", ";
      cout << endl;

    }
    cout << endl;
  */
  /*
    const double min_val = -0.0005;
    const double max_val = 0.0005;

    std::uniform_real_distribution<double> dist(min_val, max_val);
    std::mt19937_64 rng(1234);
    std::cout << "INIT SPINS" << std::endl;

    for (int i = 0; i < N; i++) {
        random_array[i] = 0.01 * dist(rng);
        std::cout << random_array[i] << ", ";
    }
    std::cout << std::endl;
  */

  /*
    data_type_x random_array[matrix_size]= {0.0 ,-0.000488281,0.0, 0.0, -0.000488281, 0.0, 0.0, 0.0, -0.000488281, 0.0, 0.0, -0.000488281, -0.000488281,
            0.0, -0.000488281, 0.0, -0.000488281, -0.000488281, -0.000488281, 0.0, -0.000488281, 0.0, -0.000488281,
            0.0, -0.000488281, 0.0, 0.0, -0.000488281, -0.000488281, 0.0, 0.0, -0.000488281, 0.0, 0.0, 0.0, -0.000488281,
            -0.000488281, 0.0, 0.0, 0.0, -0.000488281, 0.0, 0.0, 0.0, -0.000488281, -0.000488281, -0.000488281,
            -0.000488281, 0.0, 0.0, 0.0, 0.0, -0.000488281, 0.0, -0.000488281, -0.000488281, -0.000488281, -0.000488281,
            -0.000488281, 0.0, 0.0, 0.0, 0.0, -0.000488281, 0.0};
  */

  data_type_x random_array[matrix_size]={
    3.58964059e-05,1.63794645e-04,1.48891121e-05,4.44594756e-04
    ,8.65550405e-05,4.03401915e-04,-3.62525296e-04,-3.60723653e-04
    ,3.07391289e-04,-1.02323163e-04,-3.34645803e-04,4.27508580e-04
    ,-1.52234140e-04,2.50812103e-04,2.25997985e-04,3.83306091e-04
    ,1.23672207e-04,2.50942434e-04,-1.51101658e-04,-2.30072108e-04
    ,3.95886218e-04,-7.19088101e-05,4.64840047e-04,1.63441498e-04
    ,1.21695720e-04,-3.85254027e-04,4.49489259e-04,-5.00878665e-05
    ,7.83896144e-05,-9.18631972e-05,-2.62973020e-04,4.03379521e-04
    ,7.36794867e-05,-4.97129673e-04,1.17144914e-04,-1.73355098e-04
    ,2.70581023e-05,3.85942099e-04,-1.42730240e-04,4.08535151e-04
    ,1.23360116e-04,-4.84178757e-04,4.29437234e-04,1.90896918e-04
    ,4.97322850e-04,-3.27659492e-04,-3.62864250e-04,4.32595463e-04
    ,1.96818161e-04,-4.33999827e-04,2.55463053e-04,2.53876188e-04
    ,4.23024536e-04,2.11524759e-04,-3.75729038e-04,-4.80119866e-04
    ,-4.73789013e-04,-4.71693512e-04,-2.53788932e-04,3.60027949e-04
    ,3.88310643e-05,5.28219787e-05,3.42030892e-04,-3.75826685e-04
    ,-2.20816321e-04
  };

  //cout << "INIT SPINS -----------------" << endl;

  // Read the initial spins from the file
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

  // Define the number of problems to solve
  int numProblems = 10;
  int startNum = 0;

  // cout << "HELLO" << endl;
  //  Do something with the matrix values here...

  // ifstream input_file("DI_MIMO_J_Conv_Fixed_Mult0.txt");

  std::string file_directory_instances = data_path+"di_mimo\\instances\\Nt16_Nr16_16QAM\\30\\";
  for (int k = startNum; k < startNum + numProblems; k++){
    std::string file_name = "DI_MIMO_J_Conv_6_19_12_2_" + std::to_string(k) + ".txt";
    ifstream input_file1(file_directory_instances + file_name);

    if (!input_file1.is_open()){
      cout << "Error opening file!" << endl;
      return 1;
    }

    #ifdef DEBUG
        cout << "INIT SPINS" << endl;
        for (int z = 0; z < matrix_size; z++){
          cout << random_array[z] << " ";
        }
        cout << endl;
    #endif
    /*
      for(int j=0; j < 20; ++j){
        for(int i=0; i < N; ++i){
          cout << x_init_arrays[j][i] << " ";
        }
        cout << std::endl;
      }

      */

    for (int i = 0; i < matrix_size; i++){
      for (int j = 0; j < matrix_size; j++){
        int hex_value;
        input_file1 >> hex >> hex_value;
        // cout << "HEX val" << hex_value << endl;
        // matrix[i][j] = static_cast<float>(hex_value) / 1048576.0f;
        // matrix[i][j] = static_cast<float>(hex_value) / 262144.0f; //18_2
        // matrix[i][j] = static_cast<float>(hex_value) / 8388608.0f; //23_2

        matrix[i][j] = static_cast<float>(hex_value) / 4096.0f; // 12_2

        // matrix[i][j] = static_cast<float>(hex_value) / 4096.0f;

        // #ifdef DEBUG
        //   cout << matrix[i][j] << "  ";
        // #endif
      }
      // #ifdef DEBUG
      //	cout << std::endl;
      // #endif DEBUG
    }
    top_func_test(matrix, x_init_arrays, k);
    input_file1.close();
  }

  /*
  std::string file_directory_instances = "C:\\Users\\ari\\Documents\\di_mimo\\instances\\Nt16_Nr16_16QAM\\30\\";
  std::string file_name = "DI_MIMO_J_Conv_18_2_0.txt";
  ifstream input_file0(file_directory_instances + file_name);

   if (!input_file0.is_open()) {
     cout << "Error opening file!" << endl;
     return 1;
   }
  for (int i = 0; i < matrix_size; i++) {
    for (int j = 0; j < matrix_size; j++) {
       int hex_value;
      input_file0 >> hex >> hex_value;
      //cout << "HEX val" << hex_value << endl;
      //matrix[i][j] = static_cast<float>(hex_value) / 1048576.0f;
      matrix[i][j] = static_cast<float>(hex_value) / 262144.0f; //18_2

     // matrix[i][j] = static_cast<float>(hex_value) / 4096.0f;

      cout << matrix[i][j] << "  ";

    }
    cout << std::endl;
  }
  top_func_test(matrix, x_init_arrays, 0);

  input_file0.close();
*/
  return 0; // Return 0 if everything executed properly
}
