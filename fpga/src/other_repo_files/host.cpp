//
//  ksat_host.cpp
//  Ising-k-SAT
//

// Xilinx libraries
#include <ap_fixed.h>
#include <xcl2.hpp>

// host code libraries
#include <math.h>
#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <vector>

// printf for debugging
#include <stdint.h>
#include <inttypes.h>

// command line parser library
#include "cxxopts.hpp"

// CNF parser library
#include "CNFparser.hpp"

// constant definitions
#include "solver_compile_consts.hpp"

void printfixvec( std::vector<fixed_point_32t, aligned_allocator<fixed_point_32t>> const &input, std::string name )
{
  std::cout << "\nPrinting vector " << name << ": ";
  for ( unsigned int i = 0; i < input.size(); i++ ) {
    std::cout << (float)(input.at( i )) << ' ';
  }

  std::cout << "\n"
            << std::flush;
}

int sat_a_num(fixed_point_32t num){
  if ((float)(num) > 0)
    return 1;
  return -1;
}

void printsol( std::vector<fixed_point_32t, aligned_allocator<fixed_point_32t>> const &input, std::string name )
{
  std::cout << "\n" << name << ": ";
  for ( unsigned int i = 0; i < input.size(); i++ ) {
    std::cout << (sat_a_num(input.at( i ))) << ' ';
  }

  std::cout << "\n"
            << std::flush;
}

template <typename T>
inline int sgn( T val )
{
  return ( T( 0 ) < val ) - ( val < T( 0 ) );
}

int main( int argc, char **argv )
{
  cl_int err;

  //--------------------------------------
  // initializations
  //--------------------------------------

  // deal with command-line arguments
  cxxopts::Options options("AnalogkSat", "Analog kSat solver with bounded algorithm (McMahon Lab @ Cornell Univ)");

  options.add_options()
    ("p,problem_path", "problem matrix directory", cxxopts::value<std::string>())
    ("x,xclbin_path", "path of xclbin", cxxopts::value<std::string>())
    ("s,num_steps", "number of steps between re-initialization of state vectors", cxxopts::value<int>()->default_value("10000"))
    ("q,quantum", "number of steps between solution checking", cxxopts::value<int>()->default_value("500"))
    ("w,is_wcnf", "1 for WCNF format and 0 for CNF format", cxxopts::value<int>()->default_value("0"))
    ("e,early_stop", "whether to stop early if all clauses are satisfied", cxxopts::value<int>()->default_value("1"))
    ("t,timeout", "how long to run before stopping (in us)", cxxopts::value<int>()->default_value("300000000")) // 300 seconds
    ("d,debug", "enable debug mode (takes longer in runtime)", cxxopts::value<int>()->default_value("0"))
    ("v,verbose", "0: default, 1: print trial", cxxopts::value<int>()->default_value("0"))
    ("i,seed", "random seed for reproducibility", cxxopts::value<int>()->default_value("0"))
    ("a, a_param", "hyperparameter A", cxxopts::value<float>()->default_value("1.8"))
    ("b, b_param", "hyperparameter B", cxxopts::value<float>()->default_value("2.0"));

  auto result = options.parse(argc, argv);

  // command-line arguments
  std::string binaryFile = result["xclbin_path"].as<std::string>();
  std::string problem_matrix_path = result["problem_path"].as<std::string>();
  ParsedVector parsed_result = CNFparser( problem_matrix_path, false );
  float A_float = result["a_param"].as<float>();
  float B_float = result["b_param"].as<float>();
  int NUM_STEPS = result["num_steps"].as<int>();
  int quantum = result["quantum"].as<int>();
  int is_wcnf = result["is_wcnf"].as<int>();
  bool early_stop = result["early_stop"].as<int>();
  int timeout = result["timeout"].as<int>();
  int seed = result["seed"].as<int>();
  int verbose = result["verbose"].as<int>();
  bool DBG_MODE = result["debug"].as<int>();

  int K = parsed_result.k;
  int NUM_VAR = parsed_result.num_var;
  int NUM_CLAUSES = parsed_result.num_clauses;

  // calculate array sizes
  size_t transposed_nnz_size_bytes = sizeof( int ) * NUM_VAR;
  size_t nnz_size_bytes = sizeof( int ) * NUM_CLAUSES;
  size_t svec_size_bytes = sizeof( fixed_point_32t ) * NUM_VAR * num_trials;
  size_t avec_size_bytes = sizeof( fixed_point_32t ) * NUM_CLAUSES * num_trials;

  // allocate arrays
  std::vector<int, aligned_allocator<int>> transposed_nnz( NUM_VAR );
  std::vector<int, aligned_allocator<int>> nnz( NUM_CLAUSES );
  std::vector<int, aligned_allocator<int>> problem_matrix;
  std::vector<int, aligned_allocator<int>> problem_matrix_T;
  std::vector<fixed_point_32t, aligned_allocator<fixed_point_32t>> solution_s( NUM_VAR * num_trials );
  std::vector<fixed_point_32t, aligned_allocator<fixed_point_32t>> solution_a( NUM_CLAUSES * num_trials );
  std::vector<fixed_point_32t, aligned_allocator<fixed_point_32t>> best_solution( NUM_VAR );
  

  std::cout << "problem matrix     : " << problem_matrix_path << "\n";
  std::cout << "number of variables: " << NUM_VAR << "\n";
  std::cout << "number of clauses  : " << NUM_CLAUSES << "\n";
  std::cout << "                  K: " << K << "\n";
  std::cout << "    number of steps: " << NUM_STEPS << "\n";
  std::cout << "                  A: " << A_float << "\n";
  std::cout << "                  B: " << B_float << "\n";
  std::cout << "            quantum: " << quantum << "\n";
  std::cout << "         early stop: " << early_stop << "\n";
  std::cout << "       timeout (us): " << timeout << "\n";
std::cout << "           DBG mode: " << DBG_MODE << "\n";
  std::cout << "timer started \n";

  //--------------------------------------
  // FileIO
  //--------------------------------------

  // read in problem matrix via FileIO
  std::ifstream ProblemFile;
  ProblemFile.open( problem_matrix_path );

  if ( ProblemFile.fail() ) {
    std::cout << "File not found!"
              << "\n";
    return -1;
  }

  std::cout << "File opened successfully!"
            << "\n";

  auto nnz_iterator = nnz.begin();
  for ( std::string line; std::getline( ProblemFile, line ); ) {
    if ( line[0] == 'c' || line[0] == 'p' || line.empty() ) {
      continue;
    }

    std::stringstream iss( line );

    int number;
    int counter = 0;
    bool first_num = true;
    while ( iss >> number ) {
      //skip first number (which is always 1)
      if ( first_num && is_wcnf ) {
        first_num = false;
        continue;
      }

      if ( number == 0 ) {
        break;
      }

      problem_matrix.push_back( number );
      counter++;
    }

    // document nnz
    *nnz_iterator = counter;
    nnz_iterator++;
  }
  ProblemFile.close();

  //--------------------------------------
  // transpose and compress problem matrix
  //--------------------------------------

  // transpose problem matrix
  std::vector<std::vector<int>> problem_matrix_temp( NUM_VAR );
  int problem_matrix_idxcnt = 0;
  for ( int i = 0; i < NUM_CLAUSES; i++ ) {
    for ( int j = 0; j < nnz[i]; j++ ) {
      int value = problem_matrix[problem_matrix_idxcnt++];
      if ( value < 0 ) {
        problem_matrix_temp[-1 * value - 1].push_back( -1 * i - 1 );
      }
      else if ( value > 0 ) {
        problem_matrix_temp[value - 1].push_back( i + 1 );
      }
    }
  }

  // calculate transposed_nnz vector
  for ( int i = 0; i < NUM_VAR; i++ ) {
    transposed_nnz[i] = int( problem_matrix_temp[i].size() );
  }

  // fill in problem_matrix_T
  for ( int i = 0; i < NUM_VAR; i++ ) {
    for ( size_t j = 0; j < problem_matrix_temp[i].size(); j++ ) {
      problem_matrix_T.push_back( problem_matrix_temp[i][j] );
    }
  }

  // sanity checks
  assert( ( problem_matrix_T.size() == problem_matrix.size() ) && "transposed problem matrix has different size then problem matrix! Something is seriously wrong!!!" );
  assert( ( problem_matrix_T.size() < MAX_TOTAL_SIZE ) && "problem matrix size exceeds maximal total size! Consider increasing K_max in compilation." );

  size_t matrix_size_bytes = sizeof( int ) * problem_matrix.size();
  size_t matrix_T_size_bytes = sizeof( int ) * problem_matrix_T.size();
  int total_size = problem_matrix.size();

  // initialize starting point
  for ( int n = 0; n < num_trials; n++ ) {
    seed += n;
    srand( seed );
    for ( int i = 0; i < NUM_VAR; i++ ) {
      solution_s[i * num_trials + n] = 0.5 * ( ( ( (float)rand() ) / (float)RAND_MAX ) * 2 - 1.0 );
    }
  };

  for ( int n = 0; n < num_trials; n++ ) {
    seed += n;
    srand( seed );
    for ( int j = 0; j < NUM_CLAUSES; j++ ) {
      solution_a[j * num_trials + n] = 1 * ( ( ( (float)rand() ) / (float)RAND_MAX ) ) + 3;
    }
  };

  //--------------------------------------
  // FPGA control
  //--------------------------------------
  std::vector<cl::Device> devices = xcl::get_xil_devices();
  cl::Device device = devices[0];
  OCL_CHECK( err,
             std::string device_name = device.getInfo<CL_DEVICE_NAME>( &err ) );
  std::cout << "Found Device=" << device_name.c_str() << std::endl;
  // Creating Context and Command Queue for selected Device
  OCL_CHECK( err, cl::Context context( device, NULL, NULL, NULL, &err ) );
  OCL_CHECK( err, cl::CommandQueue q( context, device, CL_QUEUE_PROFILING_ENABLE,
                                      &err ) );
  std::cout << "Created context and q\n";

  // Read the Binary Files
  auto fileBuf = xcl::read_binary_file(binaryFile);
  cl::Program::Binaries bins{{fileBuf.data(), fileBuf.size()}};
  std::cout << "Found binary\n";

  // Create program
  devices.resize( 1 );
  OCL_CHECK( err, cl::Program program( context, devices, bins, NULL, &err ) );

  // Allocate Buffer in Global Memory
  OCL_CHECK( err, cl::Buffer buffer_problem_matrix(
                      context, CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY,
                      matrix_size_bytes, problem_matrix.data(), &err ) );
  OCL_CHECK( err, cl::Buffer buffer_problem_matrix_T(
                      context, CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY,
                      matrix_T_size_bytes, problem_matrix_T.data(), &err ) );
  OCL_CHECK( err, cl::Buffer buffer_transposed_nnz(
                      context, CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY,
                      transposed_nnz_size_bytes, transposed_nnz.data(), &err ) );
  OCL_CHECK( err, cl::Buffer buffer_nnz(
                      context, CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY,
                      nnz_size_bytes, nnz.data(), &err ) );
  OCL_CHECK( err, cl::Buffer buffer_solution_s(
                      context, CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY,
                      svec_size_bytes, solution_s.data(), &err ) );
  OCL_CHECK( err, cl::Buffer buffer_solution_a(
                      context, CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY,
                      avec_size_bytes, solution_a.data(), &err ) );

  // create the kernel (to execute the k-sat solver)
  OCL_CHECK( err, cl::Kernel krnl_k_sat( program, "ksat", &err ) );

  // set kernel arguments
  OCL_CHECK( err, err = krnl_k_sat.setArg( 0, buffer_problem_matrix ) );
  OCL_CHECK( err, err = krnl_k_sat.setArg( 1, buffer_problem_matrix_T ) );
  OCL_CHECK( err, err = krnl_k_sat.setArg( 2, buffer_nnz ) );
  OCL_CHECK( err, err = krnl_k_sat.setArg( 3, buffer_transposed_nnz ) );
  OCL_CHECK( err, err = krnl_k_sat.setArg( 4, buffer_solution_s ) );
  OCL_CHECK( err, err = krnl_k_sat.setArg( 5, buffer_solution_a ) );
  OCL_CHECK( err, err = krnl_k_sat.setArg( 6, NUM_VAR ) );
  OCL_CHECK( err, err = krnl_k_sat.setArg( 7, NUM_CLAUSES ) );
  OCL_CHECK( err, err = krnl_k_sat.setArg( 8, quantum ) );
  OCL_CHECK( err, err = krnl_k_sat.setArg( 9, A_float ) );
  OCL_CHECK( err, err = krnl_k_sat.setArg( 10, B_float ) );
  OCL_CHECK( err, err = krnl_k_sat.setArg( 11, total_size ) );

  // timer START
  auto start = std::chrono::high_resolution_clock::now();

  // Transfer input data from host to FPGA (over PCIe)
  OCL_CHECK( err,
             err = q.enqueueMigrateMemObjects(
                 { buffer_problem_matrix, buffer_problem_matrix_T, buffer_nnz, buffer_transposed_nnz, buffer_solution_s, buffer_solution_a },
                 0 /* 0 means from host*/ ) );

  //--------------------------------------
  // Incomplete solving loop...
  //--------------------------------------
  // only two way to exit:
  // 1. Found a solution
  // 2. Timed-out

  // use this variable to keep track of the lowest cost
  int single_trial_cost = 0;
  int maxsat_COST_glob = NUM_CLAUSES; // FIXME: need change if doing weighted wcnf problems

  while ( true ) {
    for ( int c_step = quantum; c_step <= NUM_STEPS; c_step += quantum ) {

      // timeout check
      auto currentruntime = std::chrono::duration_cast<std::chrono::microseconds>( std::chrono::high_resolution_clock::now() - start );
      int durationcount = currentruntime.count();
      std::cout << "current duration is " << durationcount << " microseconds. \n";
      if ( durationcount > timeout && timeout != -1 ) {
        std::cout << "Timeout after " << durationcount << " microseconds. No solution found.\n"
                  << std::flush;
        goto SOLUTION_FOUND;
      }

      // Launch the Kernel
      // For HLS kernels global and local size is always (1,1,1). So, it is
      // recommended to always use enqueueTask() for invoking HLS kernel
      OCL_CHECK( err, err = q.enqueueTask( krnl_k_sat ) );

      // Copy Result from Device Global Memory to Host Local Memory
      OCL_CHECK( err, err = q.enqueueMigrateMemObjects(
                          { buffer_solution_s, buffer_solution_a },
                          CL_MIGRATE_MEM_OBJECT_HOST ) );
      q.finish();

      // check solution
      std::cout << "Checking solution for step " << c_step << "... ";

      // print vectors if in verbosed mode
      if (verbose){
        printfixvec(solution_s, "solution_s coming out from kernel");
        printfixvec(solution_a, "solution_a coming out from kernel");
      }

      // perform NaN checking and back tracking if DBG_MODE is enabled
      if ( DBG_MODE ) {
        printfixvec(solution_s, "solution_s");
        printfixvec(solution_a, "solution_a");
      }

      // ------------------------------------
      // SOLUTION CHECK BEGINNNNNNNN!!!!!!!
      // ------------------------------------
      //for (int idk = 0; idk < problem_matrix.size(); idk ++){
      //  std::cout << problem_matrix[idk] << " ";
      //}

      for (int n = 0; n < num_trials; n++){
        //std::cout << "going through trial " << n << "\n";
        single_trial_cost = 0;
        int pm_idx = 0;
        for (int nc = 0; nc < NUM_CLAUSES; nc++){
          //std::cout << "checking the " << nc << " clause.\n"; 
          float res = 0;
          for (int k_schk_idx = 0; k_schk_idx < nnz[nc]; k_schk_idx++){
            int value = problem_matrix[pm_idx++];
            //std::cout << "pm_idx = " << pm_idx << "\n";
            //std::cout << "pm value is " << value << "\n";
            if (value > 0){
              //std::cout << "plus " << solution_s[(value-1) * num_trials + n] << " after sat_a_num this becomes " << sat_a_num( solution_s[(value-1) * num_trials + n]) << "\n";
              res += sat_a_num( solution_s[(value-1) * num_trials + n]);
            }
            else{
              //std::cout << "minus " << solution_s[(-value-1) * num_trials + n] << " after sat_a_num this becomes " << sat_a_num( solution_s[(-value-1) * num_trials + n]) << "\n";
              res -= sat_a_num( solution_s[(-value-1) * num_trials + n]);
            }

            //std::cout << "res is " << res << "\n";
          }

          if (res > -nnz[nc]){
            continue;
          }
          else{
            //std::cout << "none sat clause is number " << nc << "!!!!\n";
            single_trial_cost ++;
          }
        }

        if (single_trial_cost < maxsat_COST_glob){
          maxsat_COST_glob = single_trial_cost;
          // copy the best trial solution over to best_solution vector
          for (int best_solution_idx = 0; best_solution_idx < NUM_VAR; best_solution_idx ++ ){
            best_solution[best_solution_idx] = solution_s[best_solution_idx * num_trials + n ];
          }
        }
      }

      if ( maxsat_COST_glob == 0 ) {
        std::cout << "Solution found! \n\n";

        // FIXME (Owen) - this is not a good programming habit. consider removing this in the future!!!
        goto SOLUTION_FOUND;
      }

      // ------------------------------------
      // SOLUTION CHECK ENDDDDDDDD!!!!!!!
      // ------------------------------------

      std::cout << "Solution not found or early stop disabled.\n";

      std::cout << "o " << maxsat_COST_glob << " timestamp: " << durationcount << " | current o " << single_trial_cost << "\n";
    }

    // no solution found after running NUM_Steps steps
    std::cout << "No solution found after running NUM_STEPS. Reinitialize vectors S and A.\n";

    // initialize starting point
    for ( int n = 0; n < num_trials; n++ ) {
      seed += n;
      srand( seed );
      for ( int i = 0; i < NUM_VAR; i++ ) {
        solution_s[i * num_trials + n] = 0.5 * ( ( ( (float)rand() ) / (float)RAND_MAX ) * 2 - 1.0 );
      }
    };

    for ( int n = 0; n < num_trials; n++ ) {
      seed += n;
      srand( seed );
      for ( int j = 0; j < NUM_CLAUSES; j++ ) {
        solution_a[j * num_trials + n] = 1 * ( ( ( (float)rand() ) / (float)RAND_MAX ) ) + 3;
      }
    };

    if (verbose){
      printfixvec(solution_s, "solution_s after resetting");
      printfixvec(solution_a, "solution_a after resetting");
    }

    OCL_CHECK( err,
               err = q.enqueueMigrateMemObjects(
                   { buffer_solution_s, buffer_solution_a },
                   0 /* 0 means from host*/ ) );
  }

SOLUTION_FOUND:
  // timer stop
  auto stop = std::chrono::high_resolution_clock::now();

  // timer stop
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>( stop - start );
  std::cout << "\n\ntimer stopped \n";

  // print timer information
  std::cout << "Computation time is " << duration.count() << " microseconds"
            << "\n";

  printsol(best_solution, "best solution found");

  std::cout << "\n";

  return 0;  // can edit this to reflect a success criterion
}
