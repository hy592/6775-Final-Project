//
//  ksat_kernel.cpp
//  Ising-k-SAT
//

#include <iostream>

#include "solver_compile_consts.hpp"

extern "C" {
extern "C++" {
// piecewise functions in the ODE
template <class T>
inline T func_f( T val )
{
#pragma HLS INLINE
  if ( val >= 1 )
    return T( 1.0 );

  else if ( val <= -1 )
    return T( -1.0 );

  else
    return val;
}

template <class T>
inline T func_g( T val )
{
#pragma HLS INLINE
  if ( val >= 1 )
    return T( 1.0 );

  else if ( val <= 0 )
    return T( 0.0 );

  else
    return val;
}
}

void ksat( const int *problem_matrix,
           const int *problem_matrix_T,
           const int *nnz,
           const int *transposed_nnz,
           fixed_point_32t *solution_s,
           fixed_point_32t *solution_a,
           const int NUM_VAR,
           const int NUM_CLAUSES,
           const int NUM_STEPS,
           float A_float,
           float B_float,
           const int total_size )
{
  // Vitis 2020 will infer AXI interface type
  // allocate local memory space for problem matrix, variables, injection terms in SRAM

  fixed_point_2t local_problem_matrix[MAX_CLAUSES * K_max];
  fixed_point_2t local_problem_matrix_T[MAX_CLAUSES * K_max];
  fixed_point_2t local_nnz[MAX_CLAUSES];
  fixed_point_2t local_transposed_nnz[MAX_VAR];
  fixed_point_32t local_solution_s[MAX_VAR][num_trials];
  fixed_point_32t local_solution_a[MAX_CLAUSES][num_trials];
  fixed_point_32t temp_solution_s[MAX_VAR][num_trials];
  fixed_point_32t temp_solution_a[MAX_CLAUSES][num_trials];
  fixed_point_32t mapped_funcf_s[MAX_VAR][num_trials];
  fixed_point_32t mapped_funcg_a[MAX_CLAUSES][num_trials];
  fixed_point_32t inject_s[MAX_CLAUSES][num_trials];
  fixed_point_32t inject_a[MAX_VAR][num_trials];
  fixed_point_32t A = A_float;
  fixed_point_32t B = B_float;

// bind-storage
// #pragma HLS bind_storage variable=temp_solution_s type=RAM_1WNR impl=LUTRAM
// #pragma HLS bind_storage variable=temp_solution_a type=RAM_1WNR impl=LUTRAM
#pragma HLS bind_storage variable=mapped_funcf_s type=RAM_1WNR impl=URAM
#pragma HLS bind_storage variable=mapped_funcg_a type=RAM_1WNR impl=URAM
#pragma HLS bind_storage variable=inject_s type=RAM_2P impl=URAM
#pragma HLS bind_storage variable=inject_a type=RAM_2P impl=URAM

// define array partition factor
#pragma HLS ARRAY_PARTITION variable = local_solution_a complete dim = 2
#pragma HLS ARRAY_PARTITION variable = local_solution_s complete dim = 2
#pragma HLS ARRAY_PARTITION variable = temp_solution_s complete dim = 2
#pragma HLS ARRAY_PARTITION variable = temp_solution_a complete dim = 2
#pragma HLS ARRAY_PARTITION variable = mapped_funcf_s complete dim = 2
#pragma HLS ARRAY_PARTITION variable = mapped_funcg_a complete dim = 2
#pragma HLS ARRAY_PARTITION variable = inject_s complete dim = 2
#pragma HLS ARRAY_PARTITION variable = inject_a complete dim = 2

  // --------------------------------------------------------------------
  // read data from DDR to local memories
  // store the problem matrix in 2 fixed point bit number
  // --------------------------------------------------------------------
  // read problem matrix

mem_read_pm:
  for ( int i = 0; i < total_size; i++ ) {
#pragma HLS loop_tripcount min = P_TRIPCOUNT max = P_TRIPCOUNT avg = P_TRIPCOUNT
    local_problem_matrix[i] = fixed_point_2t( problem_matrix[i] );
  }

mem_read_ptm:
  for ( int i = 0; i < total_size; i++ ) {
#pragma HLS loop_tripcount min = P_TRIPCOUNT max = P_TRIPCOUNT avg = P_TRIPCOUNT
    local_problem_matrix_T[i] = fixed_point_2t( problem_matrix_T[i] );
  }

mem_read_transposednnz:
  for ( int i = 0; i < NUM_VAR; i++ ) {
#pragma HLS loop_tripcount min = NUM_VAR_TRIPCOUNT max = NUM_VAR_TRIPCOUNT avg = NUM_VAR_TRIPCOUNT
    local_transposed_nnz[i] = fixed_point_2t( transposed_nnz[i] );
  }

mem_read_nnz:
  for ( int i = 0; i < NUM_CLAUSES; i++ ) {
#pragma HLS loop_tripcount min = NUM_CLAUSES_TRIPCOUNT max = NUM_CLAUSES_TRIPCOUNT avg = NUM_CLAUSES_TRIPCOUNT
    local_nnz[i] = fixed_point_2t( nnz[i] );
  }

// read solution variables
mem_read_s:
  for ( int i = 0; i < NUM_VAR; i++ ) {
#pragma HLS loop_tripcount min = NUM_VAR_TRIPCOUNT max = NUM_VAR_TRIPCOUNT avg = NUM_VAR_TRIPCOUNT
  mem_read_s_inner:
    for ( int j = 0; j < num_trials; j++ ) {
      local_solution_s[i][j] = solution_s[i * num_trials + j];
    }
  }

// read auxiliary variables
mem_read_a:
  for ( int i = 0; i < NUM_CLAUSES; i++ ) {
#pragma HLS loop_tripcount min = NUM_CLAUSES_TRIPCOUNT max = NUM_CLAUSES_TRIPCOUNT avg = NUM_CLAUSES_TRIPCOUNT
  mem_read_a_inner:
    for ( int j = 0; j < num_trials; j++ ) {
      local_solution_a[i][j] = solution_a[i * num_trials + j];
    }
  }

// --------------------------------------------------------------------
// evolution of the system
// --------------------------------------------------------------------
step_iteration:
  for ( int t = 0; t < NUM_STEPS; t++ ) {
#pragma HLS loop_tripcount min = NUM_ITERATIONS max = NUM_ITERATIONS avg = NUM_ITERATIONS

    // --------------------------------------------------------------------
    // initialization related with solution, auxiliary variables
    // injection term, array copy for variables, array of variables mapped to
    // functions
    // --------------------------------------------------------------------

  s_var_init:
    for ( int i = 0; i < NUM_VAR; i++ ) {
#pragma HLS loop_tripcount min = NUM_VAR_TRIPCOUNT max = NUM_VAR_TRIPCOUNT avg = NUM_VAR_TRIPCOUNT
      for ( int j = 0; j < num_trials; j++ ) {
#pragma HLS unroll
        inject_a[i][j] = fixed_point_32t( 0.0 );
        temp_solution_s[i][j] = local_solution_s[i][j];
        mapped_funcf_s[i][j] = func_f( local_solution_s[i][j] );

        inject_s[i][j] = fixed_point_32t( 0.0 );
        temp_solution_a[i][j] = local_solution_a[i][j];
        mapped_funcg_a[i][j] = func_g( local_solution_a[i][j] );
      }
    }

  a_var_init:
    for ( int i = NUM_VAR; i < NUM_CLAUSES; i++ ) {
#pragma HLS loop_tripcount min = NUM_CLAUSES_TRIPCOUNT max = NUM_CLAUSES_TRIPCOUNT avg = NUM_CLAUSES_TRIPCOUNT
      for ( int j = 0; j < num_trials; j++ ) {
#pragma HLS unroll
        inject_s[i][j] = fixed_point_32t( 0.0 );
        temp_solution_a[i][j] = local_solution_a[i][j];
        mapped_funcg_a[i][j] = func_g( local_solution_a[i][j] );
      }
    }

    // --------------------------------------------------------------------
    // matrix-vector product
    // compute injection term in the dynamics of solution, auxiliary variable,
    // matrix vector multiplication interchange the order of two for loops
    // --------------------------------------------------------------------

    int scounter = 0;
  mat_vec_souter:
    for ( int i = 0; i < NUM_CLAUSES; i++ ) {
#pragma HLS loop_tripcount min = NUM_CLAUSES_TRIPCOUNT max = NUM_CLAUSES_TRIPCOUNT avg = NUM_CLAUSES_TRIPCOUNT

    mat_vec_sinner:
      for ( int j = 0; j < local_nnz[i]; j++ ) {
#pragma HLS loop_tripcount min = NNZ_MEAN_TRIPCOUNT max = NNZ_MEAN_TRIPCOUNT avg = NNZ_MEAN_TRIPCOUNT
#pragma HLS PIPELINE

        int value = int( local_problem_matrix[scounter] );
        scounter++;

      mat_vec_strial:
        for ( int trial = 0; trial < num_trials; trial++ ) {
          if ( value > 0 ) {
            inject_s[i][trial] += mapped_funcf_s[value - 1][trial];
          }
          else {
            inject_s[i][trial] -= mapped_funcf_s[-value - 1][trial];
          }
        }
      }
    }

    int counter = 0;
  mat_vec_aouter:
    for ( int i = 0; i < NUM_VAR; i++ ) {
#pragma HLS loop_tripcount min = NUM_VAR_TRIPCOUNT max = NUM_VAR_TRIPCOUNT avg = NUM_VAR_TRIPCOUNT

    mat_vec_ainner:
      for ( int j = 0; j < transposed_nnz[i]; j++ ) {
#pragma HLS loop_tripcount min = NNZ_MEAN_TRIPCOUNT max = NNZ_MEAN_TRIPCOUNT avg = NNZ_MEAN_TRIPCOUNT
#pragma HLS PIPELINE

        int value = int( local_problem_matrix_T[counter] );
        counter++;

      mat_vec_atrial:
        for ( int trial = 0; trial < num_trials; trial++ ) {
          if ( value > 0 ) {
            inject_a[i][trial] += mapped_funcg_a[value - 1][trial];
          }
          else {
            inject_a[i][trial] -= mapped_funcg_a[-value - 1][trial];
          }
        }
      }
    }

    // --------------------------------------------------------------------
    // update solution vectors
    // --------------------------------------------------------------------

  s_loop_update:
    for ( int i = 0; i < NUM_VAR; i++ ) {
#pragma HLS loop_tripcount min = NUM_VAR_TRIPCOUNT max = NUM_VAR_TRIPCOUNT avg = NUM_VAR_TRIPCOUNT
      int k = local_nnz[i];
      for ( int j = 0; j < num_trials; j++ ) {
#pragma HLS unroll
        local_solution_s[i][j] = temp_solution_s[i][j] + time_step * ( -temp_solution_s[i][j] + A * mapped_funcf_s[i][j] + inject_a[i][j] );
        local_solution_a[i][j] = temp_solution_a[i][j] + time_step * ( -temp_solution_a[i][j] + B * mapped_funcg_a[i][j] - inject_s[i][j] + fixed_point_32t( 1 - k ) );
      }
    }

  a_loop_update:
    for ( int i = NUM_VAR; i < NUM_CLAUSES; i++ ) {
#pragma HLS loop_tripcount min = NUM_CLAUSES_TRIPCOUNT max = NUM_CLAUSES_TRIPCOUNT avg = NUM_CLAUSES_TRIPCOUNT
      int k = local_nnz[i];
      for ( int j = 0; j < num_trials; j++ ) {
#pragma HLS unroll
        local_solution_a[i][j] = temp_solution_a[i][j] + time_step * ( -temp_solution_a[i][j] + B * mapped_funcg_a[i][j] - inject_s[i][j] + fixed_point_32t( 1 - k ) );
      }
    }
  }

// --------------------------------------------------------------------
// write data from local memories to the DDR
// --------------------------------------------------------------------
mem_write_s:
  for ( int i = 0; i < NUM_VAR; i++ ) {
#pragma HLS loop_tripcount min = NUM_VAR_TRIPCOUNT max = NUM_VAR_TRIPCOUNT avg = NUM_VAR_TRIPCOUNT
  mem_write_s_inner:
    for ( int j = 0; j < num_trials; j++ ) {
      solution_s[i * num_trials + j] = local_solution_s[i][j];
    }
  }

// read auxiliary variables
mem_write_a:
  for ( int i = 0; i < NUM_CLAUSES; i++ ) {
#pragma HLS loop_tripcount min = NUM_CLAUSES_TRIPCOUNT max = NUM_CLAUSES_TRIPCOUNT avg = NUM_CLAUSES_TRIPCOUNT
  mem_write_a_inner:
    for ( int j = 0; j < num_trials; j++ ) {
      solution_a[i * num_trials + j] = local_solution_a[i][j];
    }
  }
}
}
