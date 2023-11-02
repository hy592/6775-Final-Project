#if !defined(SOLVER_CONSTANTS_H)
#define SOLVER_CONSTANTS_H 1

#include <ap_fixed.h>

typedef ap_fixed<12, 12, AP_RND, AP_SAT> fixed_point_2t;
typedef ap_fixed<32, 6, AP_RND, AP_SAT> fixed_point_32t;

// compile-time parameters
// tunables
const int MAX_VAR = 512;
const int MAX_CLAUSES = 2048;
const int K_max = 3;
const int num_trials = 192;

const fixed_point_32t time_step = 0.05;

// for kernel estimation
const int MAX_TOTAL_SIZE = K_max * MAX_CLAUSES;
const int NUM_VAR_TRIPCOUNT = MAX_VAR;
const int NUM_CLAUSES_TRIPCOUNT = MAX_CLAUSES;
const int NNZ_MEAN_TRIPCOUNT = MAX_CLAUSES * K_max / MAX_VAR;
const int P_TRIPCOUNT = MAX_CLAUSES * K_max;
const int NUM_ITERATIONS = 500;

#endif