#include "ap_fixed.h"
#include "hls_half.h"
#include "ap_int.h"
#include <hls_stream.h>

typedef ap_uint<32> bit32_t;

#define N 65
#define MAX_WIDTH 20
#define intBits 2
#define numProblems 10

typedef ap_fixed<MAX_WIDTH, intBits> data_type_J;       // weights matrix
typedef ap_fixed<MAX_WIDTH, intBits+1> data_type_x;     // spain vector
typedef ap_fixed<MAX_WIDTH, intBits+4> data_type_e;     // energy vector
typedef half data_t;  // Use data-type half
typedef ap_int<2> spin_sign;

// Maximize the number of MIMO channels N we can decode per FPGA, 
// subject to a constraint on the maximum allowable error rate.

#define sign_bit(x) (x > 0) ? 1 : ((x==0.0 | x==0) ? 0 : -1)

class AHC{
    public:	
        AHC();

        data_type_e bestEnergySpins(spin_sign bestSpins[N]);

        void ahc_solver(data_type_x x_init[N]);

        void matmul();
        void update();
        void setSpins();
        void IsingEnergy();

        void reset();

        void writeDebug(int index);
        void updateX(
            data_type_x x_init[N]
            // ap_fixed<MAX_WIDTH, 2> coupling_strength_new, 
            // ap_fixed<MAX_WIDTH,2> new_mu
        );
        void updateJ(data_type_J J_init[N][N]);

    private:
        ap_fixed<MAX_WIDTH, 3> dt;
        ap_fixed<MAX_WIDTH, 3> r;
        ap_fixed<MAX_WIDTH, 3> beta;
        ap_fixed<MAX_WIDTH, 3> coupling_strength;

        ap_fixed<MAX_WIDTH, 3> mu;
        ap_fixed<MAX_WIDTH, 3> target_a_baseline;
        ap_fixed<MAX_WIDTH, 3> target_a;

        int num_time_steps;

        data_type_J J[N][N];
        data_type_x x[N];
        data_type_x xx[N];

        data_type_x MVM_out[N];
        data_type_e de[N];
        data_type_e e[N];

        data_type_e bestEnergy = 10;
        spin_sign bestSpins[N];
        spin_sign lastSpins[N];
};

// void ahc_top(
//     data_type_J J_matrix[N][N], 
//     data_type_x x_init[N], 
//     spin_sign bestSpinsOut[N]
// );

void dut(hls::stream<bit32_t> &strm_in, hls::stream<bit32_t> &strm_out);
