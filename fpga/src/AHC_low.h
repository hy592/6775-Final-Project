#include "ap_fixed.h"
#include "hls_half.h"

#define N 65
#define MAX_WIDTH 16
#define intBits 2

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
	ap_fixed<MAX_WIDTH, 2> dt;
	ap_fixed<MAX_WIDTH, 2> r;
	ap_fixed<MAX_WIDTH, 2> beta;
	ap_fixed<MAX_WIDTH, 2> coupling_strength;

    ap_fixed<MAX_WIDTH, 3> mu;
    ap_fixed<MAX_WIDTH, 2> target_a_baseline;
    ap_fixed<MAX_WIDTH, 2> target_a;

    int num_time_steps;

    data_type_J J[N][N];
    data_type_x x[N];
    data_type_x MVM_out[N];
    data_type_e e[N];

    data_type_x xx[N], dx[N], dx2[N];

	data_type_e de[N];

    /*
    float dt;
    float r;
    float beta;
    float coupling_strength;

    float mu;
    float target_a_baseline;
    float  target_a;
    //int num_time_steps;

    float J[N][N];
    float x[N];
    float MVM_out[N];
    float e[N];

    float xx[N], dx[N], dx2[N];

    float de[N];

    //float bestEnergy = 100;
    //float bestSpins[N];
    */

	data_type_e bestEnergy = 10;
	spin_sign bestSpins[N];
	spin_sign lastSpins[N];

    AHC(data_type_x x_init[N], data_type_J J_init[N][N]);

    void ahc_solver();
    void square();
    void matmul();
    void update();
    void setSpins();
    void reset();

    void writeDebug(int index);
    void updateSpins(
        data_type_x x_init[N], 
        ap_fixed<MAX_WIDTH, 2> coupling_strength_new, 
        ap_fixed<MAX_WIDTH,2> new_mu
    );

    data_type_e IsingEnergy();
};

void ahc_top(
    data_type_J J_matrix[N][N], 
    data_type_x x_init[N], 
    spin_sign bestSpinsOut[N]
    );
