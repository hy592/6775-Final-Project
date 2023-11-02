//#######################################
// Header file for AHC algorithm implementation in HLS
// For the MIMO decoding project
//#######################################

#ifndef __MIMO_AHC_H__
#define __MIMO_AHC_H__

// #define __USE_FLOAT__
#define __DEBUG__

#include <stdio.h>
#include <ap_fixed.h>

#define N 4
#define sign_bit(x) (x>=0)?1:-1

#ifndef __USE_FLOAT__
typedef ap_fixed<10, 7> data_type_J;
typedef ap_fixed<18, 4> data_type_x;
typedef ap_fixed<16, 10> data_type_e;
typedef ap_fixed<16, 2> data_type;
#endif


#ifdef __USE_FLOAT__
typedef float data_type_J;
typedef float data_type_x;
typedef float data_type_e;
typedef float data_type;
#endif

using namespace std;

// Random generator seed
static int rnd_seed = 7;

static int random_integer()
{
	unsigned int hi,lo;

	hi = 16807 * (rnd_seed >> 16);
	lo = 16807 * (rnd_seed & 0xFFFF);
	lo += (hi & 0x7FFF) << 16;
	lo += hi >> 15;
	if (lo > 2147483647)
		lo -= 2147483647;
	rnd_seed = lo;
	return rnd_seed;
}

class AHC{
    public:
    
    data_type dt;
    data_type r;  
    data_type beta;
    data_type coupling_strength;

#ifndef __USE_FLOAT__
    ap_fixed<6, 3> mu;
    ap_fixed<6, 2> target_a_baseline, target_a;
#endif

#ifdef __USE_FLOAT__
   float mu;
   float target_a_baseline, target_a;
#endif

    data_type_J J[N][N];
    data_type_x x[N];
    data_type_x MVM[N];
    data_type_e e[N];

    data_type_x xx[N], dx[N], dx2[N];

	data_type_e de[N];

    int num_time_steps;

    AHC();
    void ahc_solver();
    void square();
    void matmul();
    void update();
    void reset();

#ifdef __DEBUG__
    void debug(int);
#endif

};

extern "C"{
void ahc(float J_matrix[N*N], float out_vector[N], float value, short option);
}
#endif
