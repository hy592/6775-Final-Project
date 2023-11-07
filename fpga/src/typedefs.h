#include "ap_fixed.h"
#include "hls_half.h"

#define N 64
#define MAX_WIDTH 16
#define intBits 2

typedef ap_fixed<MAX_WIDTH, intBits> data_type_J;       // weights matrix
typedef ap_fixed<MAX_WIDTH, intBits+1> data_type_x;     // spain vector
typedef ap_fixed<MAX_WIDTH, intBits+4> data_type_e;     // energy vector
typedef half data_t;  // Use data-type �half�
typedef ap_int<2> spin_sign;
