n_int_bits = 2
n_frac_bits = 16

import math


snr = 30


#N_t x N_r MIMO

N_t = 16                #Number Of Transmit antennas
N_r = 16             #Number of Receive antennas
modulation = 16        #Type of modulation
numInstances = 200    #Number of Ising Instances to generate

bits_per_symbol = int(math.log2(modulation))

bitsPerEntry = int(0.5*bits_per_symbol)

# maxQI = 2*np.sqrt(modulation/4) - 1

# _qam = comm.QAMModem(modulation)

totalBits = bits_per_symbol*N_t

snr_list = [1,5,10,15,20,25,30] ## in dB

# H = np.zeros((N_r,N_t),dtype=np.complex128)

# ## Seed for verification

instance_type = "Nt"+str(N_t)+"_Nr"+str(N_r)+"_"+str(modulation)+"QAM/"

path_to_instances = path_to_all_instances + instance_type + str(snr) + "/"

# Create a fixed-point array
#initialize empty 2D array

k = 0
for k in range(numInstances):
    file = open(path_to_instances + "DI_MIMO_J_Quad_" + str(k) + ".txt", "r")
    J = np.loadtxt(file)
    file.close()
    N = len(J)
    print(J)


    J_quad_fixed = [[] for _ in range(len(J))]


    for i in range(len(J)):
        for j in J[i]:

            val = (int)(j*(2**n_frac_bits))
            # Handle negative values
            if val < 0:
                val = (val + (1 << (n_int_bits+n_frac_bits))) % (1 << (n_int_bits+n_frac_bits))
            hex_val = hex(val)

            #3res = int(hex_val,16) # convert hex to float
            #if(i==1):
            #   print((res))
            
        #  hex_val = hex_val.zfill(4)

            J_quad_fixed[i].append(hex_val)
    #print(J_quad)
    #print(J_quad_fixed)
    #print(J_quad_fixed_back)

    file = open(path_to_instances + "DI_MIMO_J_Conv_6_19_12_2_" + str(k) + ".txt", "w")
    for row in J_quad_fixed:
        file.write(' '.join(row) + '\n')
    file.close()
#=============================================================================
## @brief: A Tcl script for synthesizing the base project
#
## Project name
#set hls_prj Ising_Machine.prj
#
## Open/reset the project
#open_project ${hls_prj} -reset
#
## Top function of the design
#set_top ahc_top
#
## Add design and testbench files
#add_files AHC_low.cpp -cflags "-std=c++11"
#
## Add testbench files
#add_files -tb AHC_test_simple.cpp -cflags "-std=c++11"
#
## open_solution "solution1" -flow_target vivado
#open_solution "solution1" 
#
## Use UltraScale device
## set_part {xcvu9p_CIV-fsgd2104-3-e}
#
## Use Zynq device
#set_part {xc7z020clg484-1}
#
## Target clock period is 10ns
#create_clock -period 10 -name default
#
#### You can insert your own directives here ###
## source "./directives.tcl"
#
## Run C simulation
#csim_design -clean -O
#
## Run synthesis
#csynth_design
#
## Run RTL cosimulation
## cosim_design
#
## Run RTL simulation
## export_design -format ip_catalog
#
#exit
