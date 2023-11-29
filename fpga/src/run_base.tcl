#=============================================================================
# run_base.tcl 
#=============================================================================
# @brief: A Tcl script for synthesizing the base project

# Project name
set hls_prj Ising_Machine.prj

# Open/reset the project
open_project ${hls_prj} -reset

# Top function of the design
set_top ahc_top

# Add design and testbench files
# add_files AHC_low.cpp -cflags "-std=c++11"
add_files AHC_new.cpp -cflags "-std=c++11"

# Add testbench files
add_files -tb AHC_test_simple.cpp -cflags "-std=c++11"

# open_solution "solution1" -flow_target vivado
open_solution "solution1" 

# Use UltraScale device
# set_part {xcvu9p_CIV-fsgd2104-3-e}

# Use Zynq device
set_part {xc7z020clg484-1}

# Target clock period is 10ns
create_clock -period 10 -name default

### You can insert your own directives here ###
# source "./directives.tcl"

# Run C simulation
csim_design -clean -O

# Run synthesis
csynth_design

# Run RTL cosimulation
# cosim_design

# Run RTL simulation
# export_design -format ip_catalog

exit