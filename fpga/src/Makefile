#==========================================================================
# Makefile
#==========================================================================
# @brief: A makefile the compiles and synthesizes the digitrec program
#
# @desc: 1. "make" runs csim by default
#        2. "make csim" compiles & executes the fixed-point implementation
#        3. "make clean" cleans up the directory


# Extract Vivado HLS include path
XILINX_VIVADO?=/opt/xilinx/Xilinx_Vivado_vitis_2019.2/Vivado/2019.2
XIL_HLS=source $(XILINX_VIVADO)/settings64.sh; vivado_hls
VHLS_INC=$(XILINX_VIVADO)/include
# Specify compilation flags
CFLAGS=-g -I${VHLS_INC} -DHLS_NO_XIL_FPO_LIB -std=c++11 -D VERSION=low

TCL_SCRIPT=run_base.tcl

.PHONY: all csim bitstream clean

all: csim

digitrec: digitrec.cpp digitrec_test.cpp
	g++ ${CFLAGS} $^ -o $@ -lrt

result/digitrec_csim.txt: digitrec
	@echo "Running AHC sim..."
	mkdir -p result
	./$< | tee $@

csim: result/digitrec_csim.txt
	@echo "Result reorded to $<"

xillydemo.bit:
	@echo "=================================================================" 
	@echo "Synthesizing digitrec and creating bitstream for $(TCL_SCRIPT)..."
	@echo "=================================================================" 
	$(XIL_HLS) -f $(TCL_SCRIPT)
	./run_bitstream.sh

bitstream: xillydemo.bit
	@echo "Bitstream saved to $<"

AHC-fpga: host.cpp
	@echo "Compiling host program"
	g++ ${CFLAGS} $^ -o $@
	@echo "Make sure bitstream is loaded!"

fpga: AHC-fpga
	@echo "Running digitrec accelerator ..."
	./$<

test: AHC_test
	g++ ${CFLAGS} -o AHC_test AHC_test.cpp AHC_low.cpp
	./AHC_test > test_output_low.txt

AHC_test: AHC_test.cpp
	g++ ${CFLAGS} -o AHC_test AHC_test.cpp AHC_low.cpp


clean:
	rm -rf digitrec *.dat *.prj *.log
	rm -rf zedboard_project* xillydemo.bit test_output_low.txt
	rm -rf digitrec-fpga
