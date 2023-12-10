#==========================================================================
# Makefile
#==========================================================================
# @brief: A makefile the compiles and runs the digitrec program
#
# @desc: 1. "make" or "make sw" runs software execution by default
#        2. "make fpga" invokes the HW accelerator
#        3. "make clean" cleans up the directory


INC_PATH=/usr/include/vivado_hls
CFLAGS = -I${INC_PATH}

.PHONY: all sw fpga

all: sw

digitrec-arm: digitrec.cpp digitrec_test.cpp
	@echo "Compiling & executing digitrec software program on ARM ..."
	g++ ${CFLAGS} $^ -o $@

result/digitrec_arm_sim.txt: digitrec-arm
	mkdir -p result
	./$< | tee $@

sw: result/digitrec_arm_sim.txt
	@echo "Result saved to $@"

digitrec-fpga: host.cpp
	@echo "Compiling host program"
	g++ ${CFLAGS} $^ -o $@
	@echo "Make sure bitstream is loaded!"

fpga: digitrec-fpga
	@echo "Running digitrec accelerator ..."
	./$<

# Clean up the directory before submission
clean:
	@echo "Clean up output files"
	rm -rf digitrec-sw vivado_hls.log *.prj result out.dat *~
	rm -rf digitrec-fpga
