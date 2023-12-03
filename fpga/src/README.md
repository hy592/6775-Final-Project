# Final Project README

## Run C Simulation

## Run C Synthesis
```bash
source /classes/ece6775/setup-ece6775.sh # setup env
vitis_hls -f run_base.tcl
```
<!-- You will implement and evaluate the performance for three designs:

- A baseline digitrec design that does not use any HLS optimization directives (`vivado_hls -f run_base.tcl`).
- An unrolled digitrec design which is similar to what you did in Lab 2 where unrolling and array partitioning are applied (`vivado_hls -f run_unroll.tcl`).
- A pipelined digitrec design which applies loop pipelining in addition to the previous optimizations (`vivado_hls -f run_pipeline.tcl`). -->

## Coding and Debugging
- ecelinux

To run the software on ecelinux, first go to the `lab3/ecelinux` directory, then run `make`

```bash
# in ecelinux
unzip lab3.zip # unzip the archive
cd lab3/ecelinux
source /classes/ece6775/setup-ece6775.sh # setup env
make # builds and runs the csim
```
- Zedboard

To run the software on a ZedBoard (i.e., ARM CPU), first copy lab3 files to a ZedBoard, 

```bash
# in ecelinux
zip -r lab3.zip lab3  # zip the file. -x to ignore unwanted file
scp lab3.zip <user>@zhang-zedboard-xx.ece.cornell.edu:~
ssh <user>@zhang-zedboard-xx.ece.cornell.edu # log in to a ZedBoard
```

then go to the directory `lab3/zedboard`, finally run `make sw` .
```bash
# in zedboard
unzip lab3.zip # unzip the archive
cd lab3/zedboard

make sw # builds and runs software
```

## Generating the Bitstream

1. The first step is to generate verilog from HLS
```bash
# In ecelinux
cd lab3/ecelinux
source /classes/ece6775/setup-ece6775.sh # setup env
make                  # builds and runs the csim

# generate Verilog
vivado_hls -f run_base.tcl
vivado_hls -f run_unroll.tcl
vivado_hls -f run_pipeline.tcl
```

You can examine the directory cordic.prj/solution1/syn/verilog to check that the Verilog files have been generated.

2. The next step is to implement the hardware circuit described by the Verilog on an the FPGA.

```bash
source run_bitstream.sh # generate FPGA bitstream
```

## Programming the FPGA

Before you attempt to login, please first
check the occupancy of the available boards using the following link:
https://www.csl.cornell.edu/courses/ece6775/zedboard.html

```bash
zhang-zedboard-xx.ece.cornell.edu
```
where xx can be 01, 02, ..., 11.

All students have an account on each Zedboard with
their NetID as the username and their case-sensitive last name as their initial password.

```bash
# In ecelinux,
# Copy the bitstream from ecelinux to Zedboard
scp xillydemo.bit hy592@zhang-zedboard-09.ece.cornell.edu:~

# Login to the Zedboard
ssh hy592@zhang-zedboard-09.ece.cornell.edu

# In Zedboard
# Mount SD card. Try rebooting if there is a device busy error
mount /mnt/sd

# Copy the bitstream file to the SD card and reboot the Zedboard
cp xillydemo.bit /mnt/sd
sudo reboot
```
The Zedboard will restart. Wait about 30 seconds and login again — the FPGA will be fully programmed after the reboot

Then, to execute the accelerator on the FPGA, you need go to directory `../zedboard` and run `make fpga`.

```bash
# Login to the Zedboard
ssh hy592@zhang-zedboard-09.ece.cornell.edu
# run FPGA program
cd lab3/zedboard
make fpga
```

```bash
# update host.cpp
scp host.bit hy592@zhang-zedboard-09.ece.cornell.edu:~
```

## 

## dut interface
```cpp
#include <hls_stream.h>

//-----------------------------------
// dut function (top module)
//-----------------------------------
// @param[in]  : strm_in - input stream
// @param[out] : strm_out - output stream 
void dut (
    hls::stream<bit32_t> &strm_in,
    hls::stream<bit32_t> &strm_out
)
{
  // define variables
  theta_type theta;
  cos_sin_type c, s;

  // ------------------------------------------------------
  // Input processing
  // ------------------------------------------------------
  // Read the two input 32-bit words (low word first)
  bit32_t input_lo = strm_in.read();
  bit32_t input_hi = strm_in.read();

  // Convert input raw bits to fixed-point representation via bit slicing
  theta(31, 0) = input_lo;
  theta(theta.length()-1, 32) = input_hi;

  // ------------------------------------------------------
  // Call CORDIC 
  // ------------------------------------------------------
  cordic( theta, s, c );

  // ------------------------------------------------------
  // Output processing
  // ------------------------------------------------------
  // Write out the cos value (low word first)
  strm_out.write( c(31, 0) );
  strm_out.write( c(c.length()-1, 32) );
  
  // Write out the sin value (low word first)
  strm_out.write( s(31, 0) );
  strm_out.write( s(s.length()-1, 32) );
}

```

##

## Testbench

To run the testbench, run ```make test``` in ```fpga/src/```. Currently, we have two implementations: ```AHC_low``` and ```AHC_new```. Modify 
- ```make test``` in ```fpga/src/Makefile```
- ```#include```s in ```fpga/src/AHC_test.cpp```

to change versions.

### Makefile

```cmake
test: AHC_test
	g++ ${CFLAGS} -o AHC_test AHC_test.cpp AHC_<low OR new>.cpp
	./AHC_test > test_output_<low OR new>.txt

AHC_test: AHC_test.cpp
	g++ ${CFLAGS} -o AHC_test AHC_test.cpp AHC_<low OR new>.cpp  
```